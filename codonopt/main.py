import os
import csv
import math
import tempfile
import hashlib
from typing import Dict, Any, Iterable, List, Optional, Tuple

from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from codonopt.cli import get_cli_args
from codonopt.logging_utils import setup_logger
from codonopt.exceptions import CodonoptError, InputFormatError, ConstraintError
from codonopt.utils import (
    is_protein_sequence,
    calculate_gc,
    max_homopolymer_length,
    verify_backtranslation_matches_protein,
)
from codonopt.core.optimizer import optimize_sequence
from codonopt.io.codon_tables import load_codon_table_from_xlsx


# ---------------------------
# Batch table helpers
# ---------------------------

def iter_batch_table(path: str) -> Iterable[Dict[str, Any]]:
    delimiter = "\t" if path.lower().endswith(".tsv") else ","
    with open(path, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            raise InputFormatError(f"Batch table has no header row: {path}")

        normalized = [h.strip() if h else h for h in reader.fieldnames]
        reader.fieldnames = normalized

        required = {"sequence", "codon_table"}
        missing = required - set(normalized)
        if missing:
            raise InputFormatError(
                f"Batch table is missing required column(s): {', '.join(sorted(missing))}\n"
                f"Found columns: {normalized}\n"
                f"Fix the header to include at least: sequence,codon_table"
            )

        for row in reader:
            clean = {}
            for k, v in row.items():
                if k is None:
                    continue
                kk = k.strip()
                vv = v.strip() if isinstance(v, str) else v
                clean[kk] = vv
            yield clean


def _parse_list_field(v: Any) -> List[str]:
    if v is None:
        return []
    if isinstance(v, list):
        return v
    s = str(v).strip()
    if not s or s.lower() == "nan":
        return []
    return [x.strip() for x in s.split("|") if x.strip()]


def _parse_float_or_none(v: Any) -> Optional[float]:
    if v is None:
        return None
    if isinstance(v, (float, int)):
        return float(v)
    s = str(v).strip()
    if not s or s.lower() == "nan":
        return None
    return float(s)


def _parse_float_or_default(v: Any, default: float) -> float:
    if v is None:
        return float(default)
    if isinstance(v, (float, int)):
        return float(v)
    s = str(v).strip()
    if not s or s.lower() == "nan":
        return float(default)
    return float(s)


def _parse_int_or_default(v: Any, default: int) -> int:
    if v is None:
        return default
    if isinstance(v, int):
        return v
    s = str(v).strip()
    if not s or s.lower() == "nan":
        return default
    return int(s)


def _parse_boolish(v: Any, default: bool = False) -> bool:
    if v is None:
        return default
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if not s or s in ("nan", "none"):
        return default
    if s in ("1", "true", "t", "yes", "y", "on"):
        return True
    if s in ("0", "false", "f", "no", "n", "off"):
        return False
    return default


def _parse_mode(v: Any) -> str:
    if v is None:
        return "kleinbub"
    s = str(v).strip().lower()
    if not s or s in ("nan", "none"):
        return "kleinbub"
    if s not in ("kleinbub", "strict"):
        raise InputFormatError("optimization_mode must be 'kleinbub' or 'strict'")
    return s


def _normalize_sheet_name(v: Any) -> Optional[str]:
    if v is None:
        return None
    s = str(v).strip()
    if not s or s.lower() in ("nan", "none"):
        return None
    return s


# ---------------------------
# CryptKeeper integration
# ---------------------------

def _extract_site_start(site) -> Optional[int]:
    for attr in ("start", "start_pos", "start_position", "begin", "pos"):
        if hasattr(site, attr):
            try:
                return int(getattr(site, attr))
            except Exception:
                pass
    try:
        for x in site:
            if isinstance(x, int):
                return int(x)
    except Exception:
        pass
    return None


def _extract_site_score(site) -> Optional[float]:
    for attr in ("rbs_score", "score", "tir", "value"):
        if hasattr(site, attr):
            try:
                return float(getattr(site, attr))
            except Exception:
                pass
    try:
        for x in site:
            if isinstance(x, float):
                return float(x)
    except Exception:
        pass
    return None


def run_cryptkeeper_screen(
    dna_seq: str,
    rbs_score_cutoff: float,
    threads: int,
    fail_score: float,
    ignore_first_nt: int,
    ignore_last_nt: int,
    logger,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Pass/fail logic:
      FAIL if any internal predicted start site has score >= fail_score.

    "Internal" means:
      start >= ignore_first_nt
      and start < len(dna_seq) - ignore_last_nt
    """
    info = {
        "cryptkeeper_internal_site_count": 0,
        "cryptkeeper_internal_site_positions": "",
        "cryptkeeper_max_internal_score": "",
    }

    try:
        from cryptkeeper import cryptkeeper as ck
    except Exception as e:
        raise InputFormatError(
            "CryptKeeper is enabled but could not be imported. "
            "Install cryptkeeper (recommended via conda/bioconda). "
            f"Import error: {e}"
        )

    seq_len = len(dna_seq)
    start_min = max(0, int(ignore_first_nt))
    start_max_exclusive = max(0, seq_len - int(ignore_last_nt))

    with tempfile.TemporaryDirectory(prefix="cryptkeeper_") as tmpdir:
        fasta_path = os.path.join(tmpdir, "candidate.fasta")
        with open(fasta_path, "w") as f:
            f.write(">candidate\n")
            f.write(dna_seq.strip().upper() + "\n")

        result = ck(
            input_file=fasta_path,
            output=tmpdir,
            circular=False,
            name="candidate",
            threads=int(threads),
            logger=logger,
            rbs_score_cutoff=float(rbs_score_cutoff),
        )

        sites = getattr(result, "translation_sites", []) or []

        internal_positions: List[int] = []
        internal_scores: List[float] = []

        for site in sites:
            start = _extract_site_start(site)
            if start is None:
                continue

            if start < start_min:
                continue
            if start >= start_max_exclusive:
                continue

            internal_positions.append(start)

            sc = _extract_site_score(site)
            if sc is None:
                continue
            try:
                scf = float(sc)
            except Exception:
                continue
            if not math.isfinite(scf):
                continue
            internal_scores.append(scf)

        info["cryptkeeper_internal_site_count"] = len(internal_positions)
        info["cryptkeeper_internal_site_positions"] = ",".join(str(p) for p in internal_positions[:200])

        max_score = max(internal_scores) if internal_scores else None
        if max_score is not None:
            info["cryptkeeper_max_internal_score"] = max_score

        if max_score is None:
            # No finite scored internal sites => don't fail
            return True, info

        passes = float(max_score) < float(fail_score)
        return passes, info


# ---------------------------
# Job preparation
# ---------------------------

def merge_defaults(args, row: Dict[str, Any], logger) -> Dict[str, Any]:
    from codonopt.defaults import DEFAULT_FORBIDDEN_MOTIFS, ALWAYS_FORBIDDEN_MOTIFS

    job = vars(args).copy()

    for k, v in row.items():
        if v in ("", None):
            continue
        job[k] = v

    # ---------------------------
    # Parse motif & codon fields
    # ---------------------------

    job["avoid_codons"] = _parse_list_field(job.get("avoid_codons"))
    job["avoid_motifs"] = _parse_list_field(job.get("avoid_motifs"))

    def _normalize_motifs(motifs: List[str]) -> List[str]:
        """
        Normalize motifs:
        - uppercase
        - convert RNA U -> DNA T
        - strip whitespace
        """
        out = []
        for m in motifs:
            if m is None:
                continue
            m2 = str(m).strip().upper().replace("U", "T")
            if m2:
                out.append(m2)
        return out

    # Normalize user-provided motifs
    user_motifs = _normalize_motifs(job["avoid_motifs"])

    # Normalize defaults
    default_motifs = _normalize_motifs(DEFAULT_FORBIDDEN_MOTIFS)
    always_motifs = _normalize_motifs(ALWAYS_FORBIDDEN_MOTIFS)

    # Merge: user + defaults + always-forbidden
    job["avoid_motifs"] = sorted(set(user_motifs) | set(default_motifs) | set(always_motifs))

    # ---------------------------
    # Numeric & mode parsing
    # ---------------------------

    job["gc_min"] = _parse_float_or_none(job.get("gc_min"))
    job["gc_max"] = _parse_float_or_none(job.get("gc_max"))
    job["n"] = _parse_int_or_default(job.get("n"), 1)
    job["max_homopolymer"] = _parse_int_or_default(job.get("max_homopolymer"), 5)

    seed_val = job.get("seed")
    if seed_val in ("", None):
        job["seed"] = None
    else:
        try:
            job["seed"] = int(seed_val)
        except Exception:
            raise InputFormatError(f"seed must be an integer (or blank); got: {seed_val!r}")

    job["min_codon_fraction"] = _parse_float_or_default(job.get("min_codon_fraction"), 0.05)

    job["optimization_mode"] = _parse_mode(job.get("optimization_mode"))
    job["include_stop_codon"] = _parse_boolish(job.get("include_stop_codon"), False)

    job["cryptkeeper_enable"] = _parse_boolish(job.get("cryptkeeper_enable"), False)
    job["cryptkeeper_threads"] = _parse_int_or_default(job.get("cryptkeeper_threads"), 1)
    job["cryptkeeper_rbs_score_cutoff"] = _parse_float_or_default(
        job.get("cryptkeeper_rbs_score_cutoff"), 500.0
    )
    job["cryptkeeper_fail_score"] = _parse_float_or_default(
        job.get("cryptkeeper_fail_score"), 6000.0
    )
    job["cryptkeeper_ignore_first_nt"] = _parse_int_or_default(
        job.get("cryptkeeper_ignore_first_nt"), 30
    )
    job["cryptkeeper_ignore_last_nt"] = _parse_int_or_default(
        job.get("cryptkeeper_ignore_last_nt"), 30
    )

    job["cryptkeeper_pool_factor"] = _parse_float_or_default(
        job.get("cryptkeeper_pool_factor"), 3.0
    )
    job["cryptkeeper_max_pool"] = _parse_int_or_default(
        job.get("cryptkeeper_max_pool"), 200
    )

    job["max_tries_per_replicate"] = _parse_int_or_default(
        job.get("max_tries_per_replicate"), 25
    )
    job["kleinbub_search_limit"] = _parse_int_or_default(
        job.get("kleinbub_search_limit"), 200000
    )
    job["backtrack_window"] = _parse_int_or_default(
        job.get("backtrack_window"), 10
    )

    job["codon_table_sheet"] = _normalize_sheet_name(job.get("codon_table_sheet"))

    job["sequence"] = str(job.get("sequence", "")).strip()

    # IMPORTANT: don't overwrite a loaded codon table (dict) with a string
    if isinstance(job.get("codon_table"), dict):
        # already loaded; leave as-is
        pass
    else:
        job["codon_table"] = str(job.get("codon_table", "")).strip()

    if not job["sequence"]:
        raise InputFormatError("Missing required 'sequence' parameter (path to FASTA/GenBank).")
    ct = job.get("codon_table")
    if not ct:
        raise InputFormatError("Missing required 'codon_table' parameter (path to codon table XLSX).")

# If codon_table is a dict, it's loaded; if it's a string, it's a path (still OK here).
    if not isinstance(ct, (dict, str)):
        raise InputFormatError(f"codon_table must be a path string or a loaded dict; got {type(ct).__name__}")


    seq_path = job["sequence"]
    if not os.path.exists(seq_path):
        raise InputFormatError(f"Sequence file not found: {seq_path}")

    logger.debug(
        f"Prepared job: seq={seq_path}, n={job['n']}, mode={job['optimization_mode']}, "
        f"tries={job['max_tries_per_replicate']}, search_limit={job['kleinbub_search_limit']}, "
        f"backtrack_window={job['backtrack_window']}, min_codon_fraction={job['min_codon_fraction']}, "
        f"cryptkeeper={job['cryptkeeper_enable']}, pool_factor={job['cryptkeeper_pool_factor']}, "
        f"max_pool={job['cryptkeeper_max_pool']}, ck_cutoff={job['cryptkeeper_rbs_score_cutoff']}, "
        f"ck_fail_score={job['cryptkeeper_fail_score']}, "
        f"ignore_first_nt={job['cryptkeeper_ignore_first_nt']}, "
        f"ignore_last_nt={job['cryptkeeper_ignore_last_nt']}, "
        f"gc_min={job['gc_min']}, gc_max={job['gc_max']}, "
        f"max_homopolymer={job['max_homopolymer']}, "
        f"default_motifs_added={len(default_motifs)}"
    )

    return job


def pick_stop_codon(codon_table: Dict[str, Any]) -> str:
    """
    Pick highest-frequency stop codon from codon table.
    Falls back to TAA if none found.
    """
    for key in ("*", "STOP", "Stop", "stop"):
        rows = codon_table.get(key)
        if rows:
            best = max(rows, key=lambda x: float(x[1]))
            return str(best[0]).upper()
    return "TAA"


# ---------------------------
# Candidate generation helpers
# ---------------------------

def _stable_u32(s: str) -> int:
    """Stable 32-bit hash of a string (independent of Python's built-in hash randomization)."""
    h = hashlib.blake2b(s.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(h, "little") & 0xFFFFFFFF


def _candidate_seed(
    job_seed: Optional[int],
    job_idx: int,
    record_id: str,
    rep: int,
    attempt: int,
) -> Optional[int]:
    """Deterministically derive a per-candidate seed from job + record + replicate + attempt."""
    if job_seed is None:
        return None
    rid = _stable_u32(record_id)
    # Mix fields with large multipliers to reduce collisions while keeping determinism.
    return int(job_seed) + rid + (job_idx * 1_000_003) + (rep * 100_003) + int(attempt)


def _generate_candidate(
    protein: str,
    job: Dict[str, Any],
    logger,
    job_idx: int,
    record_id: str,
    rep: int,
    attempt: int,
) -> Tuple[Optional[str], Optional[str]]:

    # Ensure codon table is loaded (dict) before optimization
    ct = job.get("codon_table")
    if isinstance(ct, str):
        try:
            ct = load_codon_table_from_xlsx(ct, sheet_name=job.get("codon_table_sheet"))
            job["codon_table"] = ct  # cache it
        except Exception as e:
            return None, f"Failed to load codon table from path {job.get('codon_table')!r}: {e}"

    if not isinstance(ct, dict):
        return None, f"codon_table must be a dict after loading; got type={type(ct).__name__}"

    cand_seed = _candidate_seed(job.get("seed"), job_idx=job_idx, record_id=record_id, rep=rep, attempt=attempt)

    try:
        optimized = optimize_sequence(
            dna=None,
            protein=protein,
            codon_table=ct,  # <-- use the verified dict
            avoid_codons=job["avoid_codons"],
            avoid_motifs=job["avoid_motifs"],
            max_homopolymer=job["max_homopolymer"],
            gc_min=job.get("gc_min"),
            gc_max=job.get("gc_max"),
            logger=logger,
            optimization_mode=job["optimization_mode"],
            max_attempts=job["kleinbub_search_limit"],
            backtrack_window=job["backtrack_window"],
            min_codon_fraction=job["min_codon_fraction"],
            seed=cand_seed,
        )
        return optimized, None
    except ConstraintError as e:
        return None, str(e)


# ---------------------------
# Main
# ---------------------------

def main():
    args = get_cli_args()
    logger = setup_logger(getattr(args, "log_file", None), getattr(args, "verbose", False))
    logger.info("Starting codonopt")

    os.makedirs(args.out, exist_ok=True)

    # Determine jobs source
    if getattr(args, "batch_table", None):
        job_rows = iter_batch_table(args.batch_table)
    else:
        if not getattr(args, "sequence", None):
            raise InputFormatError("Single mode requires --sequence")
        if not getattr(args, "codon_table", None):
            raise InputFormatError("Single mode requires --codon-table (path to .xlsx)")

        job_rows = [{
            "sequence": args.sequence,
            "codon_table": args.codon_table,
            "codon_table_sheet": getattr(args, "codon_table_sheet", None),

            "optimization_mode": getattr(args, "optimization_mode", "kleinbub"),
            "avoid_codons": getattr(args, "avoid_codons", ""),
            "avoid_motifs": getattr(args, "avoid_motifs", ""),
            "gc_min": getattr(args, "gc_min", None),
            "gc_max": getattr(args, "gc_max", None),
            "seed": getattr(args, "seed", None),
            "n": getattr(args, "n", 1),
            "max_homopolymer": getattr(args, "max_homopolymer", 5),

            "max_tries_per_replicate": getattr(args, "max_tries_per_replicate", 25),
            "kleinbub_search_limit": getattr(args, "kleinbub_search_limit", 200000),
            "backtrack_window": getattr(args, "backtrack_window", 10),

            "min_codon_fraction": getattr(args, "min_codon_fraction", 0.05),

            "cryptkeeper_enable": bool(getattr(args, "cryptkeeper_enable", False)),
            "cryptkeeper_rbs_score_cutoff": getattr(args, "cryptkeeper_rbs_score_cutoff", 500.0),
            "cryptkeeper_threads": getattr(args, "cryptkeeper_threads", 1),
            "cryptkeeper_fail_score": getattr(args, "cryptkeeper_fail_score", 6000.0),
            "cryptkeeper_ignore_first_nt": getattr(args, "cryptkeeper_ignore_first_nt", 30),
            "cryptkeeper_ignore_last_nt": getattr(args, "cryptkeeper_ignore_last_nt", 30),

            "cryptkeeper_pool_factor": getattr(args, "cryptkeeper_pool_factor", 3.0),
            "cryptkeeper_max_pool": getattr(args, "cryptkeeper_max_pool", 200),
        }]

    fasta_out = os.path.join(args.out, "optimized_sequences.fasta")
    tsv_out = os.path.join(args.out, "metrics.tsv")

    all_records: List[SeqRecord] = []
    metrics_rows: List[Dict[str, Any]] = []

    input_record_count = 0

    for job_idx, row in enumerate(tqdm(job_rows, desc="Processing jobs"), start=1):
        job = merge_defaults(args, row, logger)

        records = list(SeqIO.parse(job["sequence"], "fasta"))
        if not records:
            logger.warning(f"No FASTA records in {job['sequence']}")
            continue

        for rec in records:
            input_record_count += 1
            if input_record_count > 1000:
                raise CodonoptError("Maximum of 1000 input sequences exceeded")

            raw_seq = str(rec.seq).strip().upper()
            if not raw_seq:
                continue

            # Convert to protein
            if is_protein_sequence(raw_seq):
                protein = raw_seq
            else:
                # DNA/RNA CDS input: validate, translate, reject internal stops, strip terminal stop if present
                if len(raw_seq) % 3 != 0:
                    raise InputFormatError(
                        f"DNA record '{rec.id}' length ({len(raw_seq)}) not divisible by 3"
                    )

                # Translate full length (do not stop early)
                protein = str(Seq(raw_seq).translate(to_stop=False))

                # If there are internal stops, that's not a valid CDS for codon optimization
                if "*" in protein[:-1]:
                    raise InputFormatError(
                        f"DNA record '{rec.id}' translates with internal stop codon(s). "
                        "Provide a valid CDS (no internal stops) or a protein FASTA."
                    )

                # Allow terminal stop codon and strip it to optimize the coding region only
                had_terminal_stop = False
                if protein.endswith("*"):
                    had_terminal_stop = True
                    protein = protein[:-1]


            # ---------------------------
            # Two-phase logic when CryptKeeper is enabled
            # ---------------------------
            if job["cryptkeeper_enable"]:
                pool_target = int(math.ceil(job["n"] * float(job["cryptkeeper_pool_factor"])))
                pool_target = max(pool_target, job["n"])
                pool_target = min(pool_target, int(job["cryptkeeper_max_pool"]))

                logger.info(
                    f"[{rec.id}] CryptKeeper enabled: generating pool_target={pool_target} "
                    f"(n={job['n']}, pool_factor={job['cryptkeeper_pool_factor']}, max_pool={job['cryptkeeper_max_pool']})"
                )

                # Phase 1: generate pool without CryptKeeper
                pool: List[Dict[str, Any]] = []
                gen_attempts = 0
                max_gen_attempts = pool_target * job["max_tries_per_replicate"]

                while len(pool) < pool_target and gen_attempts < max_gen_attempts:
                    gen_attempts += 1
                    pseudo_rep = (len(pool) + 1)
                    optimized, _fail = _generate_candidate(
                        protein=protein,
                        job=job,
                        logger=logger,
                        job_idx=job_idx,
                        record_id=rec.id,
                        rep=pseudo_rep,
                        attempt=gen_attempts,
                    )
                    if optimized is None:
                        continue

                    # NEW: ensure optimized DNA translates back to the same protein
                    ok, reason, aa = verify_backtranslation_matches_protein(
                                        optimized,
                                        protein,
                                        allow_terminal_stop=job.get("include_stop_codon", False)
                                    )
                    if not ok:
                        logger.debug(
                            f"[{rec.id}] Rejecting candidate (AA mismatch): reason={reason}, "
                            f"translated_len={len(aa) if aa else 'NA'}"
                        )
                        continue

                    pool.append({
                        "dna": optimized,
                        "gen_attempt": gen_attempts,
                    })

                logger.info(f"[{rec.id}] Pool generated: {len(pool)}/{pool_target} candidates (attempts={gen_attempts})")

                # Phase 2: CryptKeeper screen pool until n passing
                passing: List[Dict[str, Any]] = []
                screened = 0

                for item in pool:
                    screened += 1
                    dna = item["dna"]

                    passes, ck_info = run_cryptkeeper_screen(
                        dna_seq=dna,
                        rbs_score_cutoff=job["cryptkeeper_rbs_score_cutoff"],
                        threads=job["cryptkeeper_threads"],
                        fail_score=job["cryptkeeper_fail_score"],
                        ignore_first_nt=job["cryptkeeper_ignore_first_nt"],
                        ignore_last_nt=job["cryptkeeper_ignore_last_nt"],
                        logger=logger,
                    )

                    item["cryptkeeper_pass"] = passes
                    item.update(ck_info)

                    if passes:
                        passing.append(item)
                        if len(passing) >= job["n"]:
                            break

                logger.info(
                    f"[{rec.id}] CryptKeeper screened={screened}, passing={len(passing)}/{job['n']} "
                    f"(fail_score={job['cryptkeeper_fail_score']}, cutoff={job['cryptkeeper_rbs_score_cutoff']}, "
                    f"ignore_first_nt={job['cryptkeeper_ignore_first_nt']}, ignore_last_nt={job['cryptkeeper_ignore_last_nt']})"
                )

                # Emit FASTA + TSV rows for passing (replicates 1..n)
                for rep in range(1, job["n"] + 1):
                    seq_id = f"{rec.id}|job{job_idx:04d}|rep{rep:03d}"

                    if rep <= len(passing):
                        dna_core = passing[rep - 1]["dna"]

                        # Optionally append terminal stop codon (excluded from constraints by design)
                        if job.get("include_stop_codon", False):
                            stop = pick_stop_codon(job["codon_table"])
                            dna_out = dna_core + stop
                        else:
                            dna_out = dna_core

                        # Belt-and-suspenders AA check on the ACTUAL emitted DNA
                        ok, reason, _aa = verify_backtranslation_matches_protein(
                            dna_out,
                            protein,
                            allow_terminal_stop=job.get("include_stop_codon", False),
                        )
                        if not ok:
                            logger.error(
                                f"[{seq_id}] Internal error: candidate passed but AA check failed at emit-time: {reason}"
                            )
                            metrics_rows.append({
                                "sequence_id": seq_id,
                                "source_record": rec.id,
                                "source_file": job["sequence"],
                                "job_index": job_idx,
                                "replicate": rep,

                                "optimization_mode": job["optimization_mode"],
                                "attempts_used": passing[rep - 1].get("gen_attempt", ""),
                                "failure_reason": f"AA identity check failed at emit-time: {reason}",

                                "length": "",
                                "gc_content": "",
                                "max_homopolymer": "",

                                "gc_min": job.get("gc_min"),
                                "gc_max": job.get("gc_max"),
                                "max_homopolymer_limit": job["max_homopolymer"],
                                "avoid_codons": "|".join(job["avoid_codons"]),
                                "avoid_motifs": "|".join(job["avoid_motifs"]),

                                "min_codon_fraction": job["min_codon_fraction"],

                                "codon_table_path": job.get("codon_table_path", ""),
                                "codon_table_sheet": job.get("codon_table_sheet") or "",

                                "cryptkeeper_enable": 1,
                                "cryptkeeper_rbs_score_cutoff": job["cryptkeeper_rbs_score_cutoff"],
                                "cryptkeeper_fail_score": job["cryptkeeper_fail_score"],
                                "cryptkeeper_threads": job["cryptkeeper_threads"],
                                "cryptkeeper_ignore_first_nt": job["cryptkeeper_ignore_first_nt"],
                                "cryptkeeper_ignore_last_nt": job["cryptkeeper_ignore_last_nt"],

                                "cryptkeeper_pool_factor": job["cryptkeeper_pool_factor"],
                                "cryptkeeper_max_pool": job["cryptkeeper_max_pool"],
                                "cryptkeeper_pool_target": pool_target,
                                "cryptkeeper_pool_generated": len(pool),
                                "cryptkeeper_pool_screened": screened,

                                "cryptkeeper_internal_site_count": "",
                                "cryptkeeper_internal_site_positions": "",
                                "cryptkeeper_max_internal_score": "",
                            })
                            continue

                        # SUCCESS: write output + metrics using dna_out (not optimized)
                        all_records.append(SeqRecord(Seq(dna_out), id=seq_id, description=""))

                        metrics_rows.append({
                            "sequence_id": seq_id,
                            "source_record": rec.id,
                            "source_file": job["sequence"],
                            "job_index": job_idx,
                            "replicate": rep,

                            "optimization_mode": job["optimization_mode"],
                            "attempts_used": passing[rep - 1].get("gen_attempt", ""),
                            "failure_reason": "",

                            "length": len(dna_out),
                            "gc_content": calculate_gc(dna_out),
                            "max_homopolymer": max_homopolymer_length(dna_out),

                            "gc_min": job.get("gc_min"),
                            "gc_max": job.get("gc_max"),
                            "max_homopolymer_limit": job["max_homopolymer"],
                            "avoid_codons": "|".join(job["avoid_codons"]),
                            "avoid_motifs": "|".join(job["avoid_motifs"]),

                            "min_codon_fraction": job["min_codon_fraction"],

                            "codon_table_path": job.get("codon_table_path", ""),
                            "codon_table_sheet": job.get("codon_table_sheet") or "",

                            "cryptkeeper_enable": 1,
                            "cryptkeeper_rbs_score_cutoff": job["cryptkeeper_rbs_score_cutoff"],
                            "cryptkeeper_fail_score": job["cryptkeeper_fail_score"],
                            "cryptkeeper_threads": job["cryptkeeper_threads"],
                            "cryptkeeper_ignore_first_nt": job["cryptkeeper_ignore_first_nt"],
                            "cryptkeeper_ignore_last_nt": job["cryptkeeper_ignore_last_nt"],

                            "cryptkeeper_pool_factor": job["cryptkeeper_pool_factor"],
                            "cryptkeeper_max_pool": job["cryptkeeper_max_pool"],
                            "cryptkeeper_pool_target": pool_target,
                            "cryptkeeper_pool_generated": len(pool),
                            "cryptkeeper_pool_screened": screened,

                            "cryptkeeper_internal_site_count": passing[rep - 1].get("cryptkeeper_internal_site_count", ""),
                            "cryptkeeper_internal_site_positions": passing[rep - 1].get("cryptkeeper_internal_site_positions", ""),
                            "cryptkeeper_max_internal_score": passing[rep - 1].get("cryptkeeper_max_internal_score", ""),
                        })

                    else:
                        metrics_rows.append({
                            "sequence_id": seq_id,
                            "source_record": rec.id,
                            "source_file": job["sequence"],
                            "job_index": job_idx,
                            "replicate": rep,

                            "optimization_mode": job["optimization_mode"],
                            "attempts_used": "",
                            "failure_reason": (
                                f"CryptKeeper enabled: only {len(passing)}/{job['n']} passing sequences found "
                                f"after pool generation/screening. Increase cryptkeeper_pool_factor, "
                                f"cryptkeeper_max_pool, max_tries_per_replicate, raise cryptkeeper_fail_score, "
                                f"or relax other constraints."
                            ),

                            "length": "",
                            "gc_content": "",
                            "max_homopolymer": "",

                            "gc_min": job.get("gc_min"),
                            "gc_max": job.get("gc_max"),
                            "max_homopolymer_limit": job["max_homopolymer"],
                            "avoid_codons": "|".join(job["avoid_codons"]),
                            "avoid_motifs": "|".join(job["avoid_motifs"]),

                            "min_codon_fraction": job["min_codon_fraction"],

                            "codon_table_path": job.get("codon_table_path", ""),
                            "codon_table_sheet": job.get("codon_table_sheet") or "",

                            "cryptkeeper_enable": 1,
                            "cryptkeeper_rbs_score_cutoff": job["cryptkeeper_rbs_score_cutoff"],
                            "cryptkeeper_fail_score": job["cryptkeeper_fail_score"],
                            "cryptkeeper_threads": job["cryptkeeper_threads"],
                            "cryptkeeper_ignore_first_nt": job["cryptkeeper_ignore_first_nt"],
                            "cryptkeeper_ignore_last_nt": job["cryptkeeper_ignore_last_nt"],

                            "cryptkeeper_pool_factor": job["cryptkeeper_pool_factor"],
                            "cryptkeeper_max_pool": job["cryptkeeper_max_pool"],
                            "cryptkeeper_pool_target": pool_target,
                            "cryptkeeper_pool_generated": len(pool),
                            "cryptkeeper_pool_screened": screened,

                            "cryptkeeper_internal_site_count": "",
                            "cryptkeeper_internal_site_positions": "",
                            "cryptkeeper_max_internal_score": "",
                        })
        if not job["cryptkeeper_enable"]:
            # ---------------------------
            # Fast path when CryptKeeper is disabled
            # ---------------------------
            for rep in range(1, job["n"] + 1):
                seq_id = f"{rec.id}|job{job_idx:04d}|rep{rep:03d}"

                success = False
                optimized = None
                failure_reason = ""
                attempts_used = 0

                for attempt in range(1, job["max_tries_per_replicate"] + 1):
                    attempts_used = attempt
                    optimized, fail = _generate_candidate(
                        protein=protein,
                        job=job,
                        logger=logger,
                        job_idx=job_idx,
                        record_id=rec.id,
                        rep=rep,
                        attempt=attempt,
                    )
                    if optimized is None:
                        failure_reason = f"attempt {attempt}/{job['max_tries_per_replicate']}: {fail}"
                        continue

                    # NEW: enforce AA identity before accepting
                    ok, reason, aa = verify_backtranslation_matches_protein(
                        optimized,
                        protein,
                        allow_terminal_stop=job.get("include_stop_codon", False)
                    )
                    if not ok:
                        failure_reason = f"attempt {attempt}/{job['max_tries_per_replicate']}: AA check failed: {reason}"
                        continue

                    success = True
                    failure_reason = ""
                    break

                if success and optimized is not None:
                    if job.get("include_stop_codon", False):
                        optimized = optimized + pick_stop_codon(job["codon_table"])
                    all_records.append(SeqRecord(Seq(optimized), id=seq_id, description=""))
                    metrics_rows.append({
                        "sequence_id": seq_id,
                        "source_record": rec.id,
                        "source_file": job["sequence"],
                        "job_index": job_idx,
                        "replicate": rep,

                        "optimization_mode": job["optimization_mode"],
                        "attempts_used": attempts_used,
                        "failure_reason": "",

                        "length": len(optimized),
                        "gc_content": calculate_gc(optimized),
                        "max_homopolymer": max_homopolymer_length(optimized),

                        "gc_min": job.get("gc_min"),
                        "gc_max": job.get("gc_max"),
                        "max_homopolymer_limit": job["max_homopolymer"],
                        "avoid_codons": "|".join(job["avoid_codons"]),
                        "avoid_motifs": "|".join(job["avoid_motifs"]),

                        "min_codon_fraction": job["min_codon_fraction"],

                        "codon_table_path": job.get("codon_table_path", ""),
                        "codon_table_sheet": job.get("codon_table_sheet") or "",

                        "cryptkeeper_enable": 0,
                        "cryptkeeper_rbs_score_cutoff": job["cryptkeeper_rbs_score_cutoff"],
                        "cryptkeeper_fail_score": job["cryptkeeper_fail_score"],
                        "cryptkeeper_threads": job["cryptkeeper_threads"],
                        "cryptkeeper_ignore_first_nt": job["cryptkeeper_ignore_first_nt"],
                        "cryptkeeper_ignore_last_nt": job["cryptkeeper_ignore_last_nt"],

                        "cryptkeeper_pool_factor": job["cryptkeeper_pool_factor"],
                        "cryptkeeper_max_pool": job["cryptkeeper_max_pool"],
                        "cryptkeeper_pool_target": "",
                        "cryptkeeper_pool_generated": "",
                        "cryptkeeper_pool_screened": "",

                        "cryptkeeper_internal_site_count": "",
                        "cryptkeeper_internal_site_positions": "",
                        "cryptkeeper_max_internal_score": "",
                    })
                else:
                    logger.error(f"Failure for {seq_id}: {failure_reason}")
                    metrics_rows.append({
                        "sequence_id": seq_id,
                        "source_record": rec.id,
                        "source_file": job["sequence"],
                        "job_index": job_idx,
                        "replicate": rep,

                        "optimization_mode": job["optimization_mode"],
                        "attempts_used": attempts_used,
                        "failure_reason": failure_reason,

                        "length": "",
                        "gc_content": "",
                        "max_homopolymer": "",

                        "gc_min": job.get("gc_min"),
                        "gc_max": job.get("gc_max"),
                        "max_homopolymer_limit": job["max_homopolymer"],
                        "avoid_codons": "|".join(job["avoid_codons"]),
                        "avoid_motifs": "|".join(job["avoid_motifs"]),

                        "min_codon_fraction": job["min_codon_fraction"],

                        "codon_table_path": job.get("codon_table_path", ""),
                        "codon_table_sheet": job.get("codon_table_sheet") or "",

                        "cryptkeeper_enable": 0,
                        "cryptkeeper_rbs_score_cutoff": job["cryptkeeper_rbs_score_cutoff"],
                        "cryptkeeper_fail_score": job["cryptkeeper_fail_score"],
                        "cryptkeeper_threads": job["cryptkeeper_threads"],
                        "cryptkeeper_ignore_first_nt": job["cryptkeeper_ignore_first_nt"],
                        "cryptkeeper_ignore_last_nt": job["cryptkeeper_ignore_last_nt"],

                        "cryptkeeper_pool_factor": job["cryptkeeper_pool_factor"],
                        "cryptkeeper_max_pool": job["cryptkeeper_max_pool"],
                        "cryptkeeper_pool_target": "",
                        "cryptkeeper_pool_generated": "",
                        "cryptkeeper_pool_screened": "",

                        "cryptkeeper_internal_site_count": "",
                        "cryptkeeper_internal_site_positions": "",
                        "cryptkeeper_max_internal_score": "",
                    })

    # Write FASTA (successful only)
    SeqIO.write(all_records, fasta_out, "fasta")
    logger.info(f"Wrote {len(all_records)} sequences to {fasta_out}")

    # Write TSV
    fieldnames = [
        "sequence_id",
        "source_record",
        "source_file",
        "job_index",
        "replicate",

        "optimization_mode",
        "attempts_used",
        "failure_reason",

        "length",
        "gc_content",
        "max_homopolymer",

        "gc_min",
        "gc_max",
        "max_homopolymer_limit",
        "avoid_codons",
        "avoid_motifs",

        "min_codon_fraction",

        "codon_table_path",
        "codon_table_sheet",

        "cryptkeeper_enable",
        "cryptkeeper_rbs_score_cutoff",
        "cryptkeeper_fail_score",
        "cryptkeeper_threads",
        "cryptkeeper_ignore_first_nt",
        "cryptkeeper_ignore_last_nt",

        "cryptkeeper_pool_factor",
        "cryptkeeper_max_pool",
        "cryptkeeper_pool_target",
        "cryptkeeper_pool_generated",
        "cryptkeeper_pool_screened",

        "cryptkeeper_internal_site_count",
        "cryptkeeper_internal_site_positions",
        "cryptkeeper_max_internal_score",
    ]

    with open(tsv_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in metrics_rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})

    logger.info(f"Wrote metrics to {tsv_out}")
    logger.info("codonopt finished successfully")


if __name__ == "__main__":
    main()