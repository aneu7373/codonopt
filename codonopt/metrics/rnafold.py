import subprocess
import tempfile
import re

def rnafold_mfe(seq):
    with tempfile.NamedTemporaryFile("w+", delete=True) as f:
        f.write(seq)
        f.flush()
        result = subprocess.run(
            ["RNAfold", "--noPS"],
            stdin=open(f.name),
            capture_output=True,
            text=True
        )

    match = re.search(r"\((\-?\d+\.\d+)\)", result.stdout)
    return float(match.group(1)) if match else 0.0
