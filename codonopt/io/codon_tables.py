import pandas as pd
from collections import defaultdict

def load_codon_table_xlsx(path, sheet):
    df = pd.read_excel(path, sheet_name=sheet)

    required = {"AA", "Codon", "Fraction"}
    if not required.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required}")

    table = defaultdict(list)
    for _, row in df.iterrows():
        aa = row["AA"]
        codon = row["Codon"].upper()
        frac = float(row["Fraction"])
        table[aa].append((codon, frac))

    # Normalize fractions
    for aa in table:
        total = sum(f for _, f in table[aa])
        table[aa] = [(c, f / total) for c, f in table[aa]]

    return dict(table)
