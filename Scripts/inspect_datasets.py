import os
import sys
import gzip
from pathlib import Path

def print_header(title):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)

def read_head_text(path, n=5):
    p = Path(path)
    is_gz = p.suffix == ".gz"
    opener = gzip.open if is_gz else open
    try:
        with opener(p, "rt", encoding="utf-8", errors="replace") as f:
            lines = []
            for i in range(n + 1):
                line = f.readline()
                if not line:
                    break
                lines.append(line.rstrip("\n"))
        if not lines:
            print("empty or unreadable")
            return
        header = lines[0]
        delim = "\t" if "\t" in header else "," if "," in header else "whitespace"
        print(f"path: {p}")
        print(f"compressed: {is_gz}")
        print(f"delimiter: {delim}")
        cols = header.split("\t") if delim == "\t" else header.split(",") if delim == "," else header.split()
        print(f"columns_count: {len(cols)}")
        print("columns: " + ", ".join(cols[:50]))
        preview = lines[1:]
        for idx, row in enumerate(preview):
            print(f"row{idx+1}: {row[:500]}")
    except Exception as e:
        print(f"error reading {p}: {e}")

def read_head_binary(path, nbytes=32):
    p = Path(path)
    try:
        with open(p, "rb") as f:
            data = f.read(nbytes)
        print(f"path: {p}")
        print("binary: true")
        print("hex: " + data.hex())
    except Exception as e:
        print(f"error reading {p}: {e}")

def inspect_exposures():
    print_header("Exposure: OneK1K sc-eQTL tables")
    base = Path("My_MR_Project/Exposure")
    targets = [
        base / "cd8et_eqtl_table.tsv.gz",
        base / "cd8nc_eqtl_table.tsv.gz",
        base / "nk_eqtl_table.tsv.gz",
        base / "nkr_eqtl_table.tsv.gz",
        base / "monoc_eqtl_table.tsv.gz",
        base / "mononc_eqtl_table.tsv.gz",
        base / "cd4et_eqtl_table.tsv.gz",
        base / "cd4nc_eqtl_table.tsv.gz",
        base / "plasma_eqtl_table.tsv.gz",
        base / "cd8et_esnp_table.tsv.gz",
    ]
    for p in targets:
        if p.exists():
            read_head_text(p)
        else:
            print(f"missing: {p}")

def inspect_outcomes():
    print_header("Outcome: Lung cancer GWAS (ILCCO and FinnGen)")
    p1 = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
    p2 = Path("My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz")
    alt1 = Path("28604730-GCST004744-EFO_0000571.h.tsv.gz")
    for p in [p1, p2, alt1]:
        if p.exists():
            if p.suffix == ".gz":
                read_head_text(p)
            else:
                read_head_text(p)
        else:
            print(f"missing: {p}")

def inspect_reference():
    print_header("Reference: 1000 Genomes EUR and synonyms")
    bdir = Path("My_MR_Project/Reference")
    bed = bdir / "g1000_eur.bed"
    bim = bdir / "g1000_eur.bim"
    fam = bdir / "g1000_eur.fam"
    for p in [bed, bim, fam]:
        if p.exists():
            if p.suffix == ".bed":
                read_head_binary(p)
            else:
                read_head_text(p)
        else:
            print(f"missing: {p}")
    s = Path("g1000_eur.synonyms")
    if s.exists():
        read_head_text(s)

def inspect_gtex():
    print_header("GTEx v8: Lung and Whole Blood eQTL")
    base = Path("My_MR_Project/GTEx_Analysis_v8_eQTL")
    files = [
        base / "Lung.v8.egenes.txt.gz",
        base / "Lung.v8.signif_variant_gene_pairs.txt.gz",
        base / "Whole_Blood.v8.egenes.txt.gz",
        base / "Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",
    ]
    for p in files:
        if p.exists():
            read_head_text(p)
        else:
            print(f"missing: {p}")
    allpairs_lung = Path("My_MR_Project/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Lung.allpairs.txt.gz")
    if allpairs_lung.exists():
        read_head_text(allpairs_lung)
    else:
        print(f"missing: {allpairs_lung}")

def inspect_extra():
    print_header("Extra: OneK1K external files")
    p = Path("OneK1K_CD4_Naive.tsv.gz")
    if p.exists():
        read_head_text(p)
    else:
        print(f"missing: {p}")

def main():
    inspect_exposures()
    inspect_outcomes()
    inspect_reference()
    inspect_gtex()
    inspect_extra()

if __name__ == "__main__":
    main()

