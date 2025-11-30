import tarfile, io, pandas as pd
tar_path = "My_MR_Project/clue_result.tar.gz"
targets = [
    "my_analysis.sig_queryl1k_tool.69293c37bf1c3c00136dd55a/gsea/TAG/arfs/NORM_CS/gsea_result.gct",
    "my_analysis.sig_queryl1k_tool.69293c37bf1c3c00136dd55a/ncs.gct",
    "my_analysis.sig_queryl1k_tool.69293c37bf1c3c00136dd55a/arfs/TAG/query_result.gct",
]
with tarfile.open(tar_path, "r:gz") as tar:
    for name in targets:
        try:
            f = tar.extractfile(name)
            s = f.read().decode("utf-8", errors="replace")
            lines = s.splitlines()
            header_row = 0
            for i, line in enumerate(lines[:50]):
                if line.startswith("id") or line.startswith("cid") or line.startswith("Name\t"):
                    header_row = i
                    break
            df = pd.read_csv(io.StringIO(s), sep="\t", skiprows=header_row)
            print(name)
            print(list(df.columns))
        except Exception as e:
            print(name, "error", e)
