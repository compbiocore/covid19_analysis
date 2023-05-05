import pandas as pd
from Bio import SeqIO

day=$(date "+%Y%m%d")

ri = pd.read_csv("../0_data/gisaid.tsv", sep="\t").sort_values("date")
ri = ri.drop_duplicates(subset = "strain")

# Join Pangolin
pangolin = pd.read_csv("../3_results/${day}/pangolin/lineage_report.csv")
columns = dict((col, "pangolin." + col) for col in pangolin.columns)
pangolin.rename(columns = {'taxon':'strain'}, inplace="True")
pangolin = pangolin.drop_duplicates(subset = "strain")
ri = ri.merge(pangolin.rename(columns=columns), how="outer", on="strain", validate="1:1")

# Join NextClade
nextclade = pd.read_csv("../3_results/${day}/nextclade.tsv", sep="\t")
nextclade = nextclade.drop_duplicates(subset = "seqName")
columns = dict((col, "nextclade." + col) for col in nextclade.columns)
columns["seqName"] = "strain"
ri = ri.merge(nextclade.rename(columns=columns), how="outer", on="strain", validate="1:1")

# Join NextStrain
#nextstrain = pd.read_csv("../3_results/${day}/nextstrain-diagnostics-flagged.tsv", sep="\t", usecols=["strain", "flagging_reason"]).rename(columns={"flagging_reason": "nextstrain.flagging_reason"})
#ri = ri.merge(nextstrain, how="outer", on="strain", validate="1:1")

# Join CDC
cdc = pd.read_csv("../2_metadata/cdc-voc-vbm-variants-ncbi02142023.txt", sep="\t")
ri = ri.merge(cdc, how="left", on="pangolin.lineage")
print(ri["pangolin.qc_status"])
failed = (
    (ri["pangolin.qc_status"] != "pass") |
    (
        (
            ri["strain"].str.startswith("hCoV-19/USA/RI_RKL") |
            ri["strain"].str.startswith("hCoV-19/USA/RI-RISHL")
        ) &
        (
            (ri["nextclade.qc.overallStatus"] == "bad") |
            (ri["nextstrain.flagging_reason"].notnull())
        )
    )
)

seq_len = {}

# Filter sequences
passed = frozenset(ri[~failed]["strain"])
with open("../3_results/${day}/qc.fa", "w") as f:
    for record in SeqIO.parse("../0_data/gisaid.fasta", "fasta"):
        seq_len[record.id] = sum(1 for nt in str(record.seq).upper() if nt != "-" and nt != "N")
        if record.id in passed:
            print(">"+record.id, file=f)
            print(record.seq, file=f)

ri = ri.merge(
  pd.DataFrame({"strain": list(seq_len.keys()), "seq_len": list(seq_len.values())}),
  how="left",
  on="strain"
)

ri[failed].to_csv("../3_results/${day}/qc-failed.csv", index=False)
ri[~failed].to_csv("../3_results/${day}/qc-passed.csv", index=False)

