import pandas as pd
from Bio import SeqIO
import datetime
import sys


ri = pd.read_csv(sys.argv[2] + "/0_data/gisaid.tsv", sep="\t").sort_values("date")
ri = ri.drop_duplicates(subset = "strain")

# Join Pangolin
#pangolin = pd.read_csv(sys.argv[2] + "/3_results/" + sys.argv[1] + "/pangolin/lineage_report.csv", dtype='str')
##pangolin = pangolin.drop_duplicates(subset = "taxon")
##columns = dict((col, "pangolin." + col) for col in pangolin.columns)
##pangolin.rename(columns=columns, inplace ="True")
#pangolin.rename(columns = {'pangolin.taxon':'strain'}, inplace="True")
#pangolin.rename(columns = {'pangolin.status':'pangolin.qc_status'}, inplace="True")
#print(pangolin)
#ri = ri.merge(pangolin, how="outer", on="strain", validate="1:1")

pangolin = pd.read_csv(sys.argv[2] + "/3_results/" + sys.argv[1] + "/pangolin/lineage_report.csv", dtype='str')
columns = dict((col, "pangolin." + col) for col in pangolin.columns)
pangolin.rename(columns = {'taxon':'strain'}, inplace="True")
pangolin = pangolin.drop_duplicates(subset = "strain")
ri = ri.merge(pangolin.rename(columns=columns), how="outer", on="strain", validate="1:1")

# Join NextClade
nextclade = pd.read_csv(sys.argv[2] + "/3_results/" + sys.argv[1] + "/nextclade.tsv", sep="\t", dtype='str')
nextclade = nextclade.drop_duplicates(subset = "seqName")
columns = dict((col, "nextclade." + col) for col in nextclade.columns)
columns["seqName"] = "strain"
ri = ri.merge(nextclade.rename(columns=columns), how="outer", on="strain", validate="1:1")

# Join NextStrain
#nextstrain = pd.read_csv(sys.argv[2] + "/3_results/${day}/nextstrain-diagnostics-flagged.tsv", sep="\t", usecols=["strain", "flagging_reason"]).rename(columns={"flagging_reason": "nextstrain.flagging_reason"})
#ri = ri.merge(nextstrain, how="outer", on="strain", validate="1:1")

# Join CDC
cdc = pd.read_csv(sys.argv[2] + "/2_metadata/cdc-voc-vbm-variants-ncbi02142023.txt", sep="\t", dtype='str')
ri = ri.merge(cdc, how="left", on="pangolin.lineage")

failed = (
    (ri["pangolin.qc_status"] != "pass") |
    (
        (
            ri["strain"].str.startswith("hCoV-19/USA/RI_RKL") |
            ri["strain"].str.startswith("hCoV-19/USA/RI-RISHL")
        ) &
        (
            (ri["nextclade.qc.overallStatus"] == "bad") #|
            #(ri["nextstrain.flagging_reason"].notnull())
        )
    )
)

seq_len = {}

# Filter sequences
passed = frozenset(ri[~failed]["strain"])
with open(sys.argv[2] + "/3_results/" + sys.argv[1] + "/qc.fa", "w") as f:
    for record in SeqIO.parse(sys.argv[2] + "/0_data/gisaid.fasta", "fasta"):
        seq_len[record.id] = sum(1 for nt in str(record.seq).upper() if nt != "-" and nt != "N")
        if record.id in passed:
            print(">"+record.id, file=f)
            print(record.seq, file=f)

ri = ri.merge(
  pd.DataFrame({"strain": list(seq_len.keys()), "seq_len": list(seq_len.values())}),
  how="left",
  on="strain"
)

ri[failed].to_csv(sys.argv[2] + "/3_results/" + sys.argv[1] + "/qc-failed.csv", index=False)
ri[~failed].to_csv(sys.argv[2] + "/3_results/" + sys.argv[1] + "/qc-passed.csv", index=False)

