import json
import pandas as pd
import numpy as np

# 04242023 Updates from https://www.ecdc.europa.eu/en/covid-19/variants-concern and Omicron specific mutaions from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9376347/

concern = {
        "L5F": ["B.1.526"],
        "S13I": ["B.1.429"],
        "L18F": ["P.1"],
        "T19R": ["B.1.617.2", "B.1.617.3"],
        "T20N": ["P.1"],
        "P26S": ["P.1"],
        "A67V": ["B.1.525"],
        "D80G": ["B.1.526.1"],
        "D80A": ["B.1.351"],
        "T95I": ["B.1.526", "B.1.617.1"],
        "D138Y": ["P.1"],
        "G142D": ["B.1.617.1", "B.1.617.2", "B.1.617.3"],
        "W152C": ["B.1.429"],
        "W152R": ["BA.2.75"],
        "E154K": ["B.1.617.1"],
        "F157S": ["B.1.526.1"],
        "F157L": ["BA.2.75"],
        "R158G": ["B.1.617.2"],
        "E180V": ["XBB.1.16"],
        "R190S": ["P.1"],
        "I210V": ["BA.2.75"],
        "D215G": ["B.1.351"],
        "D253G": ["B.1.526"],
        "G257S": ["BA.2.75"],
        "D339H": ["BA.2.75"],
        "G339D": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "S371L": ["B.1.1.529", "BA.1", "XBB.1.5"],
        "S371F": ["BA.2", "BA.4", "BA.5", "XBB"],
        "S373P": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "S375F": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "K417N": ["B.1.351", "B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "K417T": ["P.1"],
        "N440K": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "K444T": ["BQ.1", "CH.1.1"],
        "G446S": ["BA.2.75", "BA.1", "XBB", "XBB.1.5"],
        "L452R": ["B.1.427", "B.1.429", "B.1.526.1", "B.1.617", "B.1.617.1", "B.1.617.2", "B.1.617.3", "CH.1.1"],
        "N460K": ["BA.2.75", "BQ.1", "XBB", "XBB.1.5"],
        "S477N": ["B.1.526", "B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "T478K": ["B.1.617.2", "B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "T478R": ["CH.1.1"],
        "E484A": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "E484K": ["B.1.526", "B.1.525", "P.2", "B.1.1.7", "P.1", "B.1.351"],
        "E484Q": ["B.1.617", "B.1.617.1", "B.1.617.3"],
        "S486P": ["XBB.1.5"],
        "F486P": ["XBB.1.16"],
        "F490S": ["XBB", "XBB.1.5"],
        "Q493": ["BA.2.75"],
        "Q493R": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "S494P": ["B.1.1.7"],
        "G496S": ["B.1.1.529", "BA.1", "XBB.1.5"],
        "Q498R": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "N501Y": ["B.1.1.7", "P.1", "B.1.351", "B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "Y505H": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "F565L": ["P.2"],
        "A570D": ["B.1.1.7"],
        "D614G": ["B.1.1.7"],
        "H655Y": ["P.1", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "Q677H": ["B.1.525"],
        "N679K": ["BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "P681H": ["B.1.1.7", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "P681R": ["B.1.617.1", "B.1.617.2", "B.1.617.3"],
        "A701V": ["B.1.526", "B.1.351"],
        "T716I": ["B.1.1.7"],
        "N764K": ["BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "T791I": ["B.1.526.1"],
        "D796Y": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "N856K": ["BA.1", "XBB.1.5"],
        "T859N": ["B.1.526.1"],
        "F888L": ["B.1.525"],
        "D950H": ["B.1.526.1"],
        "D950N": ["B.1.617.2", "B.1.617.3"],
        "Q954H": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "N969K": ["B.1.1.529", "BA.1", "BA.2", "BA.4", "BA.5", "XBB", "XBB.1.5"],
        "L981F": ["B.1.1.529", "BA.1", "XBB.1.5"],
        "S982A": ["B.1.1.7"],
        "T1027I": ["P.1"],
        "Q1071H": ["B.1.617.1"],
        "D1118H": ["B.1.1.7"],
        "V1176F": ["P.2"],
        "K1191N": ["B.1.1.7"],
        "H69-": ["B.1.525", "B.1.1.7", "B.1.1.529", "BA.1", "XBB.1.5"],
        "V70-": ["B.1.525", "B.1.1.7", "B.1.1.529", "XBB.1.5"],
        "Y144-": ["B.1.525", "B.1.526.1", "B.1.1.7"],
        "E156-": ["B.1.617.2"],
        "F157-": ["B.1.617.2"],
        "L241-": ["B.1.351"]
}
concern = dict(("S:"+key, value) for key, value in concern.items())

mutations = json.load(open("../3_results/20230313/mutations.json"))
seqs = pd.read_csv("../3_results/20230313/qc-passed.csv", usecols=["strain", "date", "pangolin.lineage", "cdc.classification"]).sort_values("date")

detail = []

for _, row in seqs.iterrows():
    for mutation in mutations.get(row.strain, []):
        if mutation in concern:
            if not row["pangolin.lineage"] in concern[mutation]:
                detail.append([row.strain, row.date, mutation])

detail = pd.DataFrame.from_records(detail, columns=["strain", "date", "mutation"])
detail.to_csv("../3_results/20230313/concern-long.csv", index=False)

detail["value"] = 1
print(detail)
detail = detail.drop_duplicates(subset = "strain")
print(detail)
detail = detail.pivot(index="strain", columns="mutation", values="value")\
               .fillna(0)\
               .astype(int)\
               .reset_index()\
               .merge(seqs[["strain", "date"]], how="left", on="strain")\
               .sort_values("date")
detail.to_csv("../3_results/20230313/concern.csv", index=False)
print(detail.describe())

