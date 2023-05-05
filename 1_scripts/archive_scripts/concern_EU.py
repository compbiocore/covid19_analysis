import json
import pandas as pd
import numpy as np

# 04242023 Only mutations listed on https://www.ecdc.europa.eu/en/covid-19/variants-concern

concern = {
        "143-": ["BA.1", "BA.3"],
        "144-": ["BA.1", "BA.3"],
        "145-": ["BA.1", "BA.3"],
        "212-": ["BA.1", "BA.3", "BA.2"],
        "69-": ["BA.1", "BA.3"],
        "70-": ["BA.1", "BA.3"],
        "A222V": ["AY.4.2"],
        "A653V": ["A.27"],
        "A67V": ["BA.1", "BA.3"],
        "A701V": ["B.1.526", "B.1.351"],
        "D339H": ["BA.2.75"],
        "D405N": ["BA.3", "BA.2"],
        "D614G": ["B.1.1.7", "B.1.427", "B.1.429", "B.1.616", "B.1.525", "P.3", "B.1.617.1", "B.1.620", "B.1.617.3", "B.1.214.2", "C.16", "B.1.526", "B.1.526.1", "B.1.526.2", "P.2", "B.1.1.519", "AV.1", "AT.1", "B.1.621", "C.37", "AY.4.2", "B.1.1.318", "C.1.2", "B.1.351", "P.1", "B.1.640", "B.1.617.2", "BA.1", "BA.3", "BA.2"],
        "D796Y": ["BA.1", "BA.3", "BA.2"],
        "E180V": ["XBB.1.16"],
        "E484A": ["BA.1", "BA.3", "BA.2"],
        "E484K": ["B.1.525", "P.3", "B.1.620", "A.28", "B.1.526", "P.2", "AV.1", "AT.1", "B.1.621", "B.1.1.318", "C.1.2", "B.1.351", "P.1"],
        "E484Q": ["B.1.617.1", "B.1.617.3"],
        "F157L": ["BA.2.75"],
        "F486P": ["XAY", "XBB.1.16"],
        "F486V": ["BF.7", "BA.4", "BA.5"],
        "F490R": ["B.1.640"],
        "F490S": ["C.37", "XBB", "XBB.1.5"],
        "F490S,": ["BN.1"],
        "G142D": ["BA.1", "BA.2"],
        "G257S": ["BA.2.75"],
        "G339D": ["BA.1", "BA.3", "BA.2"],
        "G446S": ["BA.1", "BA.3", "BA.2.75"],
        "G496S": ["BA.1"],
        "G669S": ["B.1.616"],
        "H655Y": ["B.1.616", "A.27", "A.28", "C.1.2", "P.1", "BA.1", "BA.3", "BA.2"],
        "I210V": ["BA.2.75"],
        "K356T": ["BN.1"],
        "K417N": ["B.1.351", "BA.1", "BA.3", "BA.2"],
        "K417T": ["P.1"],
        "K444R": ["BA.2.3.20"],
        "K444T": ["CH.1.1", "BQ.1"],
        "L452M": ["BA.2.3.20"],
        "L452Q": ["C.37"],
        "L452R": ["B.1.427", "B.1.429", "B.1.617.1", "B.1.617.3", "A.27", "C.16", "B.1.526.1", "AY.4.2", "B.1.617.2", "BA.4", "BA.5", "CH.1.1"],
        "L981F": ["BA.1"],
        "N211I": ["BA.1", "BA.3", "BA.2"],
        "N394S": ["B.1.640"],
        "N439K": ["AV.1"],
        "N440K": ["BA.1", "BA.3", "BA.2"],
        "N450K": ["B.1.214.2"],
        "N460K": ["BA.2.3.20", "BA.2.75", "BQ.1", "XBB", "XBB.1.5"],
        "N501T": ["A.28"],
        "N501Y": ["B.1.1.7", "P.3", "A.27", "B.1.621", "C.1.2", "B.1.351", "P.1", "B.1.640", "BA.1", "BA.3", "BA.2"],
        "N679K": ["AT.1", "C.1.2", "BA.1", "BA.3", "BA.2"],
        "N764K": ["BA.1", "BA.2"],
        "N856K": ["BA.1"],
        "N969K": ["BA.1", "BA.3", "BA.2"],
        "P681H": ["B.1.1.7", "P.3", "B.1.620", "AV.1", "B.1.621", "B.1.1.318", "B.1.640", "BA.1", "BA.3", "BA.2"],
        "P681R": ["B.1.617.1", "B.1.617.3", "AY.4.2", "B.1.617.2"],
        "Q414K": ["B.1.214.2"],
        "Q493R": ["BA.1", "BA.3", "BA.2"],
        "Q498R": ["BA.1", "BA.3", "BA.2"],
        "Q677H": ["B.1.525"],
        "Q954H": ["BA.1", "BA.3", "BA.2"],
        "R346K": ["B.1.621"],
        "R346S": ["B.1.640"],
        "R346T": ["BF.7", "BN.1"],
        "R408S": ["BA.2"],
        "R493Q": ["BA.4", "BA.5"],
        "S371F": ["BA.3", "BA.2"],
        "S371L": ["BA.1"],
        "S373P": ["BA.1", "BA.3", "BA.2"],
        "S375F": ["BA.1", "BA.3", "BA.2"],
        "S477N": ["B.1.620", "B.1.526.2", "BA.1", "BA.3", "BA.2"],
        "S486P": ["XBB.1.5"],
        "T376A": ["BA.2"],
        "T478K": ["B.1.1.519", "AY.4.2", "B.1.617.2", "BA.1", "BA.3", "BA.2"],
        "T478R": ["XBB.1.16"],
        "T547K": ["BA.1"],
        "T95I": ["BA.1"],
        "V213G": ["BA.2"],
        "V483A": ["B.1.616"],
        "W152R": ["BA.2.75"],
        "Y145H": ["AY.4.2"],
        "Y449H": ["C.1.2"],
        "Y449N": ["B.1.640"],
        "Y505H": ["BA.1", "BA.3", "BA.2"]
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
detail.to_csv("../3_results/20230313/concern-long_EU.csv", index=False)

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
detail.to_csv("../3_results/20230313/concern_EU.csv", index=False)
print(detail.describe())

