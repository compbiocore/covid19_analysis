import json
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

genes = ["E","M","N","ORF10","ORF14","ORF1a","ORF1b","ORF3a","ORF6","ORF7a","ORF7b","ORF8","ORF9b","S"]

genemap = pd.read_csv("genes.csv", index_col="gene") ### could change this to read in the gff once everything else is functional

nt_ref = next(SeqIO.parse("src/reference.fasta", "fasta")).seq

mutations = defaultdict(list)

for gene in genes:
    ref = nt_ref[genemap.loc[gene, "start"]-1:genemap.loc[gene, "end"]].translate()
    assert "*" not in ref[:-1], ref
    for seq in SeqIO.parse(f"results/nextalign/nextalign_gene_{gene}.translation.fasta", "fasta"):
        assert len(seq) == len(ref), seq.id
        for i, (aa0, aa1) in enumerate(zip(ref, seq.seq), start=1):
            assert aa0 != "-" and aa0 != "X"
            if aa1 != "X" and aa1 != aa0:
                mutations[seq.id].append(f"{gene}:{aa0}{i}{aa1}")

with open("results/mutations.json", "w") as f:
    json.dump(mutations, f, sort_keys=True, indent=2)