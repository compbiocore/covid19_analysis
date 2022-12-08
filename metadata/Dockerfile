FROM staphb/pangolin

#Install software
RUN pip install snakemake
RUN apt-get update; apt-get install curl -y
RUN mkdir -p bin
RUN curl -fsSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu" -o "bin/nextclade" && chmod +x bin/nextclade
RUN curl -fsSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextalign-x86_64-unknown-linux-gnu" -o "bin/nextalign" && chmod +x bin/nextalign

COPY ../scripts/run.sh /data/
COPY ri_sequences.fa /data/
COPY src /data/src
