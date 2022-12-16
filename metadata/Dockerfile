FROM staphb/pangolin

#Install software
RUN pip install snakemake
RUN apt-get update && \
    apt-get -y --no-install-recommends install --fix-missing \
    curl \
    nano \
    vim \
    git \
    unzip
RUN cd /usr/local/bin/ && curl -fsSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu" -o "./nextclade" && chmod +x ./nextclade
RUN cd /usr/local/bin/ && curl -fsSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextalign-x86_64-unknown-linux-gnu" -o "./nextalign" && chmod +x ./nextalign

