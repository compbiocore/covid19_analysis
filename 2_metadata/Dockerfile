FROM staphb/pangolin

WORKDIR /home

RUN pip install snakemake
RUN apt-get clean && apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -y --no-install-recommends install --fix-missing \
    apt-file \
    apt-utils \
    build-essential \
    bzip2 \
    ca-certificates \
    curl \
    cmake \
    default-jdk \
    default-jre\
    gcc \
    gdb \
    git \
    hdf5-helpers \
    git \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    lbzip2 \
    libbz2-dev \
    libcairo2-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libfftw3-dev \
    libgeos-dev \
    libgl1-mesa-dev \
    libglpk-dev \
    libglu1-mesa-dev \
    libgsl0-dev \
    libhdf4-alt-dev \
    libhdf5-dev \
    libjpeg-dev \
    libjq-dev \
    liblzma-dev \
    libmariadbd-dev \
    libnetcdf-dev \
    libpng-dev \
    libpq-dev \
    libproj-dev \
    libprotobuf-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libssl-dev \
    libudunits2-dev \
    libxml2-dev \
    libxt-dev \
    libxt6 \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    make \
    nano \
    netcdf-bin \
    postgis \
    protobuf-compiler \
    python3-pip \
    r-base \
    vim \
    sqlite3 \
    tk-dev \
    unixodbc-dev \
    wget \
    unzip \
    libpoppler-cpp-dev \
    zlib1g-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Update pip and install python3 libraries
RUN pip install --upgrade pip
RUN pip3 install --upgrade pip
RUN pip3 install pandas biopython

# Install nextclade and nextalign
RUN cd /usr/local/bin/ && curl -fsSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu" -o "./nextclade" && chmod +x ./nextclade
RUN cd /usr/local/bin/ && curl -fsSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextalign-x86_64-unknown-linux-gnu" -o "./nextalign" && chmod +x ./nextalign

# Install IQtree
WORKDIR /download/IQ-TREE
RUN curl -fsSL https://github.com/iqtree/iqtree2/releases/download/v2.1.2/iqtree-2.1.2-Linux.tar.gz \
  | tar xzvpf - --no-same-owner --strip-components=1 \
 && mv bin/iqtree2 /usr/local/bin/iqtree 
 
# Clean 
WORKDIR /home
RUN rm -r /download/IQ-TREE

#Install R packages
RUN R -e "install.packages(c('tidyverse', 'lubridate', 'scales', 'readxl', 'RColorBrewer'), dependencies=TRUE, lib='/usr/lib/R/library', repos='http://cran.rstudio.com/', ask = FALSE)"
RUN R -e "install.packages(c('openxlsx'), dependencies=TRUE, lib='/usr/lib/R/library', repos='http://cran.rstudio.com/', ask = FALSE)"

