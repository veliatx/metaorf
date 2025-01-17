FROM ubuntu:20.04
LABEL maintainer="Stephen Federowicz (steve@veliatx.com)"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends --fix-missing wget git unzip pip tar gzip bzip2 pigz \
	less gcc g++ vim samtools dirmngr gnupg apt-transport-https ca-certificates software-properties-common \
	libbz2-dev libcurl4-openssl-dev libgsl0-dev liblzma-dev libncurses5-dev libperl-dev libssl-dev zlib1g-dev make \
	libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff5-dev

# install miniconda
RUN wget -P /tmp \
    "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH

# install samtools
COPY samtools-1.19.tar.bz2 samtools-1.19.tar.bz2
RUN tar -xvf samtools-1.19.tar.bz2
RUN cd samtools-1.19 && make && make install

# install STAR and fastp
COPY 2.7.4a.tar.gz 2.7.4a.tar.gz
RUN tar -xzf 2.7.4a.tar.gz
RUN cd /usr/local/bin && \
    wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp
ENV PATH=$PATH:/fastp:/STAR-2.7.4a/bin/Linux_x86_64_static

# install python libs
RUN pip install git-remote-codecommit numpy==1.23.5 scikit-learn igraph tables brewer2mpl awscli umi_tools

# install R and ORFik
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt install -y r-base
RUN Rscript -e 'install.packages(c("pacman","tibble","conflicted","R.utils","tidyverse","multidplyr","ineq","BiocManager"))'
RUN Rscript -e 'BiocManager::install(c("IRanges", "GenomicRanges", "GenomicFeatures", "ORFik"))'
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install_github("Roleren/ORFik")'

WORKDIR /usr/src/
COPY src src

RUN chmod -R 777 /usr/src/
    
