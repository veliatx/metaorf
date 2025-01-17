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

# install python libs
RUN pip install git-remote-codecommit numpy scikit-learn igraph tables brewer2mpl awscli

# install R and ORFquant
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
#RUN apt list -a r-base
RUN apt-get install -y --no-install-recommends r-base r-base-dev cmake
RUN Rscript -e 'install.packages(c("pacman","tibble","conflicted","R.utils","tidyverse","multidplyr","ineq","lme4","BiocManager"))'

RUN apt-get install -y libxml2-dev
RUN Rscript -e 'install.packages("XML", repos = "http://www.omegahat.net/R")'
RUN Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/numDeriv_2016.8-1.1.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.5-0.1.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'install.packages(c("car","ggpubr","rstatix"))'
RUN Rscript -e 'BiocManager::install(c("IRanges", "GenomicRanges", "GenomicFeatures"))'

RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install_github("lcalviell/Ribo-seQC")'
RUN Rscript -e 'devtools::install_github("lcalviell/ORFquant")'

RUN apt-get install -y pandoc
RUN pip install pandas

WORKDIR /usr/src/
COPY src src

RUN chmod -R 777 /usr/src/
    
