FROM ubuntu:20.04
LABEL maintainer="Stephen Federowicz (steve@veliatx.com)"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends --fix-missing wget git unzip pip tar gzip bzip2 pigz \
	less gcc g++ vim samtools dirmngr gnupg apt-transport-https ca-certificates software-properties-common \
	libbz2-dev libcurl4-openssl-dev libgsl0-dev liblzma-dev libncurses5-dev libperl-dev libssl-dev zlib1g-dev make \
	libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff5-dev bedtools

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

# install python libs
RUN pip install git-remote-codecommit numpy awscli pandas biopython pyBigWig


WORKDIR /usr/src/
COPY src src
COPY tools tools

RUN chmod -R 777 /usr/src/
    
