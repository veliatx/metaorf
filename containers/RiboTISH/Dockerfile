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
RUN pip install pandas ribotish

WORKDIR /usr/src/
COPY src src

RUN chmod -R 777 /usr/src/
    
