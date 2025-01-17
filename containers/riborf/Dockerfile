FROM ubuntu:20.04
RUN DEBIAN_FRONTEND=noninteractive && export DEBIAN_FRONTEND && apt-get update && apt-get -y install dumb-init runit git openssl build-essential perl bowtie2 tophat-recondition r-base cpanminus make gcc g++ libpng-dev uuid-dev libmariadbclient-dev libssl-dev libmariadb-dev-compat

RUN apt-get install -y --no-install-recommends --fix-missing wget
 
#RUN mkdir /work /data && cd /work && git clone https://github.com/zhejilab/RibORF && git clone https://github.com/ucscGenomeBrowser/kent-core
#RUN cd /work/kent-core && make
#COPY GRCh38.primary_assembly.genome.fa.gz /work
#COPY gencode.v28.annotation.gtf.gz /work
#RUN mkdir /var/svc.d
#COPY container-start  /etc
#ENTRYPOINT ["dumb-init", "--single-child", "--rewrite", "15:1", "--verbose", "--"]
#CMD ["/etc/container-start"]

# install miniconda
RUN wget -P /tmp \
    "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH
RUN pip install pandas

WORKDIR /usr/src/
COPY src src
COPY RibORF RibORF

RUN chmod -R 777 /usr/src/
