FROM ubuntu:20.04
LABEL maintainer="Stephen Federowicz (steve@veliatx.com)"

RUN apt-get update && \
    apt-get install -y --no-install-recommends --fix-missing wget git unzip tar gzip bzip2 pigz vim python3 pip less default-jre make gcc libz-dev samtools g++
RUN pip install awscli pandas

RUN wget --no-check-certificate https://github.com/erhard-lab/gedi/releases/download/Gedi_1.0.5/Gedi_1.0.5.tar.gz
RUN tar -xvzf Gedi_1.0.5.tar.gz
RUN chmod +x /Gedi_1.0.5/gedi

# RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
# RUN tar -xzf 2.7.10b.tar.gz
#RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/refs/tags/2.7.4a.tar.gz
#RUN tar -xzf 2.7.4a.tar.gz
    
#RUN cd /usr/local/bin && \
#    wget http://opengene.org/fastp/fastp && \
#    chmod a+x ./fastp

# ENV PATH="${PATH}:/Gedi_1.0.5:/fastp:/STAR-2.7.10b/bin/Linux_x86_64_static"
#ENV PATH="${PATH}:/Gedi_1.0.5:/fastp:/STAR-2.7.4a/bin/Linux_x86_64_static"
ENV PATH="${PATH}:/Gedi_1.0.5"

WORKDIR /usr/src/

# RUN aws s3 cp --recursive s3://velia-analyses-dev/unit_test/input/ ./input
# RUN aws s3 cp --recursive s3://velia-analyses-dev/unit_test/output/ ./output
COPY ./src/ ./src/

# CMD [ "python", "./hello_world.py", "--message", "test_teXt", "--uppercase"]
