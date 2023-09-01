FROM debian:bookworm-20230612-slim
MAINTAINER yttria.aniseia@czbiohub.org

#WORKDIR /home

# install tools required for build
ENV PACKAGES build-essential gzip python3=3.11.* python3-biopython python3-pip wget
RUN apt-get update && apt-get install -y ${PACKAGES} && apt-get clean

# TODO: combine the tool installs for each step
# step 1 (fastq-dump)
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz -O sratoolkit.3.0.0-ubuntu64.tar.gz && \
	tar -xzvf sratoolkit.3.0.0-ubuntu64.tar.gz && mv sratoolkit.3.0.0-ubuntu64/bin/* /usr/bin/ && \
	rm -rf sratoolkit.3.0.0-ubuntu64 sratoolkit.3.0.0-ubuntu64.tar.gz

# step 2 (trimming and qc)
RUN	apt-get install -y fastp=0.23.2* && apt-get clean && \
	wget https://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz -O PriceSource140408.tar.gz && \
	tar -xzvf PriceSource140408.tar.gz && cd PriceSource140408 && make && mv PriceSeqFilter /usr/bin/PriceSeqFilter && cd .. && \
	rm -rf PriceSource140408/ PriceSource140408.tar.gz

# step 3 (map to zf)
RUN	apt-get install -y rna-star=2.7.10* samtools=1.16.1* && apt-get clean
# htseq setup from source is troublesome & requires git, pip, cython3, and swig, ends up better to just use pip
RUN pip install 'HTSeq==2.0.2' --break-system-packages

# configure (for sra toolkit)
# (run vdb-config and quit and ignore all errors and output)
RUN set +e; yes "q" | vdb-config -i > /dev/null 2>&1; set -e

# TODO: where should this go?
COPY util/median.py median.py
