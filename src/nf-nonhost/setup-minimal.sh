#!/bin/sh

DESTDIR="$PWD/"
bindir="bin"

mkdir -p ${DESTDIR}${bindir}

wget https://github.com/yttria-aniseia/fastq-lengths/archive/refs/tags/v0.2.8.tar.gz -O fastq-lengths-0.2.8.tar.gz && \
	tar -xzvf fastq-lengths-0.2.8.tar.gz && cd fastq-lengths-0.2.8 && make && mv bin/fastq-lengths ${DESTDIR}${bindir}/fastq-lengths && cd ../.. && \
	rm -rf fastq-lengths-0.2.8/ && rm fastq-lengths-0.2.8.tar.gz

wget https://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz -O PriceSource140408.tar.gz && \
	tar -xzvf PriceSource140408.tar.gz && cd PriceSource140408 && make && mv PriceSeqFilter ${DESTDIR}${bindir}/PriceSeqFilter && cd .. && \
	rm -rf PriceSource140408/ && rm PriceSource140408.tar.gz

wget https://github.com/chanzuckerberg/czid-dedup/releases/download/v0.1.2/czid-dedup-Linux -O ${DESTDIR}${bindir}/czid-dedup && \
	chmod a+x ${DESTDIR}${bindir}/czid-dedup