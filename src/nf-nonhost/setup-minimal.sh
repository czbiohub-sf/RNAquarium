#!/bin/sh

DESTDIR="$PWD/"
bindir="bin"

mkdir -p ${DESTDIR}${bindir}

wget https://github.com/yttria-aniseia/fastq-lengths/releases/download/v0.1.1/fastq-lengths-0.1.1-x86_64.gz && \
	gzip -dc fastq-lengths-0.1.1-x86_64.gz > ${DESTDIR}${bindir}/fastq-lengths && chmod a+x ${DESTDIR}${bindir}/fastq-lengths && \
	rm fastq-lengths-0.1.1-x86_64.gz

wget https://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz -O PriceSource140408.tar.gz && \
	tar -xzvf PriceSource140408.tar.gz && cd PriceSource140408 && make && mv PriceSeqFilter ${DESTDIR}${bindir}/PriceSeqFilter && cd .. && \
	rm -rf PriceSource140408/ PriceSource140408.tar.gz

wget https://github.com/chanzuckerberg/czid-dedup/releases/download/v0.1.2/czid-dedup-Linux -O ${DESTDIR}${bindir}/czid-dedup && \
	chmod a+x ${DESTDIR}${bindir}/czid-dedup