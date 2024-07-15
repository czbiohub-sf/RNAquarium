#!/bin/bash
#SBATCH --job-name=rnaquarium-deps
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

DESTDIR="$PWD/"
bindir="bin"

mkdir -p ${DESTDIR}${bindir}

wget https://github.com/yttria-aniseia/fastq-lengths/archive/refs/tags/v0.2.8.tar.gz -O fastq-lengths-0.2.8.tar.gz && \
	tar -xzvf fastq-lengths-0.2.8.tar.gz && cd fastq-lengths-0.2.8 && make && mv bin/fastq-lengths ${DESTDIR}${bindir}/fastq-lengths && cd .. && \
	rm -rf fastq-lengths-0.2.8/ && rm fastq-lengths-0.2.8.tar.gz

wget https://github.com/yttria-aniseia/fastq-namefilter/archive/refs/tags/0.2.8.tar.gz -O fastq-namefilter-0.2.8.tar.gz && \
	tar -xzvf fastq-namefilter-0.2.8.tar.gz && cd fastq-namefilter-0.2.8 && make && mv bin/fastq-namefilter ${DESTDIR}${bindir}/fastq-namefilter && cd .. && \
	rm -rf fastq-namefilter-0.2.8/ && rm fastq-namefilter-0.2.8.tar.gz

wget https://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz -O PriceSource140408.tar.gz && \
	tar -xzvf PriceSource140408.tar.gz && patch -ruN --verbose -d PriceSource140408 -i $PWD/pricesource.patch && \
	cd PriceSource140408 && make && mv PriceSeqFilter ${DESTDIR}${bindir}/PriceSeqFilter && cd .. && \
	rm -rf PriceSource140408/ && rm PriceSource140408.tar.gz
# diff -ruN bin/PriceSource140408_modified/ bin/PriceSource140408/ --exclude=.git --exclude=*.o --exclude PriceSeqFilter

wget https://github.com/chanzuckerberg/czid-dedup/releases/download/v0.1.2/czid-dedup-Linux -O ${DESTDIR}${bindir}/czid-dedup && \
	chmod a+x ${DESTDIR}${bindir}/czid-dedup

## replace "sse42" with an appropriate fallback, if necessary
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2023-12-01.tar.gz -O gmap-gsnap-2023-12-01.tar.gz && \
	tar -xzvf gmap-gsnap-2023-12-01.tar.gz && cd gmap-2023-12-01 && ./configure && make && \
	./configure --with-simd-level=sse42 --prefix=${DESTDIR} && make && make install && \
	find util -maxdepth 1 -executable -type f -exec mv {} ${DESTDIR}${bindir} \; && cd .. && \
	rm -rf gmap-2023-12-01 && rm gmap-gsnap-2023-12-01.tar.gz
