#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.myExecutor = 'slurm'
params.prefix = 'bin/'

workflow {
	install_fastq_lengths()
	install_fastq_namefilter()
	install_priceseqfilter(Channel.fromPath( 'pricesource.patch' ))
	install_czid_dedup()
	install_gsnap(Channel.fromList( ["", "sse42"]) )
}

process install_fastq_lengths {
	executor params.myExecutor
	cpus 1
	memory '2GB'
	time '10min'
	output: file("fastq-lengths")
	publishDir params.prefix, mode: 'move'

	script:
	def BINNAME="fastq-lengths"
	def URL="https://github.com/yttria-aniseia/fastq-lengths/archive/refs/tags/v0.2.8.tar.gz"
	"""
wget $URL -O fastq-lengths-0.2.8.tar.gz
tar -xzf fastq-lengths-0.2.8.tar.gz && cd fastq-lengths-0.2.8
make -j${task.cpus} && mv bin/${BINNAME} ../${BINNAME} && cd ..
rm -rf fastq-lengths-0.2.8/ && rm fastq-lengths-0.2.8.tar.gz
	"""
}

process install_fastq_namefilter {
	executor params.myExecutor
	cpus 1
	memory '2GB'
	time '10min'
	output: file("fastq-namefilter")
	publishDir params.prefix, mode: 'move'

	script:
	def BINNAME="fastq-namefilter"
	def URL="https://github.com/yttria-aniseia/fastq-namefilter/archive/refs/tags/0.2.8.tar.gz"
	"""
wget $URL -O fastq-namefilter-0.2.8.tar.gz
tar -xzf fastq-namefilter-0.2.8.tar.gz && cd fastq-namefilter-0.2.8
make -j${task.cpus} && mv bin/${BINNAME} ../${BINNAME} && cd ..
rm -rf fastq-namefilter-0.2.8/ && rm fastq-namefilter-0.2.8.tar.gz
	"""
}

process install_priceseqfilter {
	executor params.myExecutor
	cpus 1
	memory '2GB'
	time '15min'
	input: file("pricesource.patch")
	stageInMode 'copy'
	output: file("PriceSeqFilter")
	publishDir params.prefix, mode: 'move'

	script:
	def BINNAME="PriceSeqFilter"
	def URL="https://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz"
	"""
wget $URL -O PriceSource140408.tar.gz
tar -xzf PriceSource140408.tar.gz --one-top-level=PriceSource140408
patch -ruN --verbose -d PriceSource140408 -i \$PWD/pricesource.patch
cd PriceSource140408 && make -j${task.cpus} && mv ${BINNAME} ../${BINNAME} && cd ..
rm -rf PriceSource140408 && rm PriceSource140408.tar.gz
	"""
}

process install_czid_dedup {
	executor params.myExecutor
	cpus 1
	memory '1GB'
	time '5min'
	output: file("czid-dedup")
	publishDir params.prefix, mode: 'move'

	script:
	def BINNAME="czid-dedup"
	def URL="https://github.com/chanzuckerberg/czid-dedup/releases/download/v0.1.2/czid-dedup-Linux"
	"""
wget $URL -O $BINNAME
chmod a+x $BINNAME
	"""
}

process install_gsnap {
	executor params.myExecutor
	cpus 4
	memory '4GB'
	time '20min'
	tag "$simd_level"
	input: val(simd_level)
	output: file("bin/*")
	publishDir params.prefix, mode: 'move'

	script:
	def SIMD=simd_level ? "--with-simd-level=$simd_level" : ""
	def URL="http://research-pub.gene.com/gmap/src/gmap-gsnap-2023-12-01.tar.gz"
	"""
mkdir -p bin
wget $URL -O gmap-gsnap-2023-12-01.tar.gz
tar -xzf gmap-gsnap-2023-12-01.tar.gz
cd gmap-2023-12-01
./configure $SIMD --prefix=${PWD} && make -j${task.cpus}
find util -maxdepth 1 -executable -type f -exec mv {} ../bin/ \\; && cd ..
rm -rf gmap-2023-12-01 && rm gmap-gsnap-2023-12-01.tar.gz
	"""
}
