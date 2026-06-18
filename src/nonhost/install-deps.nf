#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.myExecutor = 'slurm'
params.prefix = 'bin/'

workflow {
	install_seq_detective()
	install_fastq_lengths()
	install_fastq_namefilter()
	//install_priceseqfilter(Channel.fromPath( 'pricesource.patch' ))
	install_czid_dedup()
	install_gsnap(Channel.fromList( ["", "sse42", "avx2"]) )
}

process install_seq_detective {
	executor params.myExecutor
	cpus 1
	memory '2GB'
	time '10min'
	output: file("seq-detective")
	publishDir params.prefix, mode: 'move'

	script:
	def BINNAME="seq-detective"
	def URL="https://github.com/czbiohub-sf/seq-tech-detective.git"
	"""
git clone $URL seq-tech-detective-core -b CORE --depth 1
cd seq-tech-detective-core
mamba env create -f environment.yml -n seq-detective
mamba activate seq-detective
make -j${task.cpus} && \
mv -u build/seq-detective/bin/* ../bin/ && \
mv -u build/seq-detective/libexec/ ../ && \
mv -u build/seq-detective/share/ ../ && \
cd ..
rm -rf seq-tech-detective-core/
	"""
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
	def BINNAME2="fastq-numfilter"
	def VER="0.3.0"
	def URL="https://github.com/yttria-aniseia/fastq-namefilter/archive/refs/tags/${VER}.tar.gz"
	"""
wget $URL -O fastq-namefilter-${VER}.tar.gz
tar -xzf fastq-namefilter-${VER}.tar.gz && cd fastq-namefilter-${VER}
make -j${task.cpus} && \
mv bin/${BINNAME} ../${BINNAME} && \
mv bin/${BINNAME2} ../${BINNAME2}
cd ..
rm -rf fastq-namefilter-${VER}/ && rm fastq-namefilter-${VER}.tar.gz
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
	output: file("pubbin/*")
	publishDir params.prefix, mode: 'move', saveAs: { it.replace(/pubbin\//, "") }

	script:
	def SIMD=simd_level ? "--with-simd-level=$simd_level" : ""
	def VER="2025-04-19" // "2021-12-17" // 2024-02-22
	def URL="http://research-pub.gene.com/gmap/src/gmap-gsnap-${VER}.tar.gz"
	"""
mkdir -p pubbin
wget $URL -O gmap-gsnap-${VER}.tar.gz
tar -xzf gmap-gsnap-${VER}.tar.gz
cd gmap-${VER}
./configure $SIMD --prefix=${PWD} && make -j${task.cpus}
find src -maxdepth 1 -executable -type f -exec mv {} ../pubbin/ \\;
find util -maxdepth 1 -executable -type f -exec mv {} ../pubbin/ \\;
cd ..
rm -rf gmap-${VER} && rm gmap-gsnap-${VER}.tar.gz
	"""
}

process install_hisat2 {
	executor params.myExecutor
	cpus 2
	memory '4GB'
	time '10min'
	output:
	file("hisat2")
	file("hisat2-build*")
	file("hisat2-align*")
	file("hisat2-inspect")
	file("hisat2*.py")
	publishDir params.prefix, mode: 'move'
	
	script:
	BINNAME="hisat2 hisat2-build* hisat2-align* hisat2-inspect hisat2*.py"
	def VER="fE9QCsX3NH4QwBi" // 2.2.1-source
	def URL="https://cloud.biohpc.swmed.edu/index.php/s/$VER/download"
	"""
wget $URL -O hisat2-${VER}.zip
unzip hisat2-${VER}.zip && cd hisat2-${VER}.zip
make -j${task.cpus} both EXTRA_FLAGS="-DPOPCNT_CAPABILITY -std=c++11 -mavx2 -g0" && mv $BINNAME ../
cd ..
rm -rf hisat2-${VER} && rm hisat2-${VER}.zip
	"""	
}
