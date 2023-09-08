# Non-zebrafish reads from  datasets in SRA

## Features  
- Input a list of SRA run accessions
- Output non-zebrafish reads as compressed fastq files
- QC reads and use three aligners to map/remove zebrafish reads
- Distributed computing on HPC environments (Slurm job scheduler)
- Highly automated


## Usage

### Clone the repository
```
git clone https://github.com/czbiohub-sf/Zebrafish-RNAquarium
```
### Install dependencies
The non-zebrafish reads pipeline requires the following dependencies in `PATH`:
 - BBMap=39.01 ( https://sourceforge.net/projects/bbmap/files/ )
 - bowtie2 ( https://sourceforge.net/projects/bowtie-bio/files/bowtie2/ )
 - czid-dedup=0.1.0 ( https://github.com/chanzuckerberg/czid-dedup )
 - fastp=0.23.2 ( https://github.com/OpenGene/fastp )
 - FastQC=0.12.1 ( https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ )
 - gmap-gsnap=2021-12-17 ( http://research-pub.gene.com/gmap/ )
 - HTSeq=2.0.2 ( https://htseq.readthedocs.io/en/master/install.html )
 - Price=1.2 ( https://derisilab.ucsf.edu/software/price/index.html )
 - samtools=1.16.1 ( https://www.htslib.org/download/ )
 - sratoolkit=3.0.0 ( https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit )
 - STAR=2.7.10 ( https://github.com/alexdobin/STAR )
 - Trimmomatic=0.38 ( https://github.com/usadellab/Trimmomatic )

Alternatively, most pre-installation can be skipped by running the nextflow
pipeline with either the `-with-docker` or `-with-singularity` parameters,
which will automatically pull and use containers for most steps.  
in this case, currently `fastq-lengths`, `PRICE`, and `czid-dedup` are not
available as containers.
The script `setup-minimal.sh` will download these.
e.g.:
```
./setup-minimal.sh
nextflow run -with-singularity
```


### Create a working directory and input list
```
cd Zebrafish-RNAquarium/src/nonhost
mkdir test && cp data/SRA_accession_list.test.txt SRA_accession_list.test.txt
```

### Step 1 fastqdump
```
nextflow run step.1.nf --accessions_list SRA_accession_list.test.txt --parallel_downloads 10 -resume
```
output is in `fastq/{accession}/*.fastq.gz`

### Step 2 QC
```
nextflow run step.2.nf --accessions_list SRA_accession_list.test.txt -resume
```
output is in `fastq/{accession}/*.PRICEfiltered.fastq.gz`

### Step 3 STAR
```
bash step.3.start.sh test
```
Runs STAR twice, once optimized to isolate reads unmappable to zebrafish genome, a second time optimized to count reads mapped to zebrafish genes  
Nonhost outputs:  
- pair-ended: `test/STAR_out/PE/{accession}/Unmapped.out.mate1.gz`
- pair-ended: `test/STAR_out/PE/{accession}/Unmapped.out.mate2.gz`
- single-ended: `test/STAR_out/SE/{accession}/Unmapped.out.mate1.gz`  

Gene read count output:  
- single-ended: `test/STAR_out/counts/{accession}/htseq-count.txt`

### Step 4 Bowtie
```

```

### Step 5 dedup
```

```

### Step 6 gsnap
```

```