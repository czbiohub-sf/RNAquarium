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
 - BBMap=39.01
 - bowtie2
 - czid-deup=0.1.0
 - fastp=0.23.2
 - FastQC=0.12.1
 - gmap-gsnap=2021-12-17
 - HTSeq=2.0.2
 - Price=1.2
 - samtools=1.16.1
 - sratoolkit=3.0.0
 - STAR=2.7.10
 - Trimmomatic=0.38

it is recommended to run the pipeline using the docker/singularity container.  
e.g.:
```
TODO
```

Price can be obtained from https://derisilab.ucsf.edu/software/price/index.html
czid-dedup can be obtained from https://github.com/chanzuckerberg/czid-dedup/releases/
```
wget https://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz -O PriceSource140408.tar.gz
tar -xzvf PriceSource140408.tar.gz && cd PriceSource140408 && make && cd ..
export PATH=$PATH:$PWD/PriceSource140408
rm PriceSource140408/
```

### Create a working directory and input list
```
cd Zebrafish-RNAquarium/src/nonhost
mkdir test && cp data/SRA_accession_list.test.txt test/SRA_accession_list.test.txt
```

### Step 1 fastqdump
```
bash step.1.start.sh test
```
`test` is the name of our working directory  
output is in `test/fastq/{accession}/*.fastq.gz`

### Step 2 QC
```
bash step.2.start.sh test
```
output is in `test/fastq/{accession}/*.PRICEfiltered.fastq.gz`

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