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
```
TODO:
    Conda env
    TOOLS directory
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
`test` is the name of our workint directory  
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