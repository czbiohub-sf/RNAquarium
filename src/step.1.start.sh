#!/bin/bash

SLURM_OUTDIR="slurm.out"
WorkingDir="/hpc/scratch/group.theory/jparas/zf_pipeline"
OUTDIR="${WorkingDir}/STAR_out"
log0="${WorkingDir}/logs/STAR.skip.log"
log1="${WorkingDir}/logs/STAR.process.log"
log2="${WorkingDir}/logs/SE_PE.log"

 
echo "making output directory $OUTDIR"
mkdir $OUTDIR
mkdir $OUTDIR/PE
mkdir $OUTDIR/SE
mkdir $OUTDIR/counts
mkdir $OUTDIR/other_stdout_stderr


echo "removing folder $SLURM_OUTDIR"
rm -rf $SLURM_OUTDIR
mkdir $SLURM_OUTDIR

echo "removing folder fastq"
rm -rf fastq
mkdir fastq

echo "removing folder prefetched"
rm -rf prefetched
mkdir prefetched

echo "checking logs directory"
if [ ! -d "${WorkingDir}/logs" ]
    then 
        mkdir "${WorkingDir}/logs"
fi

echo "remvoing previous log files"

if [ -e $log0 ]
    then 
        rm $log0
fi
if [ -e $log1 ]
    then 
        rm $log1
fi
if [ -e $log2 ]
    then 
        rm $log2
fi

sbatch step.1.fastqDump_trim_QCfilter_STAR.sh
