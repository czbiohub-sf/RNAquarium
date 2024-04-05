#!/usr/bin/env python

import argparse
from pathlib import Path
import subprocess
import json

SCRIPT_DESC = "Merge FASTQs for a single BioProject"
BP_ID_DESC = "BioProject ID"
MAPPING_DESC = "Mapping of BioProject IDs to sequencing runs"
INDIR_DESC = "Directory containing sequencing runs"
OUTDIR_DESC = "Output directory for merged FASTQs"

OUTFQ_MATE1 = "Unmapped.out.mate1.filteredbyBT.dedup.merged.fastq"
OUTFQ_MATE2 = "Unmapped.out.mate2.filteredbyBT.dedup.merged.fastq"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=SCRIPT_DESC)
    parser.add_argument("--ID", type=str, help=BP_ID_DESC)
    parser.add_argument('--mapping', type=str, help=MAPPING_DESC)
    parser.add_argument('--indir', type=str, help=INDIR_DESC)
    parser.add_argument('--outdir', type=str, help=OUTDIR_DESC)
    args = parser.parse_args()
    return args


def get_run_fqs(fq_dir: Path) -> list[Path]:
    fqs = list(fq_dir.glob("*.fastq.gz"))
    return fqs


def append_fq_to_file(src_fq: Path, target_fq: Path) -> None:
    cmd = ["zcat", str(src_fq)]
    with open(target_fq, "a") as f:
        subprocess.run(cmd, stdout=f)


def gzip(fq: Path) -> None:
    cmd = ["gzip", str(fq)]
    subprocess.run(cmd)


def main():
    args = parse_args()
    bp_id = args.ID
    indir = Path(args.indir)
    outdir = Path(args.outdir)

    outdir.mkdir(exist_ok=True)
    fqs_to_compress = set()

    with open(args.mapping, "r") as f:
        mapping = json.load(f)

    run_ids = mapping[args.ID]
    for run_id in run_ids:
        run_fq_dir = indir / run_id
        run_fqs = get_run_fqs(run_fq_dir)

        endedness = len(run_fqs)
        if endedness == 1:
            bp_outdir = outdir / f"{bp_id}_S"
            bp_outdir.mkdir(exist_ok=True)
            bp_outfq = bp_outdir / OUTFQ_MATE1

            append_fq_to_file(run_fqs[0], bp_outfq)

            fqs_to_compress.add(bp_outfq)
        elif endedness == 2:
            bp_outdir = outdir / f"{bp_id}_P"
            bp_outdir.mkdir(exist_ok=True)

            infq_mate1 = next(fq for fq in run_fqs if "mate1" in fq.name)
            infq_mate2 = next(fq for fq in run_fqs if "mate2" in fq.name)

            outfq_mate1 = bp_outdir / OUTFQ_MATE1
            outfq_mate2 = bp_outdir / OUTFQ_MATE2

            append_fq_to_file(infq_mate1, outfq_mate1)
            append_fq_to_file(infq_mate2, outfq_mate2)

            fqs_to_compress.add(outfq_mate1)
            fqs_to_compress.add(outfq_mate2)
        else:
            raise ValueError(f"Unknown number of FASTQs found for {run_id}!")

    for fq in fqs_to_compress:
        gzip(fq)


if __name__ == "__main__":
    main()
