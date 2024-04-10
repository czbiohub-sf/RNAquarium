#!/usr/bin/env python

import argparse
from collections import defaultdict
import json
from pathlib import Path

SCRIPT_DESC = "Create mapping between BioProject IDs and accession IDs"
TAB_DESC = "SRA accessions tab file"
ACC_DESC = "SRA accessions to extract"
OUTFILE_DESC = "Path to output mapping JSON"
UNMAPPED_DESC = "Path to output with unmapped accessions"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=SCRIPT_DESC)
    parser.add_argument("--tab", type=str, help=TAB_DESC)
    parser.add_argument('--acc', type=str, help=ACC_DESC)
    parser.add_argument('--outfile', type=str, help=OUTFILE_DESC)
    parser.add_argument('--unmapped', type=str, help=UNMAPPED_DESC)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    tabfile = Path(args.tab)
    accfile = Path(args.acc)
    outfile = Path(args.outfile)
    unmapfile = Path(args.unmapped)

    mapping_dict = defaultdict(list)

    with open(accfile, "r") as f:
        accessions = set(f.read().splitlines())

    with open(tabfile, "r") as f:
        header = f.readline().strip().split("\t")
        bioproj_idx = header.index("BioProject")
        accession_idx = header.index("Accession")

        for row in f:
            row = row.strip().split("\t")
            this_bioproj_id = row[bioproj_idx]
            this_accession_id = row[accession_idx]

            if this_bioproj_id == "-":
                continue
            if this_accession_id not in accessions:
                continue

            mapping_dict[this_bioproj_id].append(this_accession_id)
            accessions.remove(this_accession_id)

            if not accessions:  # No more accessions to map
                break

    with open(outfile, "w") as f:
        json.dump(mapping_dict, f)

    with open(unmapfile, "w") as f:
        f.write("\n".join(accessions))


if __name__ == "__main__":
    main()
