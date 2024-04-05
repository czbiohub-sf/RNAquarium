import csv
from collections import defaultdict
from dataclasses import dataclass
import json
from pathlib import Path

RNAQ_PROTO_DIR = Path(
    "/hpc/projects/balla_group/sra_experiments/RNAquarium_prototyping"
)
RNAQ_MAPPING_TSV = RNAQ_PROTO_DIR / "RNaquarium_run_to_bioproject_mapping.tsv"
OUTFILE = "test/sample_mapping.json"
RNAQ_SAMPLE_DIR = RNAQ_PROTO_DIR / "rnaquarium_unmapped_sample"


@dataclass
class MappingEntry:
    accession: str
    bioproject: str
    lib_layout: str
    avg_length: int
    size_mb: int


def main():
    mapping_dict = defaultdict(list)
    existing_runs = [run.stem for run in RNAQ_SAMPLE_DIR.glob("*")]

    with open(RNAQ_MAPPING_TSV, "r") as f:
        mapping = csv.reader(f, delimiter="\t")
        # Skip first row
        next(mapping)

        for line in mapping:
            entry = MappingEntry(*line)
            if entry.accession not in existing_runs:
                continue
            mapping_dict[entry.bioproject].append(entry.accession)

    with open(OUTFILE, "w") as f:
        json.dump(mapping_dict, f)


if __name__ == "__main__":
    main()
