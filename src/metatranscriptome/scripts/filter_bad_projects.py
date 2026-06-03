from dataclasses import dataclass
import json
from pathlib import Path
import re

HOME_DIR = Path(
    "/path/to/RNAquarium/merge_assemble/"
)
LOG_FILE = HOME_DIR / "slurm/nf-RNAquarium.out.bak"
FULL_MAP_FILE = HOME_DIR / "tmp/full_mapping.json"
FILT_MAP_FILE = HOME_DIR / "tmp/filt_mapping.json"

MERGE_FAIL_REGEX = re.compile(
    r"^\[[a-f0-9]{2}/[a-f0-9]{6}\] "
    r"NOTE: Process `MERGE_UNMAPPED \((PRJ[A-Z]+\d+)\)` terminated with an err"
    r"or exit status \((\d+)\) -- Error is ignored$"
)


@dataclass
class FailedMerge:
    bioproject_id: str
    exit_code: int


def main():
    content = LOG_FILE.read_text().splitlines()

    failing_bioprojects: list[FailedMerge] = []
    for line in content:
        reg_match = MERGE_FAIL_REGEX.match(line)
        if reg_match is None:
            continue

        bioproject_id, exit_code = reg_match.groups()
        exit_code = int(exit_code)
        failing_bioprojects.append(FailedMerge(bioproject_id, exit_code))

    failing_ids = [x.bioproject_id for x in failing_bioprojects]

    with FULL_MAP_FILE.open() as f:
        full_map = json.load(f)

    filt_map = {k: v for k, v in full_map.items() if k not in failing_ids}
    with FILT_MAP_FILE.open("w") as f:
        json.dump(filt_map, f)


if __name__ == "__main__":
    main()
