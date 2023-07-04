from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import pandas as pd
from pathlib import Path
import subprocess
import sys
from tqdm import tqdm


seq_base_dir = Path('/hpc/projects/theory_ds/internship/jacob.paras/Gsnap_out')
results_file = 'counts.csv'
mate1_file = 'Unmapped.out.mate1.filteredbyBT.dedup.gsnapFiltered.fastq.gz'
mate2_file = 'Unmapped.out.mate2.filteredbyBT.dedup.gsnapFiltered.fastq.gz'
num_cores = 100


class MyParser(ArgumentParser):
    def error(self, message):
        sys.stderr.write(f'error: {message}\n')
        self.print_help()
        sys.exit(2)


def parse_args():
    parser = MyParser(description='This script parses the STAR_out directory and returns failed runs in order to re-submit them')
    parser.add_argument('--dir', default='', type=str, help='directory you want to save the output file in', metavar='')
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return vars(config)


# MAKE SURE THESE PATHS ARE CORRECT
config_dict = parse_args()
working_dir = Path(os.path.realpath(config_dict['dir'].rstrip('/').rstrip('\\')))


def get_count(gzip_file):
    gzip = subprocess.run(['gunzip', '-c', f'{gzip_file}'], check=True, capture_output=True)
    total_count = subprocess.run(['wc', '-l'], input=gzip.stdout, capture_output=True)
    return int(total_count.stdout.decode('utf-8').strip()) // 4


def return_count(accession):
    mate1 = seq_base_dir.joinpath(accession, mate1_file)
    mate2 = seq_base_dir.joinpath(accession, mate2_file)
    mate1_count = get_count(mate1)
    if not mate2.is_file():
        return accession, mate1_count
    mate2_count = get_count(mate2)
    return accession, mate1_count + mate2_count


def main():
    accessions = [subdir for subdir in os.listdir(seq_base_dir)]
    with tqdm(total=len(accessions)) as progress:
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            results_list = []
            results = [executor.submit(return_count, accession) for accession in accessions]
            for future in as_completed(results):
                result = future.result()
                results_list.append(result)
                progress.update()
    results_df = pd.DataFrame(results_list, columns=['SRA_run_accession', 'count'])
    results_df.to_csv(results_file, index=False)


if __name__ == '__main__':
    main()
