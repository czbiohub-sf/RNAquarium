from argparse import ArgumentParser
import os
import pandas as pd
from pathlib import Path
import sys
from tqdm import tqdm


data_dir = 'data'
accessions_file = 'SRA_accession_list.1.27.23.txt'
star_dir = 'STAR_out/'
pair_dir = 'PE'
single_dir = 'SE'
other_dir = 'other_stdout_stderr'
log_file = 'Log.final.out'
readlength_file = 'readlength.txt'
saved_file = 'parsed_reads.csv'
min_file_size = 10


class MyParser(ArgumentParser):
    def error(self, message):
        sys.stderr.write(f'error: {message}\n')
        self.print_help()
        sys.exit(2)


def parse_args():
    parser = MyParser(description='This script parses the STAR_out directory for number of reads and readlengths for each run')
    parser.add_argument('--dir', default='', type=str, help='path to the parent directory of the STAR output directory', metavar='')
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return vars(config)


# MAKE SURE THESE PATHS ARE CORRECT
config_dict = parse_args()
working_dir = Path(config_dict['dir'].rstrip('/').rstrip('\\')).resolve()
accessions_path = working_dir.joinpath(data_dir, accessions_file)
star_path = working_dir.joinpath(star_dir)


def main():
    with open(accessions_path, 'r') as f:
        accessions = [line.strip() for line in f]
    df_list = []
    with tqdm(total=len(accessions)) as progress:
        for accession in accessions:
            if not ((readlength_path := star_path.joinpath(other_dir, accession, readlength_file)).is_file() and
                    readlength_path.stat().st_size > min_file_size):
                continue
            if star_path.joinpath(pair_dir, accession).is_dir():
                endedness = 'PE'
            elif star_path.joinpath(single_dir, accession).is_dir():
                endedness = 'SE'
            else:
                continue
            if not ((log_path := working_dir.joinpath(star_dir, endedness, accession, log_file)).is_file() and
                    log_path.stat().st_size > min_file_size):
                continue
            log_df = pd.read_table(log_path, header=None, converters={0: lambda s: s.strip()})
            log_df = log_df.set_index(0)
            num_reads = log_df.loc['Number of input reads |']
            with open(readlength_path) as f:
                try:
                    readlength1, *readlength2 = f.readline().strip()
                except ValueError:
                    print(readlength_path)
                    return
            readlength2 = readlength2 or None
            df_list.append((accession, num_reads, readlength1, readlength2))
            progress.update()
    df = pd.DataFrame(df_list, columns=['Run ID', 'Number of reads', 'Mate1', 'Mate2'])
    df.to_csv(saved_file)


if __name__ == '__main__':
    main()
