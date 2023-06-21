from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
import os
import sys


class MyParser(ArgumentParser):
    def error(self, message):
        sys.stderr.write(f'error: {message}\n')
        self.print_help()
        sys.exit(2)


def parse_args():
    parser = MyParser(description='This script parses the STAR_out directory and returns failed runs in order to re-submit them')
    parser.add_argument('--dir', default='', type=str, help='path to the parent directory of the STAR output directory', metavar='')
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return vars(config)


def check_run(run):
    pass


config_dict = parse_args()
working_dir = os.path.realpath(config_dict['dir'].rstrip('/').rstrip('\\'))

star_dir = 'STAR_out/'
accessions_file = 'SRA_accession_list.1.27.23.txt'

star_path = os.path.join(working_dir, star_dir)
accessions_path = os.path.join(working_dir, accessions_file)

num_files = 54189
num_cores = 60


def main():
    with open(accessions_file, 'r') as f:
        accessions = [line.strip() for line in f]
    with ProcessPoolExecutor(max_workers=num_cores) as ex:
        results = ex.map(check_run, accessions, chunksize=num_files // num_cores)


if __name__ == '__main__':
    main()
