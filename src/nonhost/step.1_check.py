from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools
import json
import os
import subprocess as sp
import sys
from tqdm import tqdm


star_dir = 'STAR_out/'
pair_dir = 'PE'
single_dir = 'SE'
counts_dir = 'counts'
data_dir = 'data'
accessions_file = 'SRA_accession_list.1.27.23.txt'
failed_file = 'SRA_failed_runs.txt'
results_file = 'results.json'
log_file = 'Log.final.out'
mate1_file = 'Unmapped.out.mate1.gz'
mate2_file = 'Unmapped.out.mate2.gz'
htseq_count_file = 'htseq-count.txt'
min_file_size = 10
num_files = 54_189
num_cores = 60


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

# MAKE SURE THESE PATHS ARE CORRECT
config_dict = parse_args()
working_dir = os.path.realpath(config_dict['dir'].rstrip('/').rstrip('\\'))
star_path = os.path.join(working_dir, star_dir)
accessions_path = os.path.join(working_dir, data_dir, accessions_file)
failed_path = os.path.join(working_dir, data_dir, failed_file)
results_path = os.path.join(working_dir, data_dir, results_file)


def check_run(accession):
    PE = False
    SE = False
    nonexistent = False
    size_fail = False
    success = True
    error_message = ''
    counts_path = os.path.join(star_path, counts_dir, accession)
    pair_path = os.path.join(star_path, pair_dir, accession)
    single_path = os.path.join(star_path, single_dir, accession)

    # creating paths using os.path.join
    counts_files = itertools.starmap(os.path.join, [(counts_path, log_file), (counts_path, htseq_count_file)])
    PE_files = itertools.starmap(os.path.join, [(pair_path, log_file), (pair_path, mate1_file), (pair_path, mate2_file)])
    SE_files = itertools.starmap(os.path.join, [(single_path, log_file), (single_path, mate1_file)])

    # creating a mask over the list of paths to see whether or not the files exist
    counts_files_existence_mask = map(os.path.isfile, counts_files)
    PE_files_existence_mask = map(os.path.isfile, PE_files)
    SE_files_existence_mask = map(os.path.isfile, SE_files)

    # checking counts directory files
    if not all(counts_files_existence_mask):
        nonexistent = True
    # checking PE files
    if any(PE_files_existence_mask):
        PE = True
        if not all(PE_files_existence_mask):
            nonexistent = True
    # checking SE files
    elif not all(SE_files_existence_mask):
        nonexistent = True

    # handling nonexistent files
    # flip the values of the masks to get all the files that didn't pass
    if nonexistent:
        success = False
        if PE:
            nonexistent_files = itertools.compress(data=itertools.chain(counts_files, PE_files),
                                                   selectors=(not f for f in itertools.chain(counts_files_existence_mask, PE_files_existence_mask)))
        else:
            nonexistent_files = itertools.compress(data=itertools.chain(counts_files, SE_files),
                                                   selectors=(not f for f in itertools.chain(counts_files_existence_mask, SE_files_existence_mask)))
        error_message = f'Nonexistent files: {*nonexistent_files,}'
        return {'id': accession, 'success': success, 'msg': error_message}

    if PE:
        # TODO: make sure htseq files are created properly in the pipeline
        PE_size_check = itertools.starmap(os.path.join, [(pair_path, mate1_file), (pair_path, mate2_file)])
        # PE_size_check = itertools.starmap(os.path.join, [(counts_path, htseq_count_file), (pair_path, mate1_file), (pair_path, mate2_file)])
        # mask to check if file size is larger than an arbitrarily small, but nonzero, size
        PE_files_size_mask = map(lambda f: os.path.getsize(f) < min_file_size, PE_size_check)
        # empty compression means files passed and will evaluate to False
        if failed_files := itertools.compress(data=PE_size_check, selectors=PE_files_size_mask):
            size_fail = True
            error_message += f'Files failed size check: {*failed_files,}\n'

        zipped_files = itertools.starmap(os.path.join, [(pair_path, mate1_file), (pair_path, mate2_file)])
        # checking if .gzip files are properly zipped
        PE_zip_mask = map(lambda f: sp.getstatusoutput(f'gzip -t {f}')[0], zipped_files)
        # empty compression means files passed and will evaluate to False
        if not_compressed := itertools.compress(data=zipped_files, selectors=PE_zip_mask):
            size_fail = True
            error_message += f'Files not compressed properly: {*not_compressed,}'
    # SE
    else:
        # TODO: make sure htseq files are created properly in the pipeline
        SE_size_check = [os.path.join(single_path, mate1_file)]
        # SE_size_check = itertools.starmap(os.path.join, [(counts_path, htseq_count_file), (single_path, mate1_file)])
        SE_files_size_mask = map(lambda file: os.path.getsize(file) < min_file_size, SE_size_check)
        if failed_files := list(itertools.compress(data=SE_size_check, selectors=SE_files_size_mask)):
            size_fail = True
            error_message += f'Files failed size check: {*failed_files,}\n'
        zipfile = os.path.join(single_path, mate1_file)
        # 1 on failure which evaluates to True, 0 on success which evaluates to False
        if sp.getstatusoutput(f'gzip -t {zipfile}')[0]:
            size_fail = True
            error_message += f'Files not compressed properly: {zipfile}'

    success = False if size_fail else True
    # error_message = None if success else error_message

    return {'id': accession, 'success': success, 'msg': error_message}


def main():
    with open(accessions_path, 'r') as f:
        accessions = [line.strip() for line in f]
    with tqdm(total=len(accessions)) as progress:
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            results_json = []
            results = [executor.submit(check_run, accession) for accession in accessions]
            with open(failed_path, 'w') as f:
                for future in as_completed(results):
                    result = future.result()
                    if not result['success']:
                        f.write(f'{result["id"]}\n{result}\n')
                    results_json.append(result)
                    progress.update()
    with open(results_path, 'w') as f:
        json.dump(results_json, f)


if __name__ == '__main__':
    main()
