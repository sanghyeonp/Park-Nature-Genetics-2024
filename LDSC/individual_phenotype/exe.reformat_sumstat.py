import code
# code.interact(local=dict(globals(), **locals()))
import csv
import argparse
import os
import random
from tqdm import tqdm
import subprocess
from util import print_args_info


def parse_args():
    parser = argparse.ArgumentParser(description="Script to reformat any summary statistics by defining column(s).")

    parser.add_argument('--sumstat', required=True,
                        help='Specify path to summary statistics.')
    parser.add_argument('--in_delimiter', required=False, default='comma',
                        help="Specify the delimiter used for input summary statistics. (Default=comma) (Choices=['comma', 'tab'])")
    parser.add_argument('--col_selection', required=True, nargs='+', 
                        help='Specify column(s) to use.')
    parser.add_argument('--col_rename', required=True, nargs='+',
                        help="Specify new names for the selected columns. It's important that the name(s) specified in the same order as selected column(s).")
    parser.add_argument('--out_filename', required=False, type=str, default=None,
                        help="Specify reformatted summary statistics filename. With default, it will have the following filename, <sumstat filename>.REFORMAT<random number>.txt. (Default=None)")
    parser.add_argument('--suffix', required=False, type=str, default='',
                        help="Specify to decorate the end of filename output. (Default='')")
    parser.add_argument('--prefix', required=False, type=str, default='',
                        help="Specify to decorate at the beginning of the filename output. (Default='')")
    parser.add_argument('--out_file_extension', required=False, type=str, default=None,
                        help="Specify file extension for the reformatted summary statistics. With default, it will be saved as .txt. (Defualt=None) (Choices=['txt', 'csv'])")
    parser.add_argument('--out_delimiter', required=False, default='comma',
                        help="Specify the delimiter used for output summary statistics. (Default=comma) (Choices=['comma', 'tab', 'whitespace'])")
    parser.add_argument('--outdir', required=False, default=None,
                        help="Directory where reformatted summary statistics will be saved. With default, it will be saved at the same location as input summary statistics. (Default=None)")
    parser.add_argument('--gunzip', action='store_true',
                        help='Specify to gunzip the reformatted summary statistics (Default=False).')
    parser.add_argument('--keep-only-gz', dest='keep_only_gz', action='store_true',
                        help='Specify to keep only gunzip file and delete other files generated (i.e., tab-delimited txt file) (Default=False).')
    parser.add_argument('--verbose', action='store_true', 
                        help='Specify to be verbose in printing (Default=False).')
    

    args = parser.parse_args()
    return args, parser


def run_bash(bash_cmd: str) -> list:
    """
    Run bash command.
    Return a list containing standard output, line by line.  
    """
    popen = subprocess.Popen(bash_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, _ = popen.communicate()
    return str(stdout, 'utf-8').strip().split('\n')


def pprint(msg):
    global verbose
    if verbose:
        print(msg)

def main(sumstat, in_delimiter, col_selection, col_rename, out_delimiter, save2):
    with open(sumstat, 'r') as f1, open(save2, 'wt', newline='') as f2:
        reader = csv.reader(f1, delimiter=in_delimiter)
        writer = csv.writer(f2, delimiter=out_delimiter)
        header = reader.__next__()
        # Check if specified delimiter is proper one
        assert len(header) > 1, "Specified delimiter for the input summary statistics seems wrong. Check `--in_delimiter`."

        # Check specified columns present in the header
        for col in col_selection:
            assert col in header, "Specified column [{}] is not present in the summary statistics. Check `--col_selection`.".format(col)    

        writer.writerow(col_rename)

        positions = [header.index(col) for col in col_selection]
        for row in tqdm(reader, leave=False, desc="Generating summary statistics", bar_format='{l_bar}{bar:20}{r_bar}{bar:-20b}'):
            new_row = [row[idx] if row[idx] else '.' for idx in positions ]
            writer.writerow(new_row)



def gunzip(sumstat):
    gz_filepath = sumstat + '.gz'
    cmd = "gzip -c {} > {}".format(sumstat, gz_filepath)
    run_bash(cmd)
    return gz_filepath


def delete_file(filename):
    run_bash('rm {}'.format(filename))


if __name__ == '__main__':
    args, parse_out = parse_args()
    verbose = args.verbose
    if verbose:
        print_args_info(args=args, parse_description=parse_out.description)

    # Set delimiter
    delimiters = {'comma':',', 'tab':'\t', 'whitespace':' '}
    assert args.in_delimiter in delimiters.keys(), "`--in_delimiter` should be one of these: {}".format(list(delimiters.keys()))
    assert args.out_delimiter in delimiters.keys(), "`--out_delimiter` should be one of these: {}".format(list(delimiters.keys()))
    args.in_delimiter = delimiters[args.in_delimiter]
    args.out_delimiter = delimiters[args.out_delimiter]

    # Check if the number of columns selected and the number of names specified match
    assert len(args.col_selection) == len(args.col_rename), "The number of columns selected and number of column names specified are different. Please check `--col_selection` and `col_rename`."

    # Format output directory
    if not args.outdir:
        args.outdir = os.path.split(args.sumstat)[0]

    # Set-up output file name
    if not args.out_filename:
        filename = '.'.join(os.path.split(args.sumstat)[-1].split(sep='.')[:-1])
        args.out_filename = args.prefix + filename + args.suffix
        # if not args.suffix:
        #     args.out_filename = filename + '.REFORMAT{}'.format(random.randint(1, 100))
        # else: # Add decorator to the filename if given
        #     args.out_filename = filename + args.suffix

    # Remove extension if it is provided in the output filename
    if args.out_filename.split(sep='.')[-1] in ['csv', 'tsv', 'txt']:
        args.out_filename = args.out_filename.split(sep='.')[0]

    # Add extension to the output file name
    if args.out_file_extension:
        if '.' in args.out_file_extension:
            args.out_filename = args.out_filename + args.out_file_extension
        else:
            args.out_filename = args.out_filename + '.' + args.out_file_extension
    else:
        args.out_filename = args.out_filename + '.txt'

    save2 = os.path.join(args.outdir, args.out_filename)

    pprint("Reformatted summary statistics will be saved at: {}".format(save2))

    main(sumstat=args.sumstat, 
        in_delimiter=args.in_delimiter, 
        col_selection=args.col_selection, 
        col_rename=args.col_rename, 
        out_delimiter=args.out_delimiter, 
        save2=save2
        )
    
    if args.gunzip:
        gz_filepath = gunzip(sumstat=save2)
        pprint("Reformatted summary statistics gz saved at: {}".format(gz_filepath))
    
        if args.keep_only_gz:
            delete_file(save2)
            pprint("Deleted {}".format(save2))