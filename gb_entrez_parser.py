import sys
import os
import time
from config import config

try:
    import argparse
    import urllib.error
    import codecs
    from Bio import SeqIO
    from Bio import Entrez
    from fuzzywuzzy import fuzz
    from fuzzywuzzy import process
    from openpyxl import Workbook
except ImportError as err:
    print("Import error!: " + err.name)
    exit(1)


Entrez.email = "bob@hotmail.com"


# Parsing for standard input arguments that user might enter for modularity
parser = argparse.ArgumentParser()
parser.add_argument('accessions')
parser.add_argument("--accessions", help="File Containing list of Accession IDs")
parser.add_argument("--recognition", help="Custom Recognition List")
parser.add_argument("--feature", help="Custom Feature List")
parser.add_argument("--header", help="Custom Headers List and data to pull out")
args = parser.parse_args()

def check_standard_input_files(stdin_file, stdin_list):
    if stdin_file is not None:
        try:
            lines = open(stdin_file).readlines()
            for line in lines:
                line = line.rstrip('\n')
                line = line.strip()
                stdin_list.append(line)
        except FileNotFoundError:
            print("Error: File was not found: " + stdin_file)
            exit(1)
    else:
        return 1


def grouper(iterable, n):
    for i in range(0, len(iterable), n):
        yield iterable[i:i + n]


def load_execute():
    """
        Function that loads all user arguments from console and processes/assigns them to variables
        setting up the environment for data execution

        To run the program, only the input file name is required. However other arguments can be input to add modularity
        """
    # If no arguments passed, exit
    if not len(sys.argv) > 1:
        print("Error: No Arguments passed")
        exit(1)

    # list of accessions that need be processed
    if len(sys.argv) == 1:
        input_file_name = sys.argv[1]
    else:
        input_file_name = args.accessions

    if not os.path.isfile(input_file_name):
        print("File could not be found: " + input_file_name)

    # the custom recognition files
    recognition_file = args.recognition
    feature_file = args.feature
    header_file = args.header

    # list to store input data
    recognition_list = []
    feature_list = []
    header_list = []

    ret = check_standard_input_files(recognition_file, recognition_list)
    if ret == 1 or (len(recognition_list) < 1):
        recognition_list = config.default_rec_list

    ret = check_standard_input_files(feature_file, feature_list)
    if ret == 1 or (len(feature_list) < 1):
        feature_list = config.default_fet_list

    ret = check_standard_input_files(header_file, header_list)
    if ret == 1 or (len(header_list) < 1):
        header_list = config.default_head_list

    input_file_name_edit = str(os.path.basename(input_file_name).split(".")[0])
    fasta_file = str("Result_" + input_file_name_edit + ".fasta")
    tsv_file = str("tsv_" + input_file_name_edit + ".tsv")
    log_file = str("log_" + input_file_name_edit + ".txt")

    execute(input_file_name, fasta_file, tsv_file, log_file, header_list, feature_list, recognition_list)



def execute(input_file_name,        # Basically list of accession IDs
            fasta_file,             # Name of the Fasta file the user wants to set
            tsv_file,               # Name of the tsv file the user wants to save data in
            log_file,               # Name of the log file the user wants to produce
            header_list,            # Header list
            feature_list,           # Feature list to look for "Hooks" to extract Sequences
            recognition_list):      # Recognition "Hook" list
    counter = 0
    long_delay = 0

 
    with codecs.open(input_file_name, 'r') as infile:
        chunks = grouper(infile.read().split("\n"), 300)
        for chunk in chunks:
            print("Processed " + str(counter) + " records...")
            while "" in chunk:
                chunk.remove('')
            acc_list = ",".join(chunk)
            try:
                fetched_gb = Entrez.efetch(db="nucleotide", id=acc_list, rettype="gbwithparts", retmode="text")
                for index, rec in enumerate(SeqIO.parse(fetched_gb, "gb")):
                    print("I was here!")


            except urllib.error.HTTPError:
                if urllib.error.HTTPError.code == 429:
                    time.sleep(5)
                    pass
                else:
                    a = False
            long_delay = long_delay + 1
            if long_delay == 20:
                time.sleep(20)
                long_delay = 0







if __name__ == '__main__':
    load_execute()
