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



def execute(input_file_name, # Basically list of accession IDs
            fasta_file,
            tsv_file,
            log_file,
            header_list,
            feature_list,
            recognition_list):
    counter = 0
    long_delay = 0

    output_file = open(fasta_file, "w")
    # tsv file for data dump based on header list
    tsv_file = open(tsv_file, "w")
    # Log file to record Errors
    log_file = open(log_file, "w")
    for header in header_list:
        tsv_file.write(header + "\t")
    tsv_file.write("\n")

    with codecs.open(input_file_name, 'r') as infile:
        accs = grouper(infile.read().split("\n"), 300)
        for chunk in accs:
            print("Processed " + str(counter) + " records...")
            while "" in chunk:
                chunk.remove('')
            acc_list = ",".join(chunk)
            try:
                fetched_gb = Entrez.efetch(db="nucleotide", id=acc_list, rettype="gbwithparts", retmode="text")
                unique_list_error_nonetype = []
                unique_list_found = []
                complete_list = []
                for index, rec in enumerate(SeqIO.parse(fetched_gb, "gb")):
                    print(rec)
                    exit(1)
                    complete_list.append(rec.id)
                    for feature in rec.features:
                        ret = parser(feature, recognition_list, rec, feature_list)
                        if ret is None:
                            if rec.id not in unique_list_error_nonetype:
                                unique_list_error_nonetype.append(rec.id)
                            break
                        gene = ret[0]
                        sequence = ret[1]
                        if gene != "No appropriate gene found!":
                            output_file.write("> " + rec.id + "|" + str(gene) + "\n" + str(sequence) + "\n")
                            tsv_file.write(rec.id)
                            for header in header_list:
                                fetch_annotation_ret = fetch_annotation(header, rec)
                                if header is "sequence":
                                    fetch_annotation_ret = str(sequence)
                                elif header is "taxonomy":
                                    try:
                                        fetch_annotation_ret = rec.annotations[header]
                                    except IndexError:
                                        fetch_annotation_ret = ""
                                elif len(fetch_annotation_ret) > 0:
                                    fetch_annotation_ret = fetch_annotation_ret[0]
                                else:
                                    fetch_annotation_ret = ""

                                if fetch_annotation_ret is not None:
                                    tsv_file.write(str(fetch_annotation_ret) + "\t")
                                else:
                                    tsv_file.write("\t")
                            tsv_file.write("\n")
                            if rec.id not in unique_list_found:
                                unique_list_found.append(rec.id)
                            counter = counter + 1
                            break
                        else:
                            pass
                output_file.close()
                tsv_file.close()

                # Printing out the IDs with Nonetype error.
                for accession_id in unique_list_error_nonetype:
                    log_file.write(
                        accession_id + "\t" + "No Data found!" + "\n")

                for check in complete_list:
                    if check not in unique_list_found:
                        if check not in unique_list_error_nonetype:
                            log_file.write(
                                check + "\t" + "No Sequence found" + "\n")

                log_file.close()
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
