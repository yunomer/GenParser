# Welcome to GenParser 👋
![Version](https://img.shields.io/badge/version-0.10-blue.svg?cacheSeconds=2592000)
[![Documentation](https://img.shields.io/badge/documentation-yes-brightgreen.svg)](https://github.com/yunomer/GenParser)

> A simple script to download genbank files and parse to extract information asked for in config files, and save them in xlsx format.
> The script can also generate fasta files and log files that contain Accession IDs that failed to produce valid data in xlsx file.

## 🚀 Usage
Make sure you have Python 3.x installed and pip version >= 9.0.1
   
First run the following command at the root of your project:
```sh
    pip3 install -r requirements.txt
```
This will ensure you have all the requirements installed on your system.

To run the script you can run the command:
```sh
    python3 gb_entrez_parser.py [ACCESSION LIST FILE]
```
Doing so will produce a xlsx file called genbankParsedData. Which will contain parsed data based in the header list you provide.
In the case above, no header list was provided and so it used a Generic Header list found in the config file. You can modify that
list directly if you plan to use the same headers over and over again.

The Available Arguments that can be passed into the script are:
```sh
    > Flags for Input Files:
    --header        [HEADER LIST FILE]
    --feature       [FEATURE LIST FILE]
    --recognition   [RECOGNITION LIST FILE]
    
    > Flags:
    -f              Produce a Fasta file for all sequences found
    -t              Convert and produce a TSV file of xlsx file produced
    -l              Produce a logs file containing all Accession IDs with no sequences
```

The main purpose of this script is to extract sequences related to Features in the Genbank files. 
To address that, currently users have to identify the sequence location in the Genbank file (Feature), and then adding a "Hook" or key to look for.
An example is provided in the config file (features.txt, recognition.txt).

In In the example, we see that the features list contains the name of the feature to look in (i.e CDS, Gene). The script requires the Feature to contain a "Hook" to latch onto and extract the sequence if present. The recognition hook can be either a key (Text before the '=') or the value (Text after the '=').
If no custom Recognition list or Feature list is provided, a default one is utilized, that can be modified in the config file. To get extract the sequence, in the custom header file use 'seq'.
However if you would like to use a custom Recognition and Feature list to identify specific sequences you can simply run the script using flags above:
```sh
    python3 gb_entrez_parser.py [ACCESSION LIST FILE] --feature [FEATURE FILE] -- recognition [RECOGNITION FILE]
```

To produce Fasta file, log file or TSV file simply add the flags at the end. 
    

## Author
👤 [@yunomer](https://github.com/yunomer)

## Show your support

Give a ⭐️ if this project helped you!