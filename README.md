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

To run the script with custom headers, you can simply add the flag and run:
```sh
    python3 gb_entrez_parser.py [ACCESSION LIST FILE] --header [HEADER FILE]
```

    

## Author
👤 **Omer Ashfaq**

* Github: [@yunomer](https://github.com/yunomer)

## Show your support

Give a ⭐️ if this project helped you!