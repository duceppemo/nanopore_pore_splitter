# nanopore_pore_splitter

## Installation

TODO

## Usage
```
usage: python nano_pore_splitter.py [-h] -i /path/to/folder/with/fastq -o
                             /path/to/output/folder [-t 16] [-p 16] -r 0-256

Split Nanopore reads based on channel number.

optional arguments:
  -h, --help            show this help message and exit
  -i /path/to/folder/with/fastq, --input /path/to/folder/with/fastq
                        Input folder with fastq file(s), gzipped or not.
  -o /path/to/output/folder, --output /path/to/output/folder
                        Output folder.
  -t 16, --threads 16   Number of CPU. Default 16
  -p 16, --parallel 16  Number of samples to process in parallel. Default 4.
  -r 0-256, --range 0-256
                        Pore range to keep. Range between 0-512
```
