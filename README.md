# longread_demux
This tool for demultiplexing long-read FASTQ files with Illumina unique dual index.

# Usage
```bash
longread_demux.pl -i input.fastq.gz -o output_prefix -d 70 -l index_list.txt
```

# Options

-h, --help: Show this help message and exit.

-o, --out: Output file prefix.

-i, --input: Input FASTQ file (supports gzipped files).

-d, --distance: Distance from the read ends to search for the index sequence.

-l, --list: A text file containing a list of index sequences to be used for demultiplexing. Please refer to the file in the index_list_example directory for an example of index sequence list.

