# sbc-ngs

## From FASTA files

To run sbc-ngs from pre-compiled FASTA files of target sequences, run:

1. `bash docker_build.sh`
2. `bash docker_run.sh [FULL_PATH_TO_DATA_DIRECTORY] 1000 128 8`
(e.g. `bash docker_run.sh /Users/username/sbc-ngs/example/fasta 1000 128 8`

This uses example data provided here in the `example` directory. The required
format of this directory follows that of `example`. Specifically, the directory
must contain the following subdirectories:

1. `data`. This directory contains both the sequence data to be analysed
(typically in fastq format, but other formats also supported), along with a file
`barcodes.csv` which defines the barcode sequences applied to each sample and the
sequence id of each sample. This file requires the headers `well`, `known_seq_id`,
`forward` and `reverse`.

2. `seqs`. This directory contains fasta files of the template sequences against
which the sequence data will be aligned. The fasta files must share the same
names as the values in the `known_seq_id` column of the `barcodes.csv` file,
described above.


## From template data held in JBEI-ICE

(Note: this is the approach typically used in the SYNBIOCHEM centre.)

To run sbc-ngs from target sequences held in JBEI-ICE, run:

1. `bash docker_build.sh`
2. `bash docker_run_ice.sh [TARGET_DIRECTORY] [ICE_URL] [ICE_USERNAME] [ICE_PASSWORD] 1000 128 8`
(e.g. `bash docker_run_ice.sh /Users/username/sbc-ngs/example/ice https://ice.mylab.com user@mylab.com password 1000 128 8`