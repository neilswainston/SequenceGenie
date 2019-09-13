# SequenceGenie

These algorithms are tuned to exploit multiple processors. Running the code with
multiple processors will result in a significant performance improvement.

Due to its use of Docker, this code is runnable (and has been tested) on cloud
computing platforms such as Google Compute Engine.

## From FASTA files

To run sbc-ngs from pre-compiled FASTA files of target sequences, run:

1. `bash docker_build.sh`
2. `bash docker_run.sh [FULL_PATH_TO_DATA_DIRECTORY] [MIN_SEQ_LEN] [MAX_SEQ_FILES]`
(e.g. `bash docker_run.sh /Users/username/sbc-ngs/example/fasta 1000 128`

The value `[MIN_SEQ_LEN]` corresponds to the minimum sequence length of a
data read to be considered in the analysis. Depending upon the length of the
template sequences to be matched to, 1000 is a sensible default.

The value `[MAX_SEQ_FILES]` corresponds to the number of data sequence files
to consider in the analysis. Note that a data sequence files (e.g. a fastq or
fasta file) is likely to contain multiple reads. This value corresponds to the
number of files to read, *not* the number of sequences. Increasing this value
(considering more data) will lead to more reliable results at the expense of
performance.

This uses example data provided here in the `example/fasta` directory. The required
format of this directory follows that of `example/fasta`. Specifically, the directory
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
2. `bash docker_run_ice.sh [TARGET_DIRECTORY] [ICE_USERNAME] [ICE_PASSWORD] [MIN_SEQ_LEN] [MAX_SEQ_FILES]`
(e.g. `bash docker_run_ice.sh /Users/username/sbc-ngs/example/ice user@mylab.com password 1000 128`
