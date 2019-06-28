docker run -d -v $1/data:/data -v $1/seqs:/seqs -v $1/results:/results sbc-ngs \
pathway \
/results \
/data \
$2 \
$3 \
8 \
0 \
/seqs/