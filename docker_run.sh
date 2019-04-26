docker run -d -v $1/data:/data -v $1/results:/results sbc-ngs \
pathway \
/results \
/data \
$4 \
$5 \
$6 \
0 \
data/sample/SBC003382.fasta