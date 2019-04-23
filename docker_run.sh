docker run -d -v $1/data:/data -v $1/results:/results sbc-ngs \
pathway \
/results \
/data \
https://ice.synbiochem.co.uk \
$2 \
$3 \
gaattcaaaagatcttttaagaag \
ttactcgagtttggatcc \
$4 \
$5 \
$6 \
True