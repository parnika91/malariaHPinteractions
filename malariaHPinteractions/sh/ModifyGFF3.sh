#! /bin/bash/


# remove all comment lines
grep -v '^#.' Pcynomolgi.gff3 > temp && mv temp Pcynomolgi.gff3

# attach prefix
sed -e 's/^/Pcy_chr/' Pcynomolgi.fa > temp && mv temp Pcynomolgi.fa

