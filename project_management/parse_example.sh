cp scripts/project_management/parse_vcf.pl temp/test.pl 
chmod u+x temp/test.pl

./temp/test.pl 1 0 100000000 data/genotypes/chr1
