cp scripts/project_management/parse_vcf.pl temp/test.pl 
chmod u+x temp/test.pl
chmod u+x scripts/project_management/parse_vcf.pl

grep full support/target_regions_filtered.tab | while read chr start end code gene; do
    
    script=cluster/submission/pro${code}.sh
    
    echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd

export PERL5LIB=\$PERL5LIB:/software/additional/lib/perl5/x86_64-linux-thread-multi:/software/additional/lib/perl5:/software/additional/lib/perl5/cpan;
export TEMP=/data_n2/vplagnol/Projects/fluidigm/association_v2/temp/Rtemp/

/usr/bin/perl ./scripts/project_management/parse_vcf.pl $chr $start $end data/genotypes/$code
#./temp/test.pl $chr $start $end data/genotypes/$code

#R CMD BATCH --no-save --no-restore --code=data/genotypes/${code} scripts/project_management/step2_process_genotypes.R cluster/R/pro_${code}.out

" > $script
    
    echo $script    
    qsub -q blades $script
    #qsub -q sunfire $script
    

done


