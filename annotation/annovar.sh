VCF=/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/fluidigm50k_final15jan2013_noIFCbarcode__has_phenotype.indiv_4377sites_41911indiv.vcf.gz


software=/data_n2/vplagnol/Software
annovar=${software}/annovar_Feb2013

annoDB=${annovar}/humandb_hg19
convert2annovar=${annovar}/convert2annovar.pl
sumAnno="${annovar}/summarize_annovar.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype gencodegene "
#sumAnno="${annovar}/summarize_annovar.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 -exonicsplicing"

VCFclean=temp/annovar/VCFclean.vcf

script=cluster/submission/annovar.sh

code=`basename $VCF .vcf.gz`
echo $code

echo "
#zcat $VCF | cut -f 1-10 >  temp/annovar/${code}.vcf
#perl scripts/annotation/split_vcf_for_annovar.pl temp/annovar/${code}.vcf $VCFclean

${convert2annovar} --format vcf4 $VCFclean > temp/annovar/${code}

grep 123377437 temp/annovar/${code}

${sumAnno} -buildver hg19 temp/annovar/${code} ${annoDB}

" > $script


echo $script
