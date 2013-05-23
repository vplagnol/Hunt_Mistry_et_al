
use strict;
#my $VCF = "temp/annovar/fluidigm50k_final15jan2013_noIFCbarcode__has_phenotype.indiv_4377sites_41911indiv.vcf";
#my $output = "temp/annovar/VCFclean.vcf";

my $VCF = $ARGV[0];
my $output = $ARGV[1];

open (INP, " < $VCF");
open (OUT, " > $output");
while (<INP>) {

    my @spl = split('\t', $_);

    if ($spl[ 4 ] =~ /,/) {
	my $alleles = $spl[4];
	
	my @altalleles = split(',', $alleles);
	my $nalleles = $#altalleles;
	#print "Number of alleles: $nalleles ".$altalleles[0]."\t".$altalleles[1]."\n";

	for (my $i = 0; $i <= $nalleles; $i++) {
	    print OUT $spl[ 0 ];
	    for (my $j = 1; $j <= 3; $j++) {print OUT "\t".$spl[ $j ];}
	    print OUT "\t".$altalleles[ $i ];
	    for (my $j = 5; $j != 10; $j++) {print OUT "\t".$spl[ $j ];}
	} 
    } else { print OUT $_;}

}

close (INP);
close (OUT);
