

#!/usr/bin/perl


use strict;
use IO::Uncompress::Gunzip;

my $VCF = "/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/fluidigm50k_final15jan2013_noIFCbarcode__has_phenotype.indiv_4377sites_41911indiv.vcf.gz";

my $hasPheno = "/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/has_phenotype.indiv";
my $noIchipButPheno = "/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/has_phenotype_noichip_control.indiv";

if ($#ARGV != 3) {die "Wrong number of arguments";}


open (ICAP, " < $hasPheno") or die "Cannot open $hasPheno";
my %icap = (); 
while (<ICAP>) {chomp $_;$icap{ $_ } = "yes";}
close (ICAP);


open (ICNP, " < $noIchipButPheno") or die "Cannot open $noIchipButPheno";
my %icnp = (); 
while (<ICNP>) {chomp $_;$icnp{ $_ } = "yes";}
close (ICNP);




my $chrQ = $ARGV[0];
my $startQ = $ARGV[1];
my $endQ = $ARGV[2];
my $code = $ARGV[3];


my $outputGeno = $code."_geno.tab";
my $outputDepth = $code."_depth.tab";
my $outputSummary = $code."_summary.tab";
my $sampleIDs = $code."_IDs.tab";

print "VCF file: $VCF\n";
print "Querying chr".$chrQ." between ".$startQ." and ".$endQ."\n";



############# start reading the VCF file
my $ptr = new IO::Uncompress::Gunzip($VCF) or die $!;
my $line;

open (OUT , " > $outputGeno") or die "Cannot open $outputGeno\n";
open (DEPTH , " > $outputDepth") or die "Cannot open $outputDepth\n";
open (SUM , " > $outputSummary") or die "Cannot open $outputSummary\n";
open (IDS , " > $sampleIDs") or die "Cannot open $sampleIDs\n";

print SUM "chromosome\tstart\tREF\tALT\tcode\tnalleles\tnhets\tnhoms\tMAF\tnhets.nic\tnhoms.nic\tMAF.nic\n";

while (defined ($line = $ptr->getline())) {
    if ($line  =~ /^.CHROM/) {last;}
}

chomp $line;
my @header = split('\t', $line);
my $nfields = $#header;
print "Number of fields $nfields\n"; 

    

################################################## parse the names of all the individuals
my @icaptab; my @icnptab;
my ($nsamplesicap, $nsamplesicnp)  = (0,0);
for (my $col = 9; $col <= $nfields; $col++) {  #for each column
    if (exists $icap {$header[ $col ] }) {  ##means it has a phenotype
	
	$icaptab[ $nsamplesicap ] = $col;
	$nsamplesicap++;
	print IDS $header[ $col ]."\n";

	if (exists $icnp {$header[ $col ] }) {  ##means a noIchip control
	    $icnptab[ $nsamplesicnp ] = $col;$nsamplesicnp++;
	} 


    }    
}

close (IDS);
print "List of individuals in $sampleIDs\n";
print "Number of individuals with phenotype, excluding the no ichip controls: ".scalar(@icaptab)."\n";
print "Number of no ichip control individuals: ".scalar(@icnptab)."\n";


################################
my $nvariant = 0;
LINE: while (defined ($line = $ptr->getline())) {
    
    my @spl = split('\t', $line);

    my $chrom = $spl[0];
    my $start = $spl[1];
    $chrom =~ s/chr//;
    
    if (($chrom eq $chrQ) && ($start > $startQ) && ($start < $endQ)) {
	$nvariant = $nvariant + 1;

	my $REF = $spl[3];
	my $ALT = $spl[4];	
	
	my @ALTALLELES = split(',', $ALT);
	my $nalleles = $#ALTALLELES + 1;
	my $code = $chrom."_".$start."_".$REF;

	print $chrom."\t".$start."\t",$REF."\t".$ALT." with ".$nalleles." alleles and code ".$code."\n";

	for (my $all = 0;  $all < $nalleles; $all++) {
	    print "Dealing with allele $all\n";
	    my $ALTLOC = $ALTALLELES[ $all ];
	    my $ALTLOCNUM = $all + 1;
	    print SUM $chrom."\t".$start."\t".$REF."\t".$ALTLOC."\t".$code."\t".$nalleles;
	 
	    ################### first parse for the samples that have ichip data
	    my $nhets = 0;
	    my $nhoms = 0;
	    my $nonmissing = 0;
	    
	    for (my $j = 0; $j <= $#icaptab; $j++) {  #for each column
		my $col = $icaptab[ $j ];
		my @geno = split(':', $spl[ $col]);
		my $locgeno = $geno[0];
			
		my @locgenosplit = split('/', $locgeno);
		my $altcount = ($locgenosplit[0] == $ALTLOCNUM) + ($locgenosplit[1] == $ALTLOCNUM);
		my $refcount = ($locgenosplit[0] == 0) + ($locgenosplit[1] == 0);

		if ($altcount + $refcount == 2) {
		    $nonmissing++;
		    print DEPTH $geno[2]."\n";

		    if ($j == 0) {print OUT $altcount;} else {print  OUT "\t".$altcount;} 
		    if ($altcount == 1)  {$nhets = $nhets + 1;}
		    if ($altcount == 2)  {$nhoms = $nhoms + 2;}
		} else {
		    if ($j == 0) {print OUT "NA";} else {print  OUT "\tNA";} 
		}
	    }
	    
	    my $MAF = ($nhets + 2*$nhoms)/(2*$nonmissing);
	    print SUM "\t".$nhets."\t".$nhoms."\t".$MAF;
	    
	    ################### then parse for the samples that do not have ichip data
	    $nhets = 0;
	    $nhoms = 0;
	    $nonmissing = 0;
	    for (my $j = 0; $j <= $#icnptab; $j++) {  #for each column
		my $col = $icnptab[ $j ];
		my @geno = split(':', $spl[ $col]);
		my $locgeno = $geno[0];
		my @locgenosplit = split('/', $locgeno);
		
		my $altcount = ($locgenosplit[0] == $ALTLOCNUM) + ($locgenosplit[1] == $ALTLOCNUM);
		my $refcount = ($locgenosplit[0] == 0) + ($locgenosplit[1] == 0);
		
		if ($altcount + $refcount == 2) {
		    $nonmissing++;
		    if ($altcount == 1)  {$nhets = $nhets + 1;}
		    if ($altcount == 2)  {$nhoms = $nhoms + 2;}
		}
	    }
	    $MAF = ($nhets + 2*$nhoms)/(2*$nonmissing);
	    print SUM "\t".$nhets."\t".$nhoms."\t".$MAF;

	    print OUT "\n";
	    print SUM "\n";
	    #if ($nvariant == 10) {last LINE;}
	} 
    }
}



$ptr->close();
close (OUT);
close (SUM);
close(DEPTH);

print "Output files are ".$outputGeno."  and ".$outputSummary."\n";



