#!/usr/bin/perl
use warnings;

if (@ARGV < 2){
	print "\n\n  usage: split.freebayes.regions.file.pl  regions.txt  data-per-file \n\n";
	exit 1;
}
my $infile = $ARGV[0];
my $dataperfile=$ARGV[1];

open (IN, $infile);

my $index=1;
my $total=0; 
my $out;
open ($out, ">regions_". "$dataperfile". "_" ."$index");

while (my $line = <IN>) {
	# scaffold19267_size1087:0-1087
	my ($name, $range) = split(":", $line);
	my ($start, $stop) = split("-", $range);
	my $size = $stop - $start; 
	$total = $total+$size; 
	if ($total > $dataperfile) {
		close ($out);
		$index++;
		open ($out, ">regions_". "$dataperfile". "_" ."$index");
		$total= $size;
		
	}
	print $out $name ."\t". $start . "\t" . $stop; # s/[:-]/\t/g'
	#print $out $line;
}

close (IN);


