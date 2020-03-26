#! /usr/bin/perl -w

use String::Unquotemeta;
use File::Basename;
use IO::Handle;
use Getopt::Std;

# Extract reads from fastq file.  

# Author: Stefan Zoller, GDC


my %options=();
getopts("f:r:",\%options);


 my $usage = <<'END_USAGE';

      usage:   extract-reads-from-fastq.pl -f <fastq file>  -r <file with full path/filenames of read-ID files>  
      example: extract-reads-from-fastq.pl -f sample1.read1.fastq  -r readID.files
 
END_USAGE

if ((not defined $options{"f"}) || (not defined $options{"r"}) ) {
                print $usage; 
                exit 1;
}


my $idfiles = $options{"r"};
my $in = $options{"f"};
my($fqfile, $directories, $suffix) = fileparse( $in );
my $outfile = "";

my %headers;
my $numreads=0;
STDOUT->printflush("\nreading fastq reads\n");

my %seqs = readFastq($fqfile);
STDOUT->printflush("number of reads found: $numreads\n");


open (IDFILES, "$idfiles");
while (my $idfile = <IDFILES>){
	chomp $idfile;
	print "file with IDs: $idfile\n";
	open (IDs, "$idfile");
	$outfile =  "extracted_reads_$fqfile.$idfile.fastq";
	open (OUT, ">$outfile");
	my $i=0;
	while ($readID = <IDs>){
		chomp $readID; 
		$readID=~s/\s+$//;
		print "looking for read: $readID \n";
		if( defined $seqs{"$readID"}){
			print "  found \n";
			print OUT $seqs{"$readID"};
			# delete $seqs{$random};
			# STDOUT->printflush("\rread $i saved"); 
		}
	}
	close IDs;
	close OUT; 
}
close IDFILES;

sub readFastq {
        my %sequences;    # make this a local variable
        my $filename = $_[0];
	my $id="";
	$numreads=0;
        open (FILE, $filename) or die "Cannot open $filename";
	
        while ($line = <FILE>) {
		$numreads++;
                chomp $line;
		($id) = $line =~ m/(.+?) .*/ ;
		#print "putting into hash read: $id \n"; 
		$seq=<FILE>;
                chomp $seq;
		$inbe=<FILE>;
                chomp $inbe;
		$qual=<FILE>;
                chomp $qual;
                $sequences{"$id"} = $line ."\n".$seq  ."\n". $inbe  ."\n". $qual   ."\n";
        }
        return %sequences;
}

