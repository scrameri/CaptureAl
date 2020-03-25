#!/usr/bin/perl

# Author: Stefan Zoller, GDC Zurich
# Submit multiple command-line commands in parallel in a controlled fashion.
# The commands have to be written one per line in a commands file.
# Do not run more than 20 commands in parallel and make sure RAM and CPUs are free.
# Usage: submit-commands-var.pl number-of-parallel-jobs file-with-one-command-per-line.

use strict;
use warnings;
use Paranoid::Process;

if (@ARGV < 2) {
         print "\nusage: submit-commands-var.pl  number-of-parallel-jobs  file-with-commands  \n" ;
     print "  where file-with-commands contains one command per line\n\n";
     exit 1;
}

my $max=$ARGV[0];
my $commandfile=$ARGV[1];

$SIG{CHLD}  = \&Paranoid::Process::sigchld;
Paranoid::Process::MAXCHILDREN = $max;

open (COMMANDFILE, "$commandfile") or die "cannot find commandfile $commandfile";
my @commands = <COMMANDFILE>;
my $pid=0;
my $command_n=0;
foreach my $command (@commands) {
	unless ($pid = Paranoid::Process::pfork() ) {
		 print "starting command: $command\n";
                system("$command");
                exit 0;
        }
}
sleep 2;
print "all commands submitted\n";
