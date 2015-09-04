#!/usr/bin/perl

use strict;
use warnings;

use Bio::PrimarySeq;
use Bio::Tools::IUPAC;
use Getopt::Long;
use Data::Dumper;

#############################################################
my $USAGE = <<USAGE;
Script to demultiplex the Illumina fastq files according to
exact match of primers with the 5'-ends of forward and reverse
reads. Script takes the R1 and R2 fastq files and the oligos file
as an input.

Oligos is a tab-delimited file of the following format:

#Col1	Col2	Col3	Col4	Col5	Col6
PRIMER	FORWARD	REVERSE	GROUP	TRIMF	TRIMR

Where FORWARD and REVERSE are the IUPAC sequences of respective primers,
and GROUP - is a group identificator for demultiplexing.
Script scans the forward and reverse reads and writes them into the output files
with GROUP prefix if R1 and R2 match the respective primer sequences.
TRIMF and TRIMR - number of 5'-end positions to trim from the R1 and R2,
respectively (optional).

Unmatched reads are written into scrap_*.fastq.

Usage: demu.pl -oligos OLIGOS -f R1.fastq -r R2.fastq [OPTIONS]
Options:
	-orient [0:1]	change R1 with R2 to bring the fragment
			orientation consistent if R1 matches REVERSE
			and R2 matches FORWARD
	-prefix [str]	prefix for output filenames

USAGE
#############################################################
my %opts;

GetOptions( \%opts, 'oligos=s', 'prefix=s', 'f=s', 'r=s', 'orient', 'help');

if($opts{help})
{
	print $USAGE;
	exit 0;
}

open(OLI, '<'.$opts{oligos}) || die("Cannot open oligos file: $!");

my %g;
my %fh;
my %c;
while(<OLI>)
{
	next if /^#/ || /^$/;
	chomp;
	my @r=split;
	my $f = Bio::PrimarySeq->new(-seq => $r[1], -alphabet => 'dna');
	my $r = Bio::PrimarySeq->new(-seq => $r[2], -alphabet => 'dna');
	my $iupac = Bio::Tools::IUPAC->new(-seq => $f);
	my $freg = $iupac->regexp();
	$iupac = Bio::Tools::IUPAC->new(-seq => $r);
	my $rreg = $iupac->regexp();
	my $gid = defined($opts{prefix}) ? $opts{prefix}.'_'.$r[3] : $r[3];
	$g{$gid} = defined($r[4]) ? [ $freg, $rreg, $r[4], $r[5] ] : [ $freg, $rreg ];
	
	open(my $fo,'>'.$gid.'.R1.fastq') || die("Cannot open output file: $!");
	open(my $ro,'>'.$gid.'.R2.fastq') || die("Cannot open output file: $!");
	$fh{$gid} = [ $fo, $ro ];
}

open(my $fo,'>'.$opts{prefix}.'_scrap.R1.fastq') || die("Cannot open output file: $!");
open(my $ro,'>'.$opts{prefix}.'_scrap.R2.fastq') || die("Cannot open output file: $!");
$fh{'scrap'} = [$fo,$ro];

open(F,"<$opts{f}") || die("Cannot open input file: $!");
open(R,"<$opts{r}") || die("Cannot open input file: $!");

while((my ($fid, $fs, $fpl, $fq) = map $_ = <F>, 1 .. 4 )[0])
{
	my ($rid, $rs, $rpl, $rq) = map $_ = <R>, 1 .. 4;
	my $gid = '';
	for(keys %g)
	{
		my ($regf, $regr) = @{$g{$_}}[0,1];
		
		if( $fs =~ /^$regf/ && $rs =~ /^$regr/ )
		{
			$gid=$_;
			last;
		}elsif($fs =~ /^$regr/ && $rs =~ /^$regf/ )
		{
			if($opts{orient})
			{
				my @t = ($fid, $fs, $fq);
				($fid, $fs, $fq) =  ($rid, $rs, $rq);
				($rid, $rs, $rq) =  @t;
				if(defined( $g{$_}[2] ))
				{
					my $t = $g{$_}[2];
					$g{$_}[3] = $g{$_}[2];
					$g{$_}[2] = $t;
				}
			}
			$gid=$_;
			last;
		}
	}
	
	$gid = 'scrap' if $gid eq '';
	($fs, $fq,  $rs, $rq) = &trim(\$g{$gid}, $fs, $fq, $rs, $rq) 
		if( $gid ne 'scrap' && (defined($g{$gid}[2]) || defined($g{$gid}[3])) );
	print { $fh{$gid}[0] } join("", $fid, $fs, $fpl, $fq);
	print { $fh{$gid}[1] } join("", $rid, $rs, $rpl, $rq);
	$c{$gid}++;
}

close F;
close R;
for(keys %fh)
{
	close $fh{$_}[0];
	close $fh{$_}[1];
}
for(keys %c){printf("%s\t%5d\n", $_, $c{$_});}
exit 0;

sub trim()
{
	my ($p,$fs,$fq,$rs,$rq) = @_;
	if(defined($$p->[2]))
	{
		$fs=substr($fs,$$p->[2]);
		$fq=substr($fq,$$p->[2]);
	}
	if(defined($$p->[3]))
	{
		$rs=substr($rs,$$p->[3]);
		$rq=substr($rq,$$p->[3]);
	}
	return ($fs, $fq, $rs, $rq);
}