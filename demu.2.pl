#!/usr/bin/perl

use strict;
use warnings;
use feature qw( say );
#use Data::Dumper;
use Getopt::Long;

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

Unmatched and ambiguously matched reads are written into scrap_*.fastq.

Usage: demu.pl -oligos OLIGOS -f R1.fastq -r R2.fastq [OPTIONS]
Options:
	-orient	[0:1]	change R1 with R2 to bring the fragment
			orientation consistent if R1 matches REVERSE
			and R2 matches FORWARD
	-prefix	[str]	prefix for output filenames
	-matchr	[int]	search for match in range from the first to
			"length of primer plus matchr" bases of read [default=0]
	-maxm	[int]	maximum number of allowed mismatches between 
			read and template [default=1]
	-mask5	[int]	mask this number of positions from 5'-end 
			of the primer sequence before matching [default=0]
	-mask3	[int]	same as mask5 but for 3'-end of the primer
	-trim5	[0:1]	trim the masked 5'-end from the read sequence if matched; 
			trim5 is substracted from TRIMF and TRIMR
			to take into account the offset [default=1]
	-help|?		get this message

USAGE
#############################################################

my %opts = (
	matchr	=> 0,
	maxm	=> 1,
	mask3	=> 0,
	mask5	=> 0,
	trim5	=> 1
);

GetOptions (
	\%opts,
	'oligos=s',
	'f=s',
	'r=s',
	'prefix=s',
	'maxm=i',
	'mask5=i',
	'mask3=i',
	'trim5=i',
	'matchr=i',
	'orient',
	'help|?'
) || die($USAGE);

die($USAGE) if $opts{help};

open(OLI, '<'.$opts{oligos}) || die("Cannot open oligos file: $!");

my %g;
my %l;
my %fh;
my %c;
while(<OLI>) {
	next if /^#/ || /^$/;
	chomp;
	my @r=split;
	my $gid = defined($opts{prefix}) ? $opts{prefix}.'_'.$r[3] : $r[3];
	$l{$gid}[0] = length($r[1]) + $opts{matchr};
	$l{$gid}[1] = length($r[2]) + $opts{matchr};
	if( $opts{mask5} > 0) {
		($r[4], $r[5]) = map { $_ - $opts{mask5} } ($r[4], $r[5])
			if defined($r[4]);
		($r[1], $r[2]) = map { substr($_, $opts{mask5}) } ($r[1], $r[2])
	}
	($r[1], $r[2]) = map { substr($_, 0, -1*$opts{mask3}) } ($r[1], $r[2])
		if $opts{mask3} > 0;
	die( join ("\n",
		'Too many bases were masked.',
		'Match pattern is too ambiguous.',
		'Decrease mask3 and/or mask5 option values.',
		"\n")
		) if length($r[1]) <= 5 || length($r[2]) <= 5;
		
	
	$g{$gid} = defined($r[4]) ? [ @r[1,2,4,5] ] : [ @r[1,2] ];
	
	open(my $fo,'>'.$gid.'.R1.fastq') || die("Cannot open output file: $!");
	open(my $ro,'>'.$gid.'.R2.fastq') || die("Cannot open output file: $!");
	$fh{$gid} = [ $fo, $ro ];
}

open(my $fo,'>'.$opts{prefix}.'_scrap.R1.fastq') || die("Cannot open output file: $!");
open(my $ro,'>'.$opts{prefix}.'_scrap.R2.fastq') || die("Cannot open output file: $!");
$fh{'scrap'} = [$fo,$ro];

open(F,"<$opts{f}") || die("Cannot open input file: $!");
open(R,"<$opts{r}") || die("Cannot open input file: $!");

while ( (my @f = map $_ = <F>, 1 .. 4 )[0] ) {
	my @r = map $_ = <R>, 1 .. 4;
	my $gid = '';
	my $isamb = 0;
	my @mp;
	for(keys %g) {
		my @mpt = ();
		$mpt[0] = match(substr($f[1], 0, $l{$_}[0]),$g{$_}[0], $opts{maxm});
		$mpt[1] = match(substr($r[1], 0, $l{$_}[1]),$g{$_}[1], $opts{maxm});
		if ($mpt[0] && $mpt[1]) {
			if($gid ne '') { $isamb = 1; last }
			$gid = $_;
			@mp = @mpt;
		}
		$mpt[0] = match(substr($f[1], 0, $l{$_}[1]), $g{$_}[1], $opts{maxm});
		$mpt[1] = match(substr($r[1], 0, $l{$_}[0]), $g{$_}[0], $opts{maxm});
		if ($mpt[0] && $mpt[1]) {
			if($gid ne '') { $isamb = 1; last }
			if($opts{orient}) {
				my @t = @r;
				@r = @f;
				@f = @t;
				@mpt = reverse(@mpt);
				@{ $g{$_} }[2,3] = @{ $g{$_} }[3,2] if(defined( $g{$_}[2] ));
			}
			$gid = $_;
			@mp = @mpt;
		}
	}
	
	$gid = 'scrap' if( $gid eq '' || $isamb == 1 );
	
	#Trim sequences according to trim5, mask5, TRIMF/R options
	if( $gid ne 'scrap' && ($opts{trim5} == 1 || defined($g{$gid}[2])) ) {
		$mp[0] += $g{$gid}[2] if defined($g{$gid}[2]);
		$mp[1] += $g{$gid}[3] if defined($g{$gid}[3]);
		@f[1,3] = map { substr($_,$mp[0]) } @f[1,3];
		@r[1,3] = map { substr($_,$mp[1]) } @r[1,3];		
	}
	print { $fh{$gid}[0] } join("", @f);
	print { $fh{$gid}[1] } join("", @r);
	$c{$gid}++;
	$c{$gid."_amb"}++ if $isamb == 1;
}

close F;
close R;
for(keys %fh) {
	close $fh{$_}[0];
	close $fh{$_}[1];
}
for(keys %c){printf("%s\t%5d\n", $_, $c{$_});}
exit 0;

sub match {
#Finding the Longest Common String problem, snip from here:
#http://stackoverflow.com/questions/6258922/efficient-substring-matching-in-perl
#with adaptations to map the ambiguous bases

	my ($s, $t, $max_x) = @_;

	# Ambiguous nucleic acid residue codes are mapped to unambiguous ones
	my %IUB = (
		A => [qw(A)],	
		C => [qw(C)],
		G => [qw(G)],
		T => [qw(T)],
		U => [qw(U)],
		M => [qw(A C)],
		R => [qw(A G)],
		S => [qw(C G)],
		W => [qw(A T)],
		Y => [qw(C T)],
		K => [qw(G T)],
		V => [qw(A C G)],
		H => [qw(A C T)],
		D => [qw(A G T)],
		B => [qw(C G T)],
		N => [qw(A C G T)],
		X => [qw(A C G T)],
	);
	%IUB = map { $_ => join('', @{ $IUB{$_} }) } keys %IUB;

	my $m = my @s = unpack('(a)*', $s);
	my $n = my @t = unpack('(a)*', $t);

	push @s, ('?')x($n-1);

	my $best_x = $max_x + 1;
	my $best_i = 0;

	OUTER:
	for my $i (0..$m-1) {
		my $x = 0;
		for my $j (0..$n-1) {
			++$x if index($IUB{ $t[$j] }, $s[$i+$j]) < 0;
			next OUTER if $x >= $best_x;
		}
		$best_x = $x;
		$best_i = $i;
		last if !$best_x;
	}   

	if ($best_x > $max_x) {
		return undef;
	} else {
		return $best_i;
	}
}

