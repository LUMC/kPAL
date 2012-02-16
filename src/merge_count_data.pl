#!/usr/local/bin/perl -w

open (REF1, $ARGV[0]) || die "cannot open file\n";
open (REF2, $ARGV[1]) || die "cannot open file\n";
open (REF3, $ARGV[2]) || die "cannot open file\n";
open (REF4, $ARGV[3]) || die "cannot open file\n";
open (REF5, $ARGV[4]) || die "cannot open file\n";
open (REF6, $ARGV[5]) || die "cannot open file\n";
open (REF7, $ARGV[6]) || die "cannot open file\n";
open (REF8, $ARGV[7]) || die "cannot open file\n";

open (OUT, ">owls_counts_merged.txt");

print OUT "$ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7]\n";

while (<REF1>){
	$_ =~ s/\n|\r//g;
	$line2 = <REF2>;
	$line2 =~ s/\n|\r//g;
	$line3 = <REF3>;
	$line3 =~ s/\n|\r//g;
	$line4 = <REF4>;
	$line4 =~ s/\n|\r//g;	
	$line5 = <REF5>;
	$line5 =~ s/\n|\r//g;
	$line6 = <REF6>;
	$line6 =~ s/\n|\r//g;
	$line7 = <REF7>;
	$line7 =~ s/\n|\r//g;
	$line8 = <REF8>;
	$line8 =~ s/\n|\r//g;
	print OUT "$_,$line2,$line3,$line4,$line5,$line6,$line7,$line8\n";
}
