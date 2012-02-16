#!/usr/local/bin/perl -w

open (IN, "owls_counts_merged.txt");

$line = <IN>;
@frags = split(",",$line);

foreach $frag (@frags){
	open ($frag, ">$frag");
}

$line = <IN>;
@temp = split(",",$line);
$kmer = $temp[0];
$line = <IN>;
$line = <IN>;

while (<IN>){
	$_ =~ s/\r|\n//g;
	@frags2 = split(",",$_);
	$flag = 0;
	foreach $frag2 (@frags2){
		if ($frag2 == 0){
			$flag = 1;
		}
	}
	if ($flag == 1){
		for ($i = 0; $i < @frags; $i++){
			$temp = $frags[$i];
			print $temp "0\n";
		}
	} else {
		for ($i = 0; $i < @frags; $i++){
			$temp = $frags[$i];
			print $temp $frags2[$i]."\n";
		}
	}
}

for ($i = 0; $i < @frags; $i++){
	close $frags[$i];
	push (@files, $frags[$i]);
	print "close $frags[$i]\n";
}
print "printed count files\n";
print "@files";
foreach $file (@files){
	print "$file\n";
	$total = 0;
	$nonZero = 0;
	open FILE, $file or die $!;
	@data = <FILE>;
	foreach (@data){
		$_ =~ s/\n|\r//g;
		$total += $_;
		if ($_ != 0){
			$nonZero += 1;
		}
	}
	
	close FILE;
	open (FILE, ">$file") or die $!;
	print FILE "$kmer\n";
	print FILE "$total\n";
	print FILE "$nonZero\n";
	for($i=0;$i<@data;$i++){
		print FILE "$data[$i]\n";
	}
	@data=();
	close FILE;
}
