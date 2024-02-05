use strict;
my($m6A,$tmp,$ipd,%hash,%hash2,@keys,$chr,$pos,$ori,$penetration,@array);
open (IN1,"<","6mAorA_inApT_doubleCutoff");
open (OUT1,">","penetration.txt");
while (defined ($m6A=<IN1>)){
	if ($m6A=~/(chr_\d+)\s+(\d+)\s+(.)\s+(.+)\s+A./){
		$tmp=join "~","$1","$2","$3";
		$ipd=$4;
		$hash{$tmp}++;
		if ($ipd=~/^6mA/){
		$hash2{$tmp}++;
		}}}
@keys=sort {$a <=> $b} keys(%hash);
foreach $_(@keys){
	@array=split "~",$_;
	$chr=$array[0];
	$pos=$array[1];
	$ori=$array[2];
	if (exists $hash2{$_}){
	$penetration=$hash2{$_}/$hash{$_};
	print OUT1 "$chr\t$pos\t$ori\t$hash2{$_}\t$hash{$_}\t$penetration\n"}
	else {print OUT1 "$chr\t$pos\t$ori\t0\t$hash{$_}\t0\n"}}
	
