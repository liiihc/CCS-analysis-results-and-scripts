#sed s/CUTOFFA/2.38/g fullhemi.pl | s/CUTOFFC/1.57/g | s/CUTOFFW/1.57/g > fullhemi.final.pl
#CUTOFFA is the A IPDr cut-off without shift;CUTOFFC and CUTOFFW are shifted.
open (IN1,"<","ccs2genome.txt");
open (IN2,"<","ccs2genome.txt");
open (IN3,"<","ccs2genome.txt");
open (IN4,"<","ccs2genome.txt");
open (OUT1,">","6mAratio");
open (OUT2,">","6mAorA_inApT_doubleCutoff");
open (OUT3,">","6mAorA_inApAGC_singleCutoff");
open (OUT4,">","6mAorA_inApO_singleCutoff");
open (OUT5,">","6mAorA_inApT_singleStrand_singleCutoff");
$num_AT=0;
while (defined ($A=<IN1>)){
	$num_A++;}
while (defined ($A=<IN2>)){
	if ($A=~/(\d+)\s+(chr_\d+)\s+(\d+)\s+(.)\s+(.+)\s+AT\s+/){
		$pos=join "~","$1","$2","$3","$4";
		$hash{$pos}=$5;
		}}
while (defined ($A=<IN3>)){
	if ($A=~/(\d+)\s+(chr_\d+)\s+(\d+)\s+(.)\s+(.+)\s+AT\s+(.+)/){
		if ($4 eq "W"){$pos=$3+1;
                	$pos1=join "~","$1","$2","$pos","C";
                	if (! exists $hash{$pos1} && $5>=CUTOFFA){print OUT5 "$1\t$2\t$3\t$4\t6mA\tAT\t$6\n";$num_6mApN++}
			elsif (! exists $hash{$pos1} && $5<CUTOFFA){print OUT5 "$1\t$2\t$3\t$4\tno6mA\tAT\t$6\n"}}
		elsif ($4 eq "C"){$pos=$3-1;
			$pos1=join "~","$1","$2","$pos","W";
                	if (! exists $hash{$pos1} && $5>=CUTOFFA){print OUT5 "$1\t$2\t$3\t$4\t6mA\tAT\t$6\n";$num_6mApN++}
                	elsif (! exists $hash{$pos1} && $5<CUTOFFA){print OUT5 "$1\t$2\t$3\t$4\tno6mA\tAT\t$6\n"}}}
	elsif ($A=~/(\d+)\s+(chr_\d+)\s+(\d+)\s+(.)\s+(.+)\s+AN\s+(.+)/){
		if ($5>=CUTOFFA){print OUT4 "$1\t$2\t$3\t$4\t6mA\tAN\t$6\n";$num_6mApN++}
		else {print OUT4 "$1\t$2\t$3\t$4\tno6mA\tAN\t$6\n"}}
	elsif ($A=~/(\d+)\s+(chr_\d+)\s+(\d+)\s+(.)\s+(.+)\s+(A[G|C|A])\s+(.+)/){
		if ($5>=CUTOFFA){print OUT3 "$1\t$2\t$3\t$4\t6mA\t$6\t$7\n";$num_6mApN++}
                else {print OUT3 "$1\t$2\t$3\t$4\tno6mA\t$6\t$7\n"}}}
while (defined ($A=<IN4>)){
	if ($A=~/(\d+)\s+(chr_\d+)\s+(\d+)\s+W\s+(.+)\s+AT\s+(.+)/){
		$pos=$3+1;
		$pos1=join "~","$1","$2","$pos","C";
		if (exists $hash{$pos1}){$num_AT+=2;
		if ($4<CUTOFFA && $hash{$pos1}<CUTOFFA){$non+=2;print OUT2 "$1\t$2\t$3\tW\tno6mA\tAT\t$5\n"}
		elsif($4>=CUTOFFA && $hash{$pos1}<CUTOFFC){$hemi_W++;$non+=1;print OUT2 "$1\t$2\t$3\tW\t6mA\tAT\t$5\n"}
		elsif($4<CUTOFFW && $hash{$pos1}>=CUTOFFA){$hemi_C++;$non+=1;print OUT2 "$1\t$2\t$3\tW\tno6mA\tAT\t$5\n"}
		else {print OUT2 "$1\t$2\t$3\tW\t6mA\tAT\t$5\n"}}}
	elsif ($A=~/(\d+)\s+(chr_\d+)\s+(\d+)\s+C\s+(.+)\s+AT\s+(.+)/){
                $pos=$3-1;
                $pos1=join "~","$1","$2","$pos","W";
                if (exists $hash{$pos1}){
                if ($4<CUTOFFA && $hash{$pos1}<CUTOFFA){print OUT2 "$1\t$2\t$3\tC\tno6mA\tAT\t$5\n"}
                elsif($4>=CUTOFFA && $hash{$pos1}<CUTOFFC){print OUT2 "$1\t$2\t$3\tC\t6mA\tAT\t$5\n"}
                elsif($4<CUTOFFW && $hash{$pos1}>=CUTOFFA){print OUT2 "$1\t$2\t$3\tC\tno6mA\tAT\t$5\n"}
                else {print OUT2 "$1\t$2\t$3\tC\t6mA\tAT\t$5\n"}}}
	}
$full=$num_AT-$non-$hemi_C-$hemi_W;
$m6ApT=$hemi_C+$hemi_W+$full;
$m6ApT_ApT=$m6ApT/$num_AT;
$m6A_A=($m6ApT+$num_6mApN)/$num_A;
$m6ApTratio=$m6ApT/($m6ApT+$num_6mApN);
$m6ApT_A=$m6ApT/$num_A;
$hemi_Cratio=$hemi_C/$m6ApT;
$hemi_Wratio=$hemi_W/$m6ApT;
$full_ratio=($full/$m6ApT)/2;
print OUT1 "full:$full\nhemiW:$hemi_W\nhemiC:$hemi_C\nnon:$non\nA:$num_A\nAT:$num_AT\n6mApT:$m6ApT\n6mApN:$num_6mApN\n";


