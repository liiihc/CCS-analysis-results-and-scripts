open (IN1,"<","CCS_hifi_reads.hifi.pbmm2.sam");#sam
open (IN2,"<","All_A_IPD.txt");#Exact IPDr from all csv files to one file.
open (IN6,"<","filtered_zmwid.txt");
open (OUT1,">","ccs2genome.txt");
while (defined ($_=<IN1>)){
		if ($_=~/\/(\d+)\/ccs\s+(\d+)\s+(.+?)\s+(\d+)\s+.+?\s+(.+?)\s+.+?\s+.+?\s+.+?\s+(\w+)\s/ && ($2==0 || $2==16)){	
		$strand{$1}=$2;
		$chr{$1}=$3;
		$pos{$1}=$4;
		$F{$1}=$5;
		$length{$1}=length($6);
		}
		elsif ($_=~/\/(\d+)\/ccs\s+(\d+)\s+(.+?)\s+(\d+)\s+.+?\s+(.+?)\s+.+?\s+.+?\s+.+?\s+(\w+)\s/ && ($2==2048 || $2==2064) && $_=~/SA:Z/){
		$strand_secondmap{$1}=$2;
                $chr_secondmap{$1}=$3;
                $pos_secondmap{$1}=$4;
                $F_secondmap{$1}=$5;
                $length_secondmap{$1}=length($6)}
		}
while (defined ($_=<IN6>)){if ($_=~/(\d+)/){
                $mac{$1}=1;}}
$la=0;$yes=0;
while (defined ($_=<IN2>)){
	if ($_=~/^(\d+)\s+/ && (($strand{$1}==0) && (exists $mac{$1})))
	{
	if ($_=~/^(\d+)\s+(\d+)\s+(\d+)\s+(A.)\s+(.+?)\s+(.+)/ && (exists $F{$1}) ){
		$pos=$2;$IPD=$5;$zmw=$1;$num=0;$strand=$3;$F=$F{$zmw};$I=0;$D=0;$motif=$4;$motif2=$6;
	while ($F=~/^(\d+[\=|I|X|S|D])(.*)/){
		$tmp=$1;$F=$2;
		if ($tmp=~/(\d+)[\=|X]/){
			$start=$num;
			$num+=$1;
			if ($pos<=$num && $pos>$start && ($tmp=~/\=/) && $strand==0){$pos1=$pos+$pos{$zmw}-$I+$D-1;print OUT1 "$zmw\t$chr{$zmw}\t$pos1\tW\t$IPD\t$motif\t$motif2\n";}
			elsif ($pos<=$num && $pos>$start && ($tmp=~/\=/) && $strand==1){$pos1=$pos+$pos{$zmw}-$I+$D-1;print OUT1 "$zmw\t$chr{$zmw}\t$pos1\tC\t$IPD\t$motif\t$motif2\n";}}
		elsif ($tmp=~/(\d+)[I|S]/){
			 $num+=$1;$I+=$1;}
		elsif ($tmp=~/(\d+)D/){
			$D+=$1;}
		}}}
	if ($_=~/^(\d+)\s+/ && (($strand{$1}==16) && (exists $mac{$1})))
        {        if ($_=~/^(\d+)\s+(\d+)\s+(\d+)\s+(A.)\s+(.+?)\s+(.+)/ && (exists $F{$1})){
                $pos=$length{$1}-$2+1;$IPD=$5;$zmw=$1;$num=0;$strand=$3;$F=$F{$zmw};$I=0;$D=0;$motif=$4;$motif2=$6;
       while ($F=~/^(\d+[\=|I|X|S|D])(.*)/){
		$tmp=$1;$F=$2;
                if ($tmp=~/(\d+)[\=|X]/){
                        $start=$num;
                        $num+=$1;;
                        if ($pos<=$num && $pos>$start && ($tmp=~/\=/) && $strand==1){$pos1=$pos+$pos{$zmw}-$I+$D-1;print OUT1 "$zmw\t$chr{$zmw}\t$pos1\tW\t$IPD\t$motif\t$motif2\n";}
                        elsif ($pos<=$num && $pos>$start && ($tmp=~/\=/) && $strand==0){$pos1=$pos+$pos{$zmw}-$I+$D-1;print OUT1 "$zmw\t$chr{$zmw}\t$pos1\tC\t$IPD\t$motif\t$motif2\n";}}
                elsif ($tmp=~/(\d+)[I|S]/){
                         $num+=$1;$I+=$1;}
                elsif ($tmp=~/(\d+)D/){
                        $D+=$1;}
		}}}
	if ($_=~/^(\d+)\s+/ && (($strand_secondmap{$1}==2048) && (exists $mac{$1})))
        {
        if ($_=~/^(\d+)\s+(\d+)\s+(\d+)\s+(A.)\s+(.+?)\s+(.+)/ && (exists $F{$1})){
                $pos=$2;$IPD=$5;$zmw=$1;$num=0;$strand=$3;$F=$F_secondmap{$zmw};$I=0;$D=0;$motif=$4;$motif2=$6;
        while ($F=~/^(\d+[\=|I|X|S|D])(.*)/){
                $tmp=$1;$F=$2;
                if ($tmp=~/(\d+)[\=|X]/){
                        $start=$num;
                        $num+=$1;
                        if ($pos<=$num && $pos>$start && ($tmp=~/\=/) && $strand==0){$pos1=$pos+$pos_secondmap{$zmw}-$I+$D-1;print OUT1 "$zmw\t$chr_secondmap{$zmw}\t$pos1\tW\t$IPD\t$motif\t$motif2\n";}
                        elsif ($pos<=$num && $pos>$start && ($tmp=~/\=/) && $strand==1){$pos1=$pos+$pos_secondmap{$zmw}-$I+$D-1;print OUT1 "$zmw\t$chr_secondmap{$zmw}\t$pos1\tC\t$IPD\t$motif\t$motif2\n";}}
                elsif ($tmp=~/(\d+)[I|S]/){
                         $num+=$1;$I+=$1;}
                elsif ($tmp=~/(\d+)D/){
                        $D+=$1;}
                }}}
        if ($_=~/^(\d+)\s+/ && (($strand_secondmap{$1}==2064) && (exists $mac{$1})))
        {        if ($_=~/^(\d+)\s+(\d+)\s+(\d+)\s+(A.)\s+(.+?)\s+(.+)/ && (exists $F_secondmap{$1})){
                $pos=$length_secondmap{$1}-$2+1;$IPD=$5;$zmw=$1;$num=0;$strand=$3;$F=$F_secondmap{$zmw};$I=0;$D=0;$motif=$4;$motif2=$6;
       while ($F=~/^(\d+[\=|I|X|S|D])(.*)/){
                $tmp=$1;$F=$2;
                if ($tmp=~/(\d+)[\=|X]/){
                        $start=$num;
                        $num+=$1;;
                        if ($pos<=$num && $pos>$start && ($tmp=~/\=/) && $strand==1){$pos1=$pos+$pos_secondmap{$zmw}-$I+$D-1;print OUT1 "$zmw\t$chr_secondmap{$zmw}\t$pos1\tW\t$IPD\t$motif\t$motif2\n";}
                        elsif ($pos<=$num && $pos>$start && ($tmp=~/\=/) && $strand==0){$pos1=$pos+$pos_secondmap{$zmw}-$I+$D-1;print OUT1 "$zmw\t$chr_secondmap{$zmw}\t$pos1\tC\t$IPD\t$motif\t$motif2\n";}}
                elsif ($tmp=~/(\d+)[I|S]/){
                         $num+=$1;$I+=$1;}
                elsif ($tmp=~/(\d+)D/){
                        $D+=$1;}
                }}}
}

