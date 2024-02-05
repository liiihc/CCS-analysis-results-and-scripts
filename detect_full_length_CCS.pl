open (IN1,"<","blastn_mac.result");#blastn to MAC
open (IN2,"<","blastn_mic.result");#blastn to MIC
open (IN3,"<","blastn_mito.result");#blastn to Mito
open (OUT1,">","MAC_full_length_zmw.txt");#MAC results bed file
open (OUT2,">","MIC_full_length_zmw.txt");#MIC results bed file
open (OUT3,">","Mito_full_length_zmw.txt");#Mito results bed file
open (OUT4,">","MAC_clean_full_length_zmw.txt");#MAC clean results bed file delta<=50
while (defined ($_=<IN3>)){
        if($_=~/^m.+?\/(\d+?)\/ccs\t(.+?)\t(.+?)\t(.+?)\t.+?\t.+?\t.+?\t.+?\t.+?\t.+?\t(.+?)\t.+?\t.+?\t(.+?)\t.+?\n/){
		$zmw=$1;#qacc
		$chr=$2;#sacc
		$nident=$4;#nident
		$score=$6;#score
		$length=$3;#length
		$qlen=$5;#qlen
		$ident=$nident/$length;
		$coverage=$length/$qlen;
		if ($coverage>=0.95 && $ident>=0.98){
			$mito{$zmw}=$score;
			print OUT3 "$zmw\n"}}}
while (defined ($_=<IN2>)){
        if($_=~/^m.+?\/(\d+?)\/ccs\t(.+?)\t(.+?)\t(.+?)\t.+?\t.+?\t.+?\t.+?\t.+?\t.+?\t(.+?)\t.+?\t.+?\t(.+?)\t.+?\n/){
                $zmw=$1;#qacc
                $chr=$2;#sacc
                $nident=$4;#nident
                $score=$6;#score
                $length=$3;#length
                $qlen=$5;#qlen
		$ident=$nident/$qlen;
                $coverage=$qlen/$length;
                if ($coverage>=0.95 && $ident>=0.98){
                        $mic{$zmw}=$score;print OUT2 "$zmw\n"}}}
while (defined ($_=<IN1>)){
        if($_=~/^m.+?\/(\d+?)\/ccs\t(.+?)\t(.+?)\t(.+?)\t.+?\t.+?\t.+?\t.+?\t.+?\t.+?\t(.+?)\t.+?\t.+?\t(.+?)\t.+?\n/){
                $zmw=$1;#qacc
                $chr=$2;#sacc
                $nident=$4;#nident
                $score=$6;#score
                $length=$3;#length
                $qlen=$5;#qlen
		$ident=$nident/$qlen;
                $coverage=$qlen/$length;
                if ($coverage>=0.95 && $ident>=0.98){print OUT1 "$zmw\n";
                        if (exists $mic{$zmw}){
			$delta=$mic{$zmw}-$score;
			if ($delta<=50){print OUT4 "$zmw\n"}}
			else {print OUT4 "$zmw\n"}}}}
