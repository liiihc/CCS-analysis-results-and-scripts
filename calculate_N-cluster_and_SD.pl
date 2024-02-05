use strict;
my (%hash,$pos,@files,%strand,%csv,$sequence,$strand,%strand,$zmw,%N,$IPD,$num,@csv,$zwm,$pos2,$tmp,$SD_sum_ws,$SD_sum_cr,$ave_ws,$ave_cr,$num_ws,$num_cr,$SD_ws,$SD_cr);
@files=glob "ipdSummary_csv/*.csv";
open (IN3,"<","CCS_hifi_reads.pbmm2.sam");
open (OUT1,">","Ncluster_test_TEST.txt");
open (OUT2,">","SD_test_TEST.txt");
while (defined ($zwm=<IN3>)){
	if ($zwm=~/\/(\d+)\/ccs\s+0\s+/ && $zwm!~/SA:Z/){
			$strand{$1}=0}
	if ($zwm=~/\/(\d+)\/ccs\s+16\s+/ && $zwm!~/SA:Z/){
                        $strand{$1}=1}}
my ($IPD_ws,$IPD_cr,$num_ws,$num_cr,@IPD_cr,@IPD_ws);
foreach $_(@files){
        if ($_=~/\/(\d+).csv/){$zwm=$1}
        open (IN1,"<","$_");
        open (IN2,"<","$_");
	@IPD_ws=();$IPD_ws=0;$num_ws=0;$ave_ws=0;$ave_cr=0;$SD_sum_ws=0;@IPD_cr=();$IPD_cr=0;$num_cr=0;$SD_sum_cr=0;
        while (defined ($sequence=<IN1>)){
		if ($sequence=~/"m.+?\/(\d+)\/ccs",(\d+),(.),A,.+?,.+?,.+?,.+?,(.+?),/ && (exists $strand{$1}) && $4<2.8){
			if ($3==$strand{$1}){
				push (@IPD_ws,$4);
				$IPD_ws+=$4;
				$num_ws+=1;}
			else{push (@IPD_cr,$4);
                                $IPD_cr+=$4;
                                $num_cr+=1;}}
                if ($sequence=~/"m.+?\/(\d+)\/ccs",(\d+),(.),[T|C|G],.+?,.+?,.+?,.+?,(.+?),/ && $4>=2.8 ){
                        $pos=join "~","$1","$2","$3";
			$N{$pos}=1;
                      }}
		if ($num_ws>0 && $num_cr>0){
			$ave_ws=$IPD_ws/$num_ws;
			$ave_cr=$IPD_cr/$num_cr;
			foreach $_(@IPD_ws){
                        $SD_sum_ws+=($_-$ave_ws)**2}
			foreach $_(@IPD_cr){
                        $SD_sum_cr+=($_-$ave_cr)**2}
                $SD_ws=sqrt($SD_sum_ws/$num_ws);
		$SD_cr=sqrt($SD_sum_cr/$num_cr);
		print OUT2 "$zwm\t$SD_ws\t$SD_cr\n";}
	while (defined ($sequence=<IN2>)){
		if ($sequence=~/"m.+?\/(\d+)\/ccs",(\d+),(.),[T|C|G],.+?,.+?,.+?,.+?,(.+?),/ && $4>=2.8 ){
			$zwm=$1;
			$pos=$2;
			$strand=$3;
			$num=0;$pos2=0;
			foreach $_(0..25){
				if ($strand==$strand{$zmw}){
				$pos2=$pos+$_;
				$tmp=join "~","$zwm","$pos2","$strand";
				if (exists $N{$tmp}){
					$num++;}
			if ($num>=3){
				$hash{$zwm}=1;
			}}
				else {
				$pos2=$pos-$_;
                                $tmp=join "~","$zwm","$pos2","$strand";
                                if (exists $N{$tmp}){
                                        $num++;}
                        if ($num>=3){
                                $hash{$zwm}=1;
				}
				}}
				}}
	close IN1;
	close IN2;
	undef %N}
foreach $_(keys %hash){if ($hash{$_}==1){print OUT1 "$_\n"}}
