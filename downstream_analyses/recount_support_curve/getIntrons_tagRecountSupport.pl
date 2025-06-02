#!usr/bin/perl
#given a file with transcripts and exons, output transcripts and introns together with recount support
#usage--> perl getIntrons_tagRecountSupport.pl ../IntronChainMT/Hv3_splicedmasterTable.gtf ../IntronChainMT/TrmapHv3_Pairs_Spliced.txt recount3_final.pass1.gff3 Hv3_splicedmasterTable_introns_recountTagged.gtf > Hv3_splicedmasterTable_allintronsSupported_recountTagged.gtf 
open(EC,$ARGV[0]);
open(ICs,$ARGV[1]);
open(RI,$ARGV[2]);
#open(IT,'>',$ARGV[3]);

while(<ICs>)
{	chomp;	#print $_;
	@Cinfo=split("\t",$_);	#print "$Cinfo[1]\n";
	$map_chains_anchs{$Cinfo[1]}=$Cinfo[0];	#anchTMs=chain
}
#-----make map for reount introns-----#
while(<RI>)
{	chomp;
	@rinfo=split("\t",$_);
	$recountINFO=$rinfo[0]."_".$rinfo[3]."_".$rinfo[4]."_".$rinfo[6];
	$RECOUNTintrons{$recountINFO}="$rinfo[8]";
}
while(<EC>)
{	chomp;	#print "$_\n";
	if($_=~"target")
	{	@ex=split("\t",$_);	
		@IDpre=split("\"",$_);
		$IDs=$IDpre[21]; 
		$transcriptInfo=$ex[8];	#print "$transcriptInfo\n";
		$currentChain=$map_chains_anchs{$IDs};##	print "###$IDs\t$map_chains_anchs{$IDs}!!!\n";
		@onlyChain=split("_",$currentChain);
		@intronChains=split(":",$onlyChain[1]);
		$totalchains=$#intronChains+1;
		$chainNum=0;$currenttrMatch=0;$allintronsdetails="";	#print "##$IDs\t$totalchains\n";
	#--sample info--#
		@samples=split(",",$IDpre[17]);
                $ont=0; $pacBio=0;
                foreach $S (@samples)
                {       @sinfo=split("",$S);
                        if($sinfo[7] eq "P"){$pacBio++;}	elsif($sinfo[7] eq "O"){$ont++;}
                }
                if($pacBio>0 && $ont>0)
                {       $tech="pacBio+ont";}elsif($pacBio>0 && $ont==0){$tech="pacBioOnly";}elsif($ont>0 && $pacBio==0){$tech="ontOnly";}
		$counter=0;
		foreach $intron(@intronChains)	#----intron level checks---#
		{	$counter++;
			@intronBounds=split(",",$intron);
			$intronstart=$intronBounds[0]+1;$intronstop=$intronBounds[1]-1;	#print"$IDs\t$intronstart\t$intronstop\n";
			$CLSintronINFO=$ex[0]."_".$intronstart."_".$intronstop."_".$ex[6];
	#--select an intron to be matched if it is present in the recount db; no other filters
			if(exists($RECOUNTintrons{$CLSintronINFO}))	
			{	$currenttrMatch++;
				$recountdetails=$RECOUNTintrons{$CLSintronINFO};
			}else{	$recountdetails="";}
	#print all current intron info TO GET INTRON TABLE
#			print IT "$ex[0]\t$ex[1]\tintron\t$intronstart\t$intronstop\t$ex[5]\t$ex[6]\t$ex[7]\t$transcriptInfo totalIntrons \"$totalchains\"; recount_details \"$recountdetails\";\n";  #matched_temp \"$currenttrMatch\";\n";	
			$allintronsdetails=$allintronsdetails."$ex[0]\t$ex[1]\tintron\t$intronstart\t$intronstop\t$ex[5]\t$ex[6]\t$ex[7]\t$transcriptInfo totalIntrons \"$totalchains\"; recount_details \"$recountdetails\"; tech \"$tech\";\n";
			if($counter == $totalchains)
			{	if($totalchains == $currenttrMatch)
				{	print "$_ totalIntrons \"$totalchains\"; tech \"$tech\"; transcriptRecountSupport \"supported\";\n$allintronsdetails";}
				else{	print "$_ totalIntrons \"$totalchains\"; tech \"$tech\"; transcriptRecountSupport \"unsupported\";\n$allintronsdetails";}
			}
		}	
	}
}
