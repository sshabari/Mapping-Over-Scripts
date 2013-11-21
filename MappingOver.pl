#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use CoGeX;
use LWP::Simple;

my ($refgenomefile,$cnsfile,$orthfile,$output_dir,$qdsgid,$sdsgid,$qversion,$sversion);
##to connect to the right server
my $connstr = 'dbi:mysql:dbname=coge;host=localhost;port=3307';
my $coge = CoGeX->connect($connstr, 'coge', '123coge321' );


GetOptions ( 
	    "ref=s"            => \$refgenomefile,
	    "cns=s"            => \$cnsfile,
	    "orth=s"             => \$orthfile,
	    "out=s"             => \$output_dir,
	    "qdsgid=i"             => \$qdsgid,
	    "sdsgid=i"             => \$sdsgid,
	    #"qversion=i"           => \$qversion,
	    #"sversion=i"           => \$sversion,
	   );

$sversion="1.1"; 
my $min="29.5";
my $cnsminlength=0;


#input tairv10 file, cns file, ortholog file, output directory

#steps

#make output directory, make "hash" directory inside output directory to hold tairv10hash, cnshash and other hashes.

#system ('mkdir $output_dir/hash/');

#1. open tairv10 make hash of all genes and start stop, chr, strand. return hash. print hash out to "hash" directory in output directory.

#my $refhash=&getrefhash($refgenomefile); # this hash was made and printed out to hash/refhash. so no reason to run this again. open this hash in #2.
#open(OUT,">$output_dir/hash/refhash");
#print OUT Dumper $refhash;
#close OUT;

#2. pass the hash from #1 into subroutine along with cns file. open cns file and make two hashes, one "gene" hash of genes for analysis, with gene start, stop, chr, strand and one array for each prime of cns start stops. if prime does not exist, put in "NA". make "cns" hash of cns id, start, chr, stop, gene, strand. return "cns" hash. print "cns" hash into the "hash" directory.

#my $reffile="$output_dir/hash/refhash";
#my ($cnshash,$refhashwithgenespaces)=&getGeneSpaceBoundaries($cnsfile,$reffile); #takes a long time to run

#open(OUT,">$output_dir/hash/cnshash");
#print OUT Dumper $cnshash;
#close OUT;

#open(OUT2,">$output_dir/hash/refhashwithgenespaces");
#print OUT2 Dumper $refhashwithgenespaces;
#close OUT2;


#3. Open ortholog file and put into hash. return hash. print hash out to the "hash" directory.

#my $orthhash=&printorthhash($orthfile);
#open (OUT,">$output_dir/hash/orthhash");
#print OUT Dumper $orthhash;
#close OUT;


#4. send second hash from #2 and hash from #3 to another subroutine. for each gene in the second hash from #2, find orthologs from hash from #3. Add the ortholog names into second hash from #2. Return second hash from #2 now with ortholog info. print second hash from #2 into the "hash" directory.

#my $refhash_gs_orth=&findOrthologs();
#open(OUT,">$output_dir/hash/refhash_gs_orth");
#print OUT Dumper $refhash_gs_orth;
#close OUT;


#5. send hash from #4 subroutine to retrieve start stop positions from CoGe. For each At gene in hash, for each ortholog, get start stop positions fromm CoGe. put the start stops for each gene into the hash. return hash. print this hash into the "hash" directory.

#my $refhash_gs_orth_ss=&findrefstartstops($sversion); #can onlye be run in system with access to CoGe DB
#open(OUT,">$output_dir/hash/refhash_gs_orth_ss");
#print OUT Dumper $refhash_gs_orth_ss;
#close OUT;

#6. send hash from #5 to subroutine to generate start stops of orthologous windows for each prime without NA to be pulled out.return hash. print hash into the "hash" directory.

#my $refhash_orthwindows=&getrefhashwindows("$output_dir/hash/refhash_gs_orth_ss");
#open(OUT,">$output_dir/hash/refhash_windows");
#print OUT Dumper $refhash_orthwindows;
#close OUT;

#7. send hash to subroutine to pull out the sequence from CoGe for the orth window for each prime of each gene. open "search" directory in output directory. for each gene, for each prime, print out the fasta sequences into one file (named "gene\_prime") into the "search" directory.

#my $file="$output_dir/hash/refhash_windows";
#&getSearchSeqs($file);

#my $file=shift;
#&getSearchSeqs2($file); #nexw subroutine after cleaning up PL3.

#8. open cns hash file (#2) from "hash" directory, and send to subroutine to retrieve the sequence for each CNS from CoGe. make "query" directory in the output directory. print cns sequence in fasta format into one file (name cnsid). return hash.

#my $file="$output_dir/hash/cnshash";
#&getQuerySeqs($file);

#9 send hash from #8 into subroutine to blast each cns to the orthologous region. open hash, foreach cns, get the gene name and the prime for the cns. make "blastoutput" directory in output directoty. then run blast with query as "query" directory/cnsid and "search" directory/genename_prime. print output into cns file name under "blastoutput" directory. return hash.

#my $queryhashfile="$output_dir/hash/cnshash";
#my $orthhashfile="$output_dir/hash/orthhash";
#my $blastouthash=&blastAway($queryhashfile,$orthhashfile);
#open(OUT,">$output_dir/hash/cns_blastoutputhash");
#print OUT Dumper $blastouthash;
#close OUT;


#10 send hash from #9 to subroutine to parse blast under "blastoutput" directory for each cns name. parse out the bit score >$min. if so, then detectable, else no.

#my $file="$output_dir/hash/cns_blastoutputhash";
#my $parsedhash=&parseBlast($file,$min);
#open(OUT,">$output_dir/hash/cns_blastparsedhash_$min");
#print OUT Dumper $parsedhash;
#close OUT;


#11 Generate sorted cns order hash for each gene
#open cnshash, for each cns, get the position and gene, make hash of sorted cnss to gene, and subgenomed orthologs print out this hash.

#my $cnshashfile="$output_dir/hash/cnshash";
#my $refhash="$output_dir/hash/refhash_windows";
#my $sortedCNShash=&sortedCNShash($cnshashfile,$refhash);
#open(OUT,">$output_dir/hash/refhashwithcnsinfo");
#print OUT Dumper $sortedCNShash;
#close OUT;

#12 Generate cartoons for each gene.
#The At cartoon.  please plot each exon as a pipe (|) and start /stop as a $.   so a 4 exon gene might look like this: +++--++-+$-+|--|-|===| $ -+++


#my $refhashfile="$output_dir/hash/refhashwithcnsinfo";
#my $blastparsehashfile="$output_dir/hash/cns_blastparsedhash_$min";
#&makeCartoons($refhashfile,$blastparsehashfile,$cnsminlength);


#13 plot bitscore versus detectability. 
 
#my $blastparsedhashfile="$output_dir/hash/cns_blastparsedhash_$min";
#my $pointshash=&getPointsForDetGraph($blastparsedhashfile);

#foreach my $key(sort keys %$pointshash)
#{
#   my $count=$pointshash->{$key}{'singlet'}{'cnsnum'};
#   my $num=$pointshash->{$key}{'singlet'}{'det'};
#   my $percent=($num/$count)*100;
#   print $key,"\t","singlet","\t",$percent,"\t",$num,"\n";   

#   my $count2=$pointshash->{$key}{'doublet'}{'cnsnum'};
#   my $num2=$pointshash->{$key}{'doublet'}{'det'};
#   my $percent2=($num2/$count2)*100;
#   print $key,"\t","doublet","\t",$percent2,"\t",$num2,"\n";   

#   my $count3=$pointshash->{$key}{'triplet'}{'cnsnum'};
#   my $num3=$pointshash->{$key}{'triplet'}{'det'};
#   my $percent3=($num3/$count3)*100;
#   print $key,"\t","triplet","\t",$percent3,"\t",$num3,"\n";   
# }

#14 plot minus minus over bit score

#my $blastparsedhashfile="$output_dir/hash/cns_blastparsedhash_.$min";
#my $pointshash=&getPointsForMinuses($blastparsedhashfile);

#foreach my $key(sort keys %$pointshash)
#{
#  my $count=$pointshash->{$key}{'cnsnum'};
#  my $mm=$pointshash->{$key}{'mm'};
#  my $percent=($mm/$count)*100;
#  print $key,"\t",$percent,"\n"   
#}


#15 plot det over subgenomes and the bins.

#my $blastparsedhashfile="$output_dir/hash/cns_blastparsedhash";
#my $orthhashfile="$output_dir/hash/orthhash";
#my $pointshash=&getPointsForDetBySubG($blastparsedhashfile,$orthhashfile);
#print Dumper $pointshash;

#foreach my $key(sort keys %$pointshash)
#{
#    my $count=$pointshash->{$key}{'I'}{'num'};
#my $num=$pointshash->{$key}{'I'}{'det'};
#    my $percent=($num/$count)*100;
#print $key,"\t","I","\t",$percent,"\t",$num,"\n";   

#my $count2=$pointshash->{$key}{'II'}{'num'};
#my $num2=$pointshash->{$key}{'II'}{'det'};
#    my $percent2=($num2/$count2)*100;
#    print $key,"\t","II","\t",$percent2,"\t",$num2,"\n";   

#    my $count3=$pointshash->{$key}{'III'}{'num'};
#    my $num3=$pointshash->{$key}{'III'}{'det'};
#    my $percent3=($num3/$count3)*100;
#    print $key,"\t","III","\t",$percent3,"\t",$num3,"\n";   
#  }

#16 Filter CNSs by bit score
#my $blastparsedhashfile="$output_dir/hash/cns_blastparsedhash";
#&filterCNSsbyBit($blastparsedhashfile,$min);

#17 make Gevo Panels and Analysis sheet
#my $cartoonFile="$output_dir/results/cartoon_$min";
#my $cartoonhash=&getCartoonHash($cartoonFile);
#my $linkhash=&makeLinks($cartoonhash);
#&makeFinalSheet($linkhash);


#18 plot cns detectability to cns length
#my $blastparsedhashfile="$output_dir/hash/cns_blastparsedhash_.$min";
#my $pointshash=&getPointsForLength2DetGraph($blastparsedhashfile);

#foreach my $key(sort keys %$pointshash)
#{
#my $count=$pointshash->{$key}{'singlet'}{'cnsnum'};
#my $num=$pointshash->{$key}{'singlet'}{'det'};
#my $percent=($num/$count)*100;
#   print $key,"\t","singlet","\t",$percent,"\t",$num,"\n";   

#my $count2=$pointshash->{$key}{'doublet'}{'cnsnum'};
#my $num2=$pointshash->{$key}{'doublet'}{'det'};
#my $percent2=($num2/$count2)*100;
#print $key,"\t","doublet","\t",$percent2,"\t",$num2,"\n";   

#   my $count3=$pointshash->{$key}{'triplet'}{'cnsnum'};
#my $num3=$pointshash->{$key}{'triplet'}{'det'};
#   my $percent3=($num3/$count3)*100;
#   print $key,"\t","triplet","\t",$percent3,"\t",$num3,"\n";   
# }


#19 
#my $file="$output_dir/hash/cns_blastoutputhash";
#&parseBlastIntoCNSList($file,$min);

#20 GEvo hack
#&blastHack();



#########################

sub blastHack
  {
    #my $link="http://genomevolution.org/CoGe/GEvo.pl?prog=blastn;show_cns=1;iw=1000;fh=20;padding=2;nt=1;cbc=0;spike_len=15;skip_feat_overlap=1;skip_hsp_overlap=1;hs=0;bnW=7;bnG=5;bnE=2;bnq=-2;bnr=1;bne=30;bnF=F;accn1=AT1G01490;dsid1=80775;dsgid1=19870;dr1up=10000;dr1down=10000;ref1=1;accn2=Bra033260;dsid2=80774;dsgid2=19869;dr2up=30000;dr2down=30000;accn3=Bra032637;dsid3=80774;dsgid3=19869;dr3up=30000;dr3down=30000;num_seqs=4;hsp_overlap_limit=0;hsp_size_limit=0;autogo=1;";
    my $url = "http://genomevolution.org/CoGe/GEvo.pl?prog=blastn;show_cns=1;iw=1000;fh=20;padding=2;nt=1;cbc=0;spike_len=15;skip_feat_overlap=1;skip_hsp_overlap=1;hs=0;bnW=7;bnG=5;bnE=2;bnq=-2;bnr=1;bne=30;bnF=F;accn1=AT1G01490;dsid1=80775;dsgid1=19870;dr1up=10000;dr1down=10000;ref1=1;accn2=Bra033260;dsid2=80774;dsgid2=19869;dr2up=30000;dr2down=30000;accn3=Bra032637;dsid3=80774;dsgid3=19869;dr3up=30000;dr3down=30000;num_seqs=4;hsp_overlap_limit=0;hsp_size_limit=0;autogo=1;";

    my $browser = LWP::UserAgent->new;    
    my $response = $browser->get($url);
    
    my $html = $response->content;
    print $html;


  }

sub parseBlastIntoCNSList
  {
    my ($file,$min)=@_;
    my $hash=do $file || print STDERR "Cannot open $file";
    
    my %parsedHash;

    print "CNS id","\t","Br ortholog","\t","length of CNS","\t","Bit Score","\t","Threshold of detectability","\t","Detectability","\n";	    

    foreach my $key(sort keys %$hash)
      {

	my @splitarray=split(/\_/,$key);
	my ($cstart,$cstop)=($splitarray[2],$splitarray[3]);
	my $cnslength=($cstop-$cstart)+1;
	my $max=100+$cnslength;
	my $bit="<29.5";

	LINE: foreach my $outputfile(@{$hash->{$key}{'outfiles'}})
	  {
	    my @array1=split("\/",$outputfile);
	    my $tag=pop(@array1);
	    my @array2=split(/\_/,$tag);
	    my $br=pop(@array2);
	    	    
	    next unless -e $outputfile;
	    	    
	    my $det="N";	   	    
	    my $grep=`grep -v "#" "$outputfile"`;
	    	    
	    foreach my $grepline (split(/\n/,$grep))
	      {
		my @splitgrep=split(/\t/,$grepline);
		my ($cnsid,$qstart,$qstop);
		($cnsid,$qstart,$qstop,$bit)=($splitgrep[0],$splitgrep[6],$splitgrep[7],$splitgrep[11]);		
		#my @splitcns=split(/\_/,$cnsid);
		#my ($cstart,$cstop)=($splitcns[2],$splitcns[3]);		
		
		if ($bit>=$min)
		  {		    		    		    
		    unless ($qstop<100 || $qstart>$max)
		      {
			$det="Y";				      
			print $key,"\t",$br,"\t",$cnslength,"\t",$bit,"\t",$min,"\t",$det,"\n";	    
			next LINE;
		      }
		  }		  
	      }	    
	    print $key,"\t",$br,"\t",$cnslength,"\t",$bit,"\t",$min,"\t",$det,"\n";
	  }	
      }
  }


sub makeFinalSheet
  {
    my $hash=shift;
    foreach my $at(sort keys %$hash)
      {	
	my $link=$hash->{$at}{'tiny'} if (exists $hash->{$at}{'tiny'});
	my $dblchk=$hash->{$at}{'dblchk'} if (exists $hash->{$at}{'dblchk'});	
	foreach my $num(sort keys %{$hash->{$at}})
	  {
	    next unless $num=~/\d+/;
	    my $input=$hash->{$at}{$num}{'input'};
	    chomp $input;
	    
	    my $newcartoon="";

	    if (exists $hash->{$at}{$num}{'newcartoon'})
	      {
		$newcartoon=join("",@{$hash->{$at}{$num}{'newcartoon'}});
	      }
	    else
	      {
		$newcartoon=join("",@{$hash->{$at}{'newcartoon'}}) if exists $hash->{$at}{'newcartoon'};
	      }
		
	    
	    my @array=split(/\t/,$input);
	    
	    my $type=$array[6];
	    next unless $type=~/Doublet/i;
	    
	    #my $pluscount=0;
	    #my $minuscount=0;

	    #$pluscount = ($cartoon =~ tr/+/+/);
	    #$minuscount = ($cartoon =~ tr/-/-/);
	    
	    print $input,"\t",$newcartoon,"\t",$dblchk,"\t",$link,"\n";	
	   }
      }
  }

sub getPointsForLength2DetGraph
  {
    my $file=shift;

    my $hash=do $file || print STDERR "Cannot open $file";

    my %newhash;

    foreach my $key(sort keys %$hash)
      {
	my $bin="0";

	my @array=split(/\_/,$key);
	my ($start,$stop)=($array[2],$array[3]);
	my $cnsl=abs($stop-$start);
	$bin="20" if $cnsl<=20;
	$bin="30" if $cnsl>20 && $cnsl<=30;
	$bin="40" if $cnsl>30 && $cnsl<=40;
	$bin="50" if $cnsl>40 && $cnsl<=50;
	$bin="75" if $cnsl>50 && $cnsl<=75;
	$bin="100" if $cnsl>75 && $cnsl<=100;
	$bin="150" if $cnsl>100 && $cnsl<=150;
	$bin="200" if $cnsl>150 && $cnsl<=200;
	$bin=">200" if $cnsl>200;

	my $array2=$hash->{$key};
	
	my $total="NA";
	my $count=0;
	my $yes=0;

	foreach my $element(@$array2)
	  {	    
	    my ($br,$det)=split(/\_/,$element);	
	    
	    if ($br=~/Br/)
	      {
		++$count;
		if ($det=~/Y/)
		  {
		    ++$yes;
		  }
	      }
	  }
	$total="singlet" if $count==1;
	$total="doublet" if $count==2;
	$total="triplet" if $count==3;
	
	if (exists $newhash{$bin}{$total})
	  {
	    my $cnsnum=$newhash{$bin}{$total}{'cnsnum'};
	    ++$cnsnum;
	    $newhash{$bin}{$total}{'cnsnum'}=$cnsnum;
	    
	    if ($yes==$count)
	      {
		my $ynum=$newhash{$bin}{$total}{'det'};
		++$ynum;
		$newhash{$bin}{$total}{'det'}=$ynum;
	      }
	  }
	else
	  {
	    $newhash{$bin}{'singlet'}{'cnsnum'}=1;
	    $newhash{$bin}{'singlet'}{'det'}=0;
	    $newhash{$bin}{'doublet'}{'cnsnum'}=1;
	    $newhash{$bin}{'doublet'}{'det'}=0;
	    $newhash{$bin}{'triplet'}{'cnsnum'}=1;
	    $newhash{$bin}{'triplet'}{'det'}=0;


	    $newhash{$bin}{$total}{'cnsnum'}=1;
	    $newhash{$bin}{$total}{'det'}=1 if $yes==$count;	    
	  }
      }
    return (\%newhash);
  }

sub getCartoonHash
  {
    my $file=shift;
    my %hash;
    open(IN,"$file");
    while(<IN>)
      {
	my $input=$_;
	chomp $input;
	my @array=split(/\t/,$input);
	
	my ($at,$type,$subg,$br,$cartoon)=($array[0],$array[6],$array[10],$array[11], $array[15]);
	
	next if $cartoon=~/NA/;

	#get minus minus counts for doublets
	if ($type=~/doublet/i)
	  {	   
	    my $minusminuscount=0;
	    
	    if (exists $hash{$at}{'cartoons'})
	      {
		my $cartoon1=$hash{$at}{'cartoons'};		
		my @splitcartoon=split(//,$cartoon);
		
		my $i=0;
		
		foreach my $element(split(//,$cartoon1))
		  {
		    if ($element=~/\-/ && $splitcartoon[$i]=~/\-/)
		      {
			++$minusminuscount;
		      }
		    else
		      {	
			push (@{$hash{$at}{$subg}{'newcartoon'}},$splitcartoon[$i]);
			push (@{$hash{$at}{'newcartoon'}},$element);
		      }
		    ++$i;
		  }		
	      }	    
	    else
	      {
		$hash{$at}{'cartoons'}=$cartoon;		
	      }
	    $hash{$at}{'dblchk'}=$minusminuscount;
	  }
	$hash{$at}{$subg}{'input'}=$input;
      }
    return(\%hash);
  }


sub makeLinks
  {
    my $hash=shift;
    
    foreach my $at(sort keys %$hash)
      {	
	my $gevolink="http://genomevolution.org/CoGe/GEvo.pl?prog=blastn;show_cns=1;iw=1000;fh=20;padding=2;nt=1;cbc=0;spike_len=15;skip_feat_overlap=1;skip_hsp_overlap=1;hs=0;bnW=7;bnG=5;bnE=2;bnq=-2;bnr=1;bne=30;bnF=F;accn1=$at;dsid1=80775;dsgid1=19870;dr1up=20000;dr1down=20000;ref1=1;";     
	 
	my $count=2;	
	my $type;

	foreach my $num(sort keys %{$hash->{$at}})
	  {	    	    
	    ++$count;
	    next unless $num=~/\d+/;
	    my $input=$hash->{$at}{$num}{'input'};
	    chomp $input;
	    my @array=split(/\t/,$input);
	    
	    my $br=$array[11];
	    $type=$array[6];
	    my $aa=$array[1];
	    $gevolink.="accn2=$aa;dsid2=80498;dsgid2=19624;dr2up=20000;dr2down=20000;";
	    
	    next unless $type=~/Doublet/i;	    
	    next unless $br=~/Br/i; 
	    
	    my $accnum="accn".($count);
	    my $dsidnum="dsid".($count);
	    my $dsgidnum="dsgid".($count);
	    my $drupnum="dr".$count."up";
	    my $drdownnum="dr".$count."down";
	    
	    $gevolink.="$accnum=$br;$dsidnum=80774;$dsgidnum=19869;$drupnum=30000;$drdownnum=30000;";	   
	  }	
	next unless $type=~/Doublet/i;
	
	$gevolink.="num_seqs=4;hsp_overlap_limit=0;hsp_size_limit=0;autogo=1;";
	my $tiny = LWP::Simple::get("http://genomevolution.org/r/yourls-api.php?signature=d57f67d3d9&action=shorturl&format=simple&url=$gevolink");
	$tiny="=HYPERLINK(\"$tiny\",\"$at\")";	
	print STDERR $tiny,"\n";
	$hash->{$at}{'tiny'}=$tiny;
      }
    return($hash);
  }

  
#TBD
# how to resolve order?
# for undetectable CNSs, pull out the eapected orthologous regions, then do global alignment and calculate p-value comparing to 10,000 randomized runs.

sub filterCNSsbyBit
  {
    my ($hashfile,$min)=@_;
    
    my $hash=do $hashfile || print STDERR "Cannot open $hashfile";
    
    foreach my $key(sort keys %$hash)
      {
	my @array=split(/\_/,$key);
	my $bit=pop(@array);	
	next unless $bit>=$min;	
	print $key,"\n";
      }    
  }



sub getPointsForDetBySubG
  {
    my ($file1,$file2)=@_;
    
    my $hash=do $file1 || print STDERR "Cannot open $file1";
    my $orthhash=do $file2 || print STDERR "Cannot open $file2";
    my %newhash;
    foreach my $key(sort keys %$hash)
      {
	my $bin="29.5";

	my @array=split(/\_/,$key);
	my $bit=pop(@array);
	my $atgene=shift(@array);
	
	my $brline=$orthhash->{$atgene};
	
	my ($at,$br1,$br2,$br3)=split(/\t/,$brline);
	

	$bin="29.5" if $bit==29.5;
	$bin="30" if $bit==30;
	$bin="35" if $bit>30 && $bit<=35;
	$bin="40" if $bit>35 && $bit<=40;
	$bin="50" if $bit>40 && $bit<=50;
	$bin="75" if $bit>50 && $bit<=75;
	$bin="100" if $bit>75 && $bit<=100;
	$bin=">100" if $bit>100;
	
	my $array2=$hash->{$key};
		
	foreach my $element(@$array2)
	  {	    
	    my $subg;
	    my ($br,$det)=split(/\_/,$element);	
	    
	    if ($br=~/Br/)
	      {
		$subg="I" if $br=~/$br1/;
		$subg="II" if $br=~/$br2/;
		$subg="III" if $br=~/$br3/;
		
		
		if (exists $newhash{$bin}{$subg})
		  {
		    my $num=$newhash{$bin}{$subg}{'num'};
		    ++$num;
		    $newhash{$bin}{$subg}{'num'}=$num;
		    if ($det=~/Y/)
		      {
			my $det=$newhash{$bin}{$subg}{'det'};
			++$det;
			$newhash{$bin}{$subg}{'det'}=$det;
		      }
		  }
		else
		  {
		    $newhash{$bin}{$subg}{'num'}=1;
		    $newhash{$bin}{$subg}{'det'}=0;
		    if ($det=~/Y/)
		      {
			$newhash{$bin}{$subg}{'det'}=1;
		      }
		  }
	      }
	  }
      }
    return(\%newhash);
  }


sub getPointsForMinuses
  {
    my $file=shift;

    my $hash=do $file || print STDERR "Cannot open $file";

    my %newhash;

    foreach my $key(sort keys %$hash)
      {
	my $bin="29.5";

	my @array=split(/\_/,$key);
	my $bit=pop(@array);
	
	$bin="29.5" if $bit==29.5;
	$bin="30" if $bit==30;
	$bin="35" if $bit>30 && $bit<=35;
	$bin="40" if $bit>35 && $bit<=40;
	$bin="50" if $bit>40 && $bit<=50;
	$bin="75" if $bit>50 && $bit<=75;
	$bin="100" if $bit>75 && $bit<=100;
	$bin=">100" if $bit>100;
	
	my $array2=$hash->{$key};
	
	my $total="NA";
	my $count=0;
	my $no=0;

	foreach my $element(@$array2)
	  {	    
	    my ($br,$det)=split(/\_/,$element);	
	    
	    if ($br=~/Br/)
	      {
		++$count;
		if ($det=~/N/)
		  {
		    ++$no;
		  }
	      }
	  }
	if (exists $newhash{$bin})
	  {
	    my $cnsnum=$newhash{$bin}{'cnsnum'};
	    ++$cnsnum;
	    $newhash{$bin}{'cnsnum'}=$cnsnum;
	    
	    if ($no>0)
	      {
		my $mm=$newhash{$bin}{'mm'};
		++$mm;
		$newhash{$bin}{'mm'}=$mm;
	      }
	  }
	else
	  {
	    $newhash{$bin}{'cnsnum'}=1;
	    $newhash{$bin}{'mm'}=0;
	    $newhash{$bin}{'mm'}=1 if $no>0;	    
	  }
       
      }
    return (\%newhash);
  }


sub getPointsForDetGraph
  {
    my $file=shift;

    my $hash=do $file || print STDERR "Cannot open $file";

    my %newhash;

    foreach my $key(sort keys %$hash)
      {
	my $bin="29.5";

	my @array=split(/\_/,$key);
	my $bit=pop(@array);
	
	$bin="29.5" if $bit==29.5;
	$bin="30" if $bit==30;
	$bin="35" if $bit>30 && $bit<=35;
	$bin="40" if $bit>35 && $bit<=40;
	$bin="50" if $bit>40 && $bit<=50;
	$bin="75" if $bit>50 && $bit<=75;
	$bin="100" if $bit>75 && $bit<=100;
	$bin=">100" if $bit>100;
	
	my $array2=$hash->{$key};
	
	my $total="NA";
	my $count=0;
	my $yes=0;

	foreach my $element(@$array2)
	  {	    
	    my ($br,$det)=split(/\_/,$element);	
	    
	    if ($br=~/Br/)
	      {
		++$count;
		if ($det=~/Y/)
		  {
		    ++$yes;
		  }
	      }
	  }
	$total="singlet" if $count==1;
	$total="doublet" if $count==2;
	$total="triplet" if $count==3;
	
	if (exists $newhash{$bin}{$total})
	  {
	    my $cnsnum=$newhash{$bin}{$total}{'cnsnum'};
	    ++$cnsnum;
	    $newhash{$bin}{$total}{'cnsnum'}=$cnsnum;
	    
	    if ($yes==$count)
	      {
		my $ynum=$newhash{$bin}{$total}{'det'};
		++$ynum;
		$newhash{$bin}{$total}{'det'}=$ynum;
	      }
	  }
	else
	  {
	    $newhash{$bin}{'singlet'}{'cnsnum'}=1;
	    $newhash{$bin}{'singlet'}{'det'}=0;
	    $newhash{$bin}{'doublet'}{'cnsnum'}=1;
	    $newhash{$bin}{'doublet'}{'det'}=0;
	    $newhash{$bin}{'triplet'}{'cnsnum'}=1;
	    $newhash{$bin}{'triplet'}{'det'}=0;


	    $newhash{$bin}{$total}{'cnsnum'}=1;
	    $newhash{$bin}{$total}{'det'}=1 if $yes==$count;	    
	  }
      }
    return (\%newhash);
  }




sub makeCartoons
  {
    my ($reffile,$blastfile,$cnsmin)=@_;
    

    my $refhash=do $reffile || print STDERR "Cannot open $reffile";
    my $blasthash=do $blastfile || print STDERR "Cannot open $blastfile";


    foreach my $key(sort keys %$refhash)
      {	
	#my $cdscartoon;
	next unless exists $refhash->{$key}{'cds'};
	my $orthhash=$refhash->{$key}{'orths'};
	my $cds=$refhash->{$key}{'cds'};
       
	my $cdsnum=@$cds;
	my $i=1;
	
	#while ($i<=$cdsnum)
	# {
	#$cdscartoon.="\|";
	#   ++$i;
	# }
	my $aa=$refhash->{$key}{'aa'};
	my $str=$refhash->{$key}{'str'};
	my $chr=$refhash->{$key}{'chr'};
	my $start=$refhash->{$key}{'start'};
	my $stop=$refhash->{$key}{'stop'};

	my $orthcount=0;
	foreach my $key2(sort keys %$orthhash)
	  {		
	    my $br=$orthhash->{$key2}{'gene'};
	    if ($br=~/Br/)
	      {	
		++$orthcount;
	      }
	  }
	my $orthtype="No Br ortholog";
	$orthtype="Singlet" if $orthcount==1;
	$orthtype="Doublet" if $orthcount==2;
	$orthtype="Triplet" if $orthcount==3;
	
	
	foreach my $key2(sort keys %$orthhash)
	  {		
	    my $br=$orthhash->{$key2}{'gene'};
	    my $fnum=0;
	    my $tnum=0;
	    my $inum=0;
	   
	    if ($br=~/Br/)
	      {			  		    
		my $cartoon="";
		my $fcartoon="";
		my $tcartoon="";		
		my $icartoon="";		

		if (exists $refhash->{$key}{'5primecnss'})
		  {
		    my $fhash=$refhash->{$key}{'5primecnss'};
		    foreach my $key3(sort keys %$fhash)
		      {
			++$fnum;

			my $cnsid=$fhash->{$key3};		  
			$cnsid=~s/\|/\_/g;
			my @cnssplit=split(/\_/,$cnsid);
			my $bit=pop(@cnssplit);
		       	my ($start,$stop)=($cnssplit[2],$cnssplit[3]);
			my $length=abs($stop-$start);		
			
			next unless $length>=$cnsmin;

			if (exists $blasthash->{$cnsid})
			  {
			    my $array=$blasthash->{$cnsid};
			    
			    foreach my $element(@$array)
			      {
				my ($br2,$det)=split(/\_/,$element);
				my $call;
				$call="+" if $det=~/Y/;
				$call="-" if $det=~/N/;
				if ($br=~/$br2/)
				  {
				    $fcartoon.=$call;				    
				  }
			      }
			  }
		      }
		  }
		if (exists $refhash->{$key}{'3primecnss'})
		  {
		    my $thash=$refhash->{$key}{'3primecnss'};
		    
		    foreach my $key3(sort keys %$thash)
		      {
			++$tnum;
			my $cnsid=$thash->{$key3};		  
			$cnsid=~s/\|/\_/g;
			my @cnssplit=split(/\_/,$cnsid);
			my $bit=pop(@cnssplit);
		       	my ($start,$stop)=($cnssplit[2],$cnssplit[3]);
			my $length=abs($stop-$start);		
			
			next unless $length>=$cnsmin;

			if (exists $blasthash->{$cnsid})
			  {
			    my $array=$blasthash->{$cnsid};
			    
			    foreach my $element(@$array)
			      {
				my ($br2,$det)=split(/\_/,$element);
				my $call;
				$call="+" if $det=~/Y/;
				$call="-" if $det=~/N/;
				if ($br=~/$br2/)
				  {
				    $tcartoon.=$call;
				  }
			      }
			  }
		      }
		  }
		if (exists $refhash->{$key}{'introncnss'})
		  {
		    my $ihash=$refhash->{$key}{'introncnss'};
		    
		    foreach my $key3(sort keys %$ihash)
		      {
			++$inum;
			my $cnsid=$ihash->{$key3};		  
			$cnsid=~s/\|/\_/g;
			my @cnssplit=split(/\_/,$cnsid);
			my $bit=pop(@cnssplit);
		       	my ($start,$stop)=($cnssplit[2],$cnssplit[3]);
			my $length=abs($stop-$start);		
			
			next unless $length>=$cnsmin;
			
			if (exists $blasthash->{$cnsid})
			  {
			    my $array=$blasthash->{$cnsid};
			    
			    foreach my $element(@$array)
			      {
				my ($br2,$det)=split(/\_/,$element);
				$det="+" if $det=~/Y/;
				$det="-" if $det=~/N/;
				if ($br=~/$br2/)
				  {
				    $icartoon.=$det;
				  }
			      }
			  }
		      }
		  }
		if ($str=~/\+/)
		{
		  $cartoon.=$fcartoon;
		  $cartoon.="\$";
		   if ($cdsnum==1)
		      {
			$cartoon.="\|";			  
			$cartoon.=$icartoon;
		      }
		    elsif($cdsnum>1)
		      {	
			$cartoon.="\|";			  
			$cartoon.=$icartoon;
			$cartoon.="\|";			  
		      }
		  $cartoon.="\$";
		  $cartoon.=$tcartoon;
		}
		elsif($str=~/\-/)
		  {
		    $tcartoon=reverse($tcartoon);
		    $icartoon=reverse($icartoon);
		    $fcartoon=reverse($fcartoon);
		    
		    $cartoon.=$fcartoon;
		    $cartoon.="\$";
		    if ($cdsnum==1)
		      {
			$cartoon.="\|";			  
			$cartoon.=$icartoon;
		      }
		    elsif($cdsnum>1)
		      {	
			$cartoon.="\|";			  
			$cartoon.=$icartoon;
			$cartoon.="\|";			  
		      }
		    $cartoon.="\$";		    
		    $cartoon.=$tcartoon;
		  }
		print $key,"\t";
		print $aa,"\t";
		print $chr,"\t";
		print $start,"\t";
		print $stop,"\t";
		print $str,"\t";
		print $orthtype,"\t";
		print $fnum,"\t";
		print $tnum,"\t";
		print $inum,"\t";
		print $key2,"\t";
		print $br,"\t";
		print $fcartoon,"\t";
		print $icartoon,"\t";
		print $tcartoon,"\t";
		print $cartoon,"\n";
	      }
	    else
	      {
		print $key,"\t";
		print $aa,"\t";
		print $chr,"\t";
		print $start,"\t";
		print $stop,"\t";
		print $str,"\t";
		print $orthtype,"\t";
		print $fnum,"\t";
		print $tnum,"\t";
		print $inum,"\t";
		print $key2,"\t";
		print $br,"\t";
		print "NA","\t";
		print "NA","\t";
		print "NA","\t";
		print "NA","\n";
	      }
		
	  }		
      }	    
    return($refhash);
  }    



sub sortedCNShash
  {
    my ($cns,$ref)=@_;

    my $cnshash=do $cns || print STDERR "Cannot open $cns";
    my $refhash=do $ref || print STDERR "Cannot open $ref";
    
    foreach my $key(sort keys %$cnshash)
      {
	my $val=$cnshash->{$key};
	my @array=split(/\t/,$val);
	my ($cnsid,$gene,$start,$stop,$bit,$prime)=($array[0],$array[1],$array[3],$array[4],$array[11],$array[12]);
	print $cnsid,"\n";		
	print $gene,"\n";

	my @cnsidarray=split(/\|/,$cnsid);
	my $aagene=$cnsidarray[5];
	$refhash->{$gene}{'aa'}=$aagene;

	if (exists $refhash->{$gene})
	  {	    
	    if ($prime=~/5/)
	      {							
		if (exists $refhash->{$gene}{'5primecnss'})
		  {
		    my $hash=$refhash->{$gene}{'5primecnss'};
		    $hash->{$start}=$cnsid;
		    $refhash->{$gene}{'5primecnss'}=$hash;
		  }
		else
		  {
		    my $hash={};
		    $hash->{$start}=$cnsid;
		    $refhash->{$gene}{'5primecnss'}=$hash;
		  }
		  
	      }
	    if ($prime=~/3/)
	      {
		if (exists $refhash->{$gene}{'3primecnss'})
		  {
		    my $hash=$refhash->{$gene}{'3primecnss'};
		    $hash->{$start}=$cnsid;
		    $refhash->{$gene}{'3primecnss'}=$hash;
		  }
		else
		  {
		    my $hash={};
		    $hash->{$start}=$cnsid;
		    $refhash->{$gene}{'3primecnss'}=$hash;
		  }		
	      }
	    if ($prime=~/Intron/i)
	      {
		if (exists $refhash->{$gene}{'introncnss'})
		  {
		    my $hash=$refhash->{$gene}{'introncnss'};
		    $hash->{$start}=$cnsid;
		    $refhash->{$gene}{'introncnss'}=$hash;
		  }
		else
		  {
		    my $hash={};
		    $hash->{$start}=$cnsid;
		    $refhash->{$gene}{'introncnss'}=$hash;
		  }		
	      }
	  }
	else
	  {
	    next;
	  }
      }
    return($refhash);
  }

sub parseBlast
  {    
    my ($file,$min)=@_;
    my $hash=do $file || print STDERR "Cannot open $file";
    
    my %parsedHash;

    foreach my $key(sort keys %$hash)
      {
      LINE: foreach my $outputfile(@{$hash->{$key}{'outfiles'}})
	  {
	    my @array1=split("\/",$outputfile);
	    my $tag=pop(@array1);
	    my @array2=split(/\_/,$tag);
	    my $br=pop(@array2);
	    	    
	    next unless -e $outputfile;
	    	    
	    my $det="N";
	    
	    my $grep=`grep -v "#" "$outputfile"`;
	    	    
	    foreach my $grepline (split(/\n/,$grep))
	      {
		my @splitgrep=split(/\t/,$grepline);
		my ($cnsid,$qstart,$qstop,$bit)=($splitgrep[0],$splitgrep[6],$splitgrep[7],$splitgrep[11]);
			
		if ($bit>=$min)
		  {
		    my @splitcns=split(/\_/,$cnsid);
		    my ($cstart,$cstop)=($splitcns[2],$splitcns[3]);
		    my $cnslength=($cstop-$cstart)+1;
		    my $max=100+$cnslength;
		    		    		    
		    unless ($qstop<100 || $qstart>$max)
		      {
			$det="Y";			
			push (@{$parsedHash{$key}},$br."_".$det);
			
			next LINE;
		      }
		  }
	      }	
	    push (@{$parsedHash{$key}},$br."_".$det);
	  }	
      }
  return(\%parsedHash);
  }
                                                                                                               

sub blastAway
  {
    my ($file1,$file2)=@_;
    my $hash1=do $file1 || print STDERR "Cannot open $file1";
    my $hash2=do $file2 || print STDERR "Cannot open $file2";
    my %cnshash;


    foreach my $key(sort keys %$hash1)
      {
	my $cnsid=$key;
	$cnsid=~s/\|/_/g;
	my $queryfile="$output_dir/query/$cnsid";

	my $input=$hash1->{$key};
	
	my @array=split(/\t/,$input);
	my ($gene,$prime)=($array[1],$array[12]);
	$prime="5prime" if $prime=~/5/;
	$prime="3prime" if $prime=~/3/;
       	$prime="Intron" if $prime=~/intron/i;
	
	if (exists $hash2->{$gene})
	  {
	    my $orths=$hash2->{$gene};
	    foreach my $element(split(/\t/,$orths))
	      {
		next unless $element=~/Br/i;
		my $searchfile="$output_dir/search/$gene\_$element";
		my $outfile="$output_dir/blastoutputdir/$cnsid\_$element";
		
		push(@{$cnshash{$cnsid}{'outfiles'}},$outfile);
		
		next unless -e $searchfile;
		
		print "Blasting $cnsid to $element\n";

		#my $blast=`/home/elyons/src/ncbi-blast-2.2.25+/bin/legacy_blast.pl bl2seq -p blastn -o "$outfile" -i "$queryfile" -j "$searchfile" -D 1 -W 7 -G 5 -E 2 -q -2 -r 1 -e 1 -F F`;
		#my $blast=`/usr/local/bin/legacy_blast.pl bl2seq -p blastn -o "$outfile" -i "$queryfile" -j "$searchfile" -D 1 -W 7 -G 5 -E 2 -q -2 -r 1 -e 1 -F F`;
		my $blast=`/usr/local/bin/legacy_blast.pl bl2seq -p blastn -o "$outfile" -i "$queryfile" -j "$searchfile" -D 1 -W 7 -G 5 -E 2 -q -2 -r 1 -e 1 -F F`;
	      }
	  }
      }
    return(\%cnshash);
  }


sub getQuerySeqs
  {
    my $file=shift;
    my $hash=do $file || print STDERR "Cannot open $file";
    my $dsg = $coge->resultset('Genome')->find($qdsgid);
    
    foreach my $key(sort keys %$hash)
      {	
	my $cnsid=$key;
	$cnsid=~s/\|/\_/g;
	
	my $input=$hash->{$key};
	my @array=split(/\t/,$input);
	my ($chr,$start,$stop,$strand)=($array[2],$array[3],$array[4],$array[5]);
	$start=$start-100;
	$stop=$stop+100;
	my $seq;
	foreach my $ds($dsg->datasets)
	  {    
	    $seq=$ds->get_genomic_sequence(start=>$start, stop=>$stop, chr=>$chr);
	  }

	if ($strand=~/\-/)
	  {
	    $seq=~ tr/ACGTacgt/TGCAtgca/;
	    $seq=reverse($seq);
	  }	
	my $file=$cnsid;
	open(OUT,">$output_dir/query/$cnsid");
	print OUT "\>$cnsid\n";
	print OUT "$seq\n";
	close OUT;	
      }
  }
  

sub getSearchSeqs
  {
    my $file=shift;
    my $hash=do $file || print STDERR "Cannot open $file";
    my $dsg = $coge->resultset('Genome')->find($sdsgid);

    my %hash2;

    foreach my $ds($dsg->datasets)
      {    
	foreach my $chr($ds->chromosomes)
	  {
	    my $end = $ds->last_chromosome_position($chr);  
	    $hash2{$chr}{'end'}=$end;
	  }    
      }     
    
    foreach my $key(sort keys %$hash)
      {	
	if (exists $hash->{$key}{'orths'})
	  {
	    foreach my $key2(sort keys %{$hash->{$key}{'orths'}}) 
	      {	      		
		my $gene=$hash->{$key}{'orths'}{$key2}{'gene'};

		next unless $gene=~/Br/i;
		my $chr=$hash->{$key}{'orths'}{$key2}{'chr'};
		my $start=$hash->{$key}{'orths'}{$key2}{'start'};
		my $stop=$hash->{$key}{'orths'}{$key2}{'stop'};
		my $strand=$hash->{$key}{'orths'}{$key2}{'strand'};
		
		my $getstart;
		my $getstop;
		
		if ($start<=30000)
		  {
		    $getstart=1;
		  }
		else
		  {
		    $getstart=$start-30000; 
		  }
		$getstop=$stop+30000;
		$getstop=$hash2{$chr}{'end'} if ($getstop>=$hash2{$chr}{'end'});
		my $seq;
		
		foreach my $ds($dsg->datasets)
		  { 
		    $seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		  }
		
		if ($strand=~/\-/)
		  {
		    $seq=~ tr/ACGTacgt/TGCAtgca/;
		    $seq=reverse($seq);
		  }

		my $file=$key."_".$gene;
		open(OUT,">$output_dir/search/$file");
		print OUT "\>$key\_$gene\_$chr\_$getstart\_$getstop\_$strand\n";
		print OUT "$seq\n";
		close OUT;
		
		#if (exists $hash->{$key}{'5primewindow'})
		#{
		#my $fprimewindow=$hash->{$key}{'5primewindow'};
		#   $fprimewindow=30000 if $fprimewindow<30000;		
		#   if ($strand=~/-1/)
		#{
		#	my $getstart=$stop;
		#	my $getstop=$stop+$fprimewindow; #unless its greater than the chr length;
		#$getstop=$hash2{$chr}{'end'} if $getstop>=$hash2{$chr}{'end'};
		#my $seq;
		#foreach my $ds($dsg->datasets)
		#{    
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#}
		#my $file=$key."_".$gene."_5prime";
		#	open(OUT,">$output_dir/search/$file");
		#	print OUT "\>$key\_$gene\_5prime\n";
		#	print OUT "$seq\n";
		#	close OUT;
		#}
		#   else
		#     {
		#	my $getstop=$start;
		#my $getstart=$start-$fprimewindow;
		#	$getstart=1 if $start<=30000;
		#my $seq;
		#	foreach my $ds($dsg->datasets)
		#{
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#	  }
		#my $file=$key."_".$gene."_5prime";
		#	open(OUT,">$output_dir/search/$file");
		#	print OUT "\>$key\_$gene\_5prime\n";
		#	print OUT "$seq\n";
		#close OUT;
		#     }
		#  }
		#if (exists $hash->{$key}{'3primewindow'})
		#{
		#my $tprimewindow=$hash->{$key}{'3primewindow'};
		#$tprimewindow=30000 if $tprimewindow<30000;		
		#if ($strand=~/-1/)
		#{
		#	my $getstop=$start;
		#my $getstart=$start-$tprimewindow; 
		#	$getstart=1 if $start<=30000;
		#my $seq;
		#foreach my $ds($dsg->datasets)
		#{
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#}
		#	my $file=$key."_".$gene."_3prime";
		#	open(OUT,">$output_dir/search/$file");
		#print OUT "\>$key\_$gene\_3prime\n";
		#print OUT "$seq\n";
		#close OUT;
		#}
		#   else
		#     {
		#	my $getstart=$stop;
		#	my $getstop=$stop+$tprimewindow;
		#$getstop=$hash2{$chr}{'end'} if $getstop>=$hash2{$chr}{'end'};
		#my $seq;
		#	foreach my $ds($dsg->datasets)
		#	  {
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#}		
		#my $file=$key."_".$gene."_3prime";		    
		#open(OUT,">$output_dir/search/$file");
		#	print OUT "\>$key\_$gene\_3prime\n";
		#print OUT "$seq\n";		    
		#close OUT;
		#}
		#}
		
		#if (exists $hash->{$key}{'intronwindow'})
		#{
		#my $iprimewindow=$hash->{$key}{'intronwindow'};
		#my ($getstart,$getstop)=split(/\_/,$iprimewindow);
		    
		    
		#my $seq;
		#foreach my $ds($dsg->datasets)
		#{
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#}
		#my $file=$key."_".$gene."_Intron";
		#open(OUT,">$output_dir/search/$file");
		#print OUT "\>$key\_$gene\_Intron\n";
		#print OUT "$seq\n";
		#close OUT;
		#}		
	      }
	  }
      }
  }

sub getSearchSeqs2
  {
    my $file=shift;

    #open CNS file, foreach cns, find reassigned gene, for each reassigned At gene, find its iortholog in Br. Then get start stops for the Br gene, add 30k to each end and then retrieve that.
    
    my $dsg = $coge->resultset('Genome')->find($sdsgid);
    
    my %hash2;

    foreach my $ds($dsg->datasets)
      {    
	foreach my $chr($ds->chromosomes)
	  {
	    my $end = $ds->last_chromosome_position($chr);  
	    $hash2{$chr}{'end'}=$end;
	  }    
      }     
    
    foreach my $key(sort keys %$hash)
      {	
	if (exists $hash->{$key}{'orths'})
	  {
	    foreach my $key2(sort keys %{$hash->{$key}{'orths'}}) 
	      {	      		
		my $gene=$hash->{$key}{'orths'}{$key2}{'gene'};

		next unless $gene=~/Br/i;
		my $chr=$hash->{$key}{'orths'}{$key2}{'chr'};
		my $start=$hash->{$key}{'orths'}{$key2}{'start'};
		my $stop=$hash->{$key}{'orths'}{$key2}{'stop'};
		my $strand=$hash->{$key}{'orths'}{$key2}{'strand'};
		
		my $getstart;
		my $getstop;
		
		if ($start<=30000)
		  {
		    $getstart=1;
		  }
		else
		  {
		    $getstart=$start-30000; 
		  }
		$getstop=$stop+30000;
		$getstop=$hash2{$chr}{'end'} if ($getstop>=$hash2{$chr}{'end'});
		my $seq;
		
		foreach my $ds($dsg->datasets)
		  { 
		    $seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		  }
		
		if ($strand=~/\-/)
		  {
		    $seq=~ tr/ACGTacgt/TGCAtgca/;
		    $seq=reverse($seq);
		  }

		my $file=$key."_".$gene;
		open(OUT,">$output_dir/search/$file");
		print OUT "\>$key\_$gene\_$chr\_$getstart\_$getstop\_$strand\n";
		print OUT "$seq\n";
		close OUT;
		
		#if (exists $hash->{$key}{'5primewindow'})
		#{
		#my $fprimewindow=$hash->{$key}{'5primewindow'};
		#   $fprimewindow=30000 if $fprimewindow<30000;		
		#   if ($strand=~/-1/)
		#{
		#	my $getstart=$stop;
		#	my $getstop=$stop+$fprimewindow; #unless its greater than the chr length;
		#$getstop=$hash2{$chr}{'end'} if $getstop>=$hash2{$chr}{'end'};
		#my $seq;
		#foreach my $ds($dsg->datasets)
		#{    
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#}
		#my $file=$key."_".$gene."_5prime";
		#	open(OUT,">$output_dir/search/$file");
		#	print OUT "\>$key\_$gene\_5prime\n";
		#	print OUT "$seq\n";
		#	close OUT;
		#}
		#   else
		#     {
		#	my $getstop=$start;
		#my $getstart=$start-$fprimewindow;
		#	$getstart=1 if $start<=30000;
		#my $seq;
		#	foreach my $ds($dsg->datasets)
		#{
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#	  }
		#my $file=$key."_".$gene."_5prime";
		#	open(OUT,">$output_dir/search/$file");
		#	print OUT "\>$key\_$gene\_5prime\n";
		#	print OUT "$seq\n";
		#close OUT;
		#     }
		#  }
		#if (exists $hash->{$key}{'3primewindow'})
		#{
		#my $tprimewindow=$hash->{$key}{'3primewindow'};
		#$tprimewindow=30000 if $tprimewindow<30000;		
		#if ($strand=~/-1/)
		#{
		#	my $getstop=$start;
		#my $getstart=$start-$tprimewindow; 
		#	$getstart=1 if $start<=30000;
		#my $seq;
		#foreach my $ds($dsg->datasets)
		#{
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#}
		#	my $file=$key."_".$gene."_3prime";
		#	open(OUT,">$output_dir/search/$file");
		#print OUT "\>$key\_$gene\_3prime\n";
		#print OUT "$seq\n";
		#close OUT;
		#}
		#   else
		#     {
		#	my $getstart=$stop;
		#	my $getstop=$stop+$tprimewindow;
		#$getstop=$hash2{$chr}{'end'} if $getstop>=$hash2{$chr}{'end'};
		#my $seq;
		#	foreach my $ds($dsg->datasets)
		#	  {
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#}		
		#my $file=$key."_".$gene."_3prime";		    
		#open(OUT,">$output_dir/search/$file");
		#	print OUT "\>$key\_$gene\_3prime\n";
		#print OUT "$seq\n";		    
		#close OUT;
		#}
		#}
		
		#if (exists $hash->{$key}{'intronwindow'})
		#{
		#my $iprimewindow=$hash->{$key}{'intronwindow'};
		#my ($getstart,$getstop)=split(/\_/,$iprimewindow);
		    
		    
		#my $seq;
		#foreach my $ds($dsg->datasets)
		#{
		#$seq=$ds->get_genomic_sequence(start=>$getstart, stop=>$getstop, chr=>$chr);
		#}
		#my $file=$key."_".$gene."_Intron";
		#open(OUT,">$output_dir/search/$file");
		#print OUT "\>$key\_$gene\_Intron\n";
		#print OUT "$seq\n";
		#close OUT;
		#}		
	      }
	  }
      }
  }

sub getrefhashwindows
  {
    my $file=shift;
    my $hash=do $file || print STDERR "Cannot open $file";
    foreach my $k (sort keys %$hash)
      {	
	my $genestart=$hash->{$k}{'start'};
	my $genestop=$hash->{$k}{'stop'};
	
	if (exists $hash->{$k}{'5prime'})
	  {
	    my @array=@{$hash->{$k}{'5prime'}};
	    my @sortarray=sort(@array);
	    	    
	    if ($hash->{$k}{'str'}=~/\+/)
	      {
		my $bstart=shift(@sortarray);
		my $bstop=$hash->{$k}{'start'}-1;		
		my $diff=($bstop-$bstart)+1;
		$hash->{$k}{'5primewindow'}=$diff;
	      }
	    elsif($hash->{$k}{'str'}=~/\-/)
	      {
		my $bstart=$hash->{$k}{'stop'}+1;
		my $bstop=pop(@sortarray);
		my $diff=($bstop-$bstart)+1;
		$hash->{$k}{'5primewindow'}=$diff;
	      }
	  }
	if (exists $hash->{$k}{'3prime'})
	  {
	    my @array=@{$hash->{$k}{'3prime'}};
	    my @sortarray=sort(@array);
	    
	    if ($hash->{$k}{'str'}=~/\+/)
	      {
		my $bstart=$hash->{$k}{'stop'}+1;		
		my $bstop=pop(@sortarray);
		my $diff=($bstop-$bstart)+1;
		$hash->{$k}{'3primewindow'}=$diff;
	      }
	    elsif($hash->{$k}{'str'}=~/\-/)
	      {
		my $bstart=shift(@sortarray);
		my $bstop=$hash->{$k}{'start'}-1;
		my $diff=($bstop-$bstart)+1;
		$hash->{$k}{'3primewindow'}=$diff;
	      }
	  }
	if (exists $hash->{$k}{'Intron'})
	  {
	    my @array=@{$hash->{$k}{'Intron'}};
	    my @sortarray=sort(@array);
	    
	    $hash->{$k}{'intronwindow'}=$genestart."_".$genestop;
	  }
      }
    return($hash);
  }

sub getrefhash
  {
    #this is for AtTAIR 10 from Phytozome
    my $file=shift;
    my %hash;
    open(IN,"$file");
    <IN>;
    while(<IN>)
      {
	my $input=$_;
	chomp $input;
	if ($input=~/mRNA/i)
	  {
	    my @array=split(/\t/,$input);
	    my ($chrtag,$start,$stop,$strand)=($array[0],$array[3],$array[4],$array[6]);
	    my $chr=$1 if $chrtag=~/Chr(\d+)/; 
	    my $last=pop(@array);
	    my ($id,$gene)=($1,$2) if $last=~/pacid\=(.*?)\;.*?Parent\=(.*)/;
	    $hash{$id}{'gene'}=$gene;
	    $hash{$id}{'str'}=$strand;
	    $hash{$id}{'chr'}=$chr;
	    $hash{$id}{'start'}=$start;
	    $hash{$id}{'stop'}=$stop;
	    
	  }
	elsif($input=~/CDS/)
	  {
	    my @array=split(/\t/,$input);
	    my $last=pop(@array);
	    my $id=$1 if $last=~/pacid\=(.*)/;
	    
	    my ($cdsstart,$cdsstop)=($array[3],$array[4]);
	    
	    if (exists $hash{$id})
	      {
		push(@{$hash{$id}{'cds'}},$cdsstart."_".$cdsstop);		
		push(@{$hash{$id}{'cdsarray'}},$cdsstart,$cdsstop);		
	      }
	  }
      } 
    return(\%hash);
  }

sub getGeneSpaceBoundaries
  {
    my ($file,$reffile)=@_;
    my %cnshash;
    my %refhash2;
    
    my $refhash=do $reffile || print STDERR "Cannot open $reffile";

    open(IN,"$file");
    <IN>;
    while(<IN>)
      {
	my $input=$_;
	chomp $input;
		
	my @array=split(/\t/,$input);
	
	my ($id,$gene,$start,$stop,$prime)=($array[0],$array[1],$array[3],$array[4],$array[12]);
	
	print $id,"\n";

	$cnshash{$id}=$input;
	
	foreach my $key(sort keys %$refhash)
	  {
	    if ($gene=~/$refhash->{$key}{'gene'}/)
	      {
		$refhash2{$gene}{'str'}=$refhash->{$key}{'str'};
		$refhash2{$gene}{'cdsarray'}=$refhash->{$key}{'cdsarray'};
		$refhash2{$gene}{'cds'}=$refhash->{$key}{'cds'};
		$refhash2{$gene}{'chr'}=$refhash->{$key}{'chr'};
		$refhash2{$gene}{'start'}=$refhash->{$key}{'start'};
		$refhash2{$gene}{'stop'}=$refhash->{$key}{'stop'};
		
		
		push(@{$refhash2{$gene}{'5prime'}},$start,$stop) if $prime=~/5/;
		push(@{$refhash2{$gene}{'3prime'}},$start,$stop) if $prime=~/3/;
		push(@{$refhash2{$gene}{'Intron'}},$start,$stop) if $prime=~/int/i;
	      }
	  }
	
      }    
    
    close IN;

    return (\%cnshash,\%refhash2);

  }

sub printorthhash
  {
    my $file=shift;
    my %hash;

    open(IN,"$file");
    <IN>;
    while(<IN>)
      {
	my $input=$_;
	chomp $input;
	my @array=split(/\t/,$input);
	my $gene=$array[0];
	
	$hash{$gene}=$input;
	
      }
    close IN;
    return(\%hash);
  }
  
sub findOrthologs 
  {      
    my $reffilewithgenespace="$output_dir/hash/refhashwithgenespaces";    
    my $orthfile="$output_dir/hash/orthhash";
    
    my $refhashwithgenespace=do $reffilewithgenespace || print STDERR "Cannot open $reffilewithgenespace";
    my $orthhash=do $orthfile || print STDERR "Cannot open $orthfile";
    
    foreach my $key(sort keys %$refhashwithgenespace)
      {	
	if (exists $orthhash->{$key})
	  {
	    my @array=split(/\t/,$orthhash->{$key});
	    my $at=shift(@array);
	    my $i=1;
	    foreach my $element(@array)
	      {
		$refhashwithgenespace->{$key}{'orths'}{$i}{'gene'}=$element;
		++$i;
	      }	    
	  }
      }
    return($refhashwithgenespace);
  }
  
sub findrefstartstops
  {
    my $version=shift;
    my $reffile_gs_orth="$output_dir/hash/refhash_gs_orth";
    my $refhash_gs_orth=do $reffile_gs_orth || print STDERR "Cannot open $reffile_gs_orth";
    foreach my $key(sort keys %$refhash_gs_orth)
      {
	print $key,"\t";

	foreach my $key2(sort keys %{$refhash_gs_orth->{$key}{'orths'}}) 
	  {	      
	    my $orth=$refhash_gs_orth->{$key}{'orths'}{$key2}{'gene'};
	    print $orth,"\t";
	    if ($orth=~/Br/i)
	      {
		foreach my $feat_name ($coge->resultset('FeatureName')->search(
									       {name=>$orth}
									      ))
		  {		      
		    my $feat = $feat_name->feature;	
		    next if $feat->type->name eq "CDS";
		    next if $feat->type->name eq "mRNA";
		    next if $feat->type->name=~/RNA/i;
		    next unless $feat->dataset->version eq $version;
		    		    
		    $refhash_gs_orth->{$key}{'orths'}{$key2}{'chr'}=$feat->chr;
		    $refhash_gs_orth->{$key}{'orths'}{$key2}{'start'}=$feat->start;
		    $refhash_gs_orth->{$key}{'orths'}{$key2}{'stop'}=$feat->stop;
		    $refhash_gs_orth->{$key}{'orths'}{$key2}{'strand'}=$feat->strand;
		  }
	      }

	  }
	print "\n";
      }
    return ($refhash_gs_orth);
  }
