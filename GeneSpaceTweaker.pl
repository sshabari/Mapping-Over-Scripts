#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use CoGeX;

##to connect to the right server
#my $coge;

my $connstr = 'dbi:mysql:dbname=coge;host=localhost;port=3307';
my $coge = CoGeX->connect($connstr, 'coge', '321coge123' );
#my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);

$|=1;


#Take every gene in at gene list
#1 remove all pseudogenes: Done
#2 remove all transposons: Done
#3 remove all M and C genes: Done
#4 remove genes that are uORFs
#5 only include the rna genes that diane wants to be included: Done
#6 MIR genes are ok to be included. Done.
#7 Remove antisense genes from the gene list: Done i think looking for overlaps
#8 Look in the at gene list to get local dups (tandems). Done.
#9 Look in haibaos list for a syntenic pair. Done.
#10 Possible tandems in Br identified. Done
#11 Is the ortholog a fragment. just make a note of this. Done
#12 Assign CNSs to At genes.
#14 cns_rna: Except for the genes that diane has called as good ones, reassign all other CNSs.Remake CNS list to indicate which ones are on valid cns-rna. Add valid cns_rna genes to gene hash. Then reassign CNSs based on new gene hash. DONE
#15 Reassign CNSs from uORFs to non-u-orf gene. The latest gene hash has the list of uORFs. See if any of the uORFs have CNSs assigned to them. If so, reassign the CNSs. DONE
#16 Reassign UTR and intron CNSs too. Done.
#16 Fix the switched prime problem with CNSs from Diane's proofing. Done. #17 Delete CNS assignments from GEvo and upload new assignments.
#18 Make CNS detectability sheet and proof. Add possible tandems and fragments in Br to this list.
#19 Make cartoon sheet and proof. Add possible Br tandems and fragments in Br to this list.
#20 rename proofing genome

my $type = "gene";
my $version;
my $hashfile=shift;
my $cnsfile=shift;
my $hash=do $hashfile;

#checkSynToNewAtGene($hash,$version,$cnsfile);
&findClosestGeneToCNS($hash,$cnsfile);

sub makeSheetWithGEvoLinks
  {
    my $file=shift;
    open(IN,"$file");
    <IN>;
    while(<IN>)
      {
	my $input=$_;
	chomp $input;
	my @array=split(/\t/,$input);
	my $link=$array[15];
	$link=~s/16746/19870/g;

	my $final=pop(@array);
	my $check=pop(@array);
	
	if ($check=~/diff/i)
	  {
	    if ($final=~/NA/)
	      {
		#Use same gevo link with 19870
	      }
	    elsif($final=~/old/i)
	      {
		#Use same gevo link with 19870
	      }
	    elsif($final=~/new/i)
	      {
		#Make new gevo link with 19870 and new At gene.
	      }
	  }
	elsif($check=~/same/i)
	  {
	    if($final=~/old/i)
	      {
		#Use same gevo link with 19870
	      }
	    elsif($final=~/new/i)
	      {
		#Use same gevo link with 19870
	      }
	  }
	
      }
  }


sub findClosestGeneToCNS
  {   
    my ($hash,$cnsfile)=@_;
    
    open(IN,"$cnsfile");      
    #<IN>;
  LINE: while(<IN>)
      {
	my $input=$_;
	chomp $input;
	print $input;
	
	my @array=split(/\t/,$input);	
	my $cns=$array[0];	  	  
	my $pos=$array[13];
	
	#print $cns,"\t",$pos,"\t";
	
	my @idarray=split(/\|/,$cns);
	my %distHash;
	my ($cgene,$cchr,$cstart,$cstop,$cstrand,$cblock)=($idarray[0],$idarray[1],$idarray[2],$idarray[3],$idarray[4],$idarray[5]);
	
	if ($array[1]=~/valid/)
	  {
	    print "\t","NA","\t","NA","\t","NA","\t","Same\n";
	    next LINE;
	  }
	
	my $mid=($cstop+$cstart)/2;
	
	foreach my $chr (sort keys %$hash)
	  {
	    if ($chr=~/^$cchr$/)
	      {		  
	      LINE2: foreach my $start(sort {$a<=>$b} keys %{$hash->{$chr}})
		  {    		      
		    my $strand=$hash->{$chr}{$start}{'strand'};
		    if ($strand=~/\-/)
		    {			  			        
		      next LINE2 unless exists $hash->{$chr}{$start}{'aa'}; 
		      my $aa=$hash->{$chr}{$start}{'aa'}; 
		      
		      my $gene=$hash->{$chr}{$start}{'gene'};   	       		
		      my $stop=$hash->{$chr}{$start}{'stop'};   
		      
		      if ($start>$stop)
			  {
			    my $start2=$start;			    
			    $start=$stop;
			    $stop=$start2;
			  }
		      
		      if ($pos=~/intron/i && $mid<=$stop && $mid>=$start)
			{
			  my $diff=0;
			  $distHash{$diff}{'gene'}=$gene;
			  $distHash{$diff}{'type'}="Intron";
			}
		      else
			{
			  my $diff1=abs($start-$mid);
			  $distHash{$diff1}{'gene'}=$gene;
			  $distHash{$diff1}{'strand'}=$cstrand;
			  
			  if ($diff1>=500)
			    {
				$distHash{$diff1}{'type'}="3-distal";
			      }
			    else
			      {
				$distHash{$diff1}{'type'}="3-proximal";
			      }
			    
			    my $diff2=abs($stop-$mid);
			    $distHash{$diff2}{'gene'}=$gene;
			    $distHash{$diff2}{'strand'}=$cstrand;
			    
			    if ($diff2>=500)
			      {
				$distHash{$diff2}{'type'}="5-distal";
			      }
			    else
			      {
				$distHash{$diff2}{'type'}="5-proximal";
			      }
			  }
		      }
		    else
		      {			
			next LINE2 unless exists $hash->{$chr}{$start}{'aa'};
			my $aa=$hash->{$chr}{$start}{'aa'}; 			
			my $stop=$hash->{$chr}{$start}{'stop'};   
			my $gene=$hash->{$chr}{$start}{'gene'};
						
			#next LINE2 unless $aa=~/$cblock/i;  
	       			
			if ($start>$stop)
			  {
			    my $start2=$start;			    
			    $start=$stop;
			    $stop=$start2;
			  }			
			
			if ($pos=~/intron/i && $mid<=$stop && $mid>=$start)
			  {
			    my $diff=0;			    
			    #print "I am in\n"; 
			    $distHash{$diff}{'gene'}=$gene;
			    $distHash{$diff}{'type'}="Intron";
			  }			
			else
			  {			   			    
			    my $diff1=abs($start-$mid);
			    $distHash{$diff1}{'gene'}=$gene;
			    $distHash{$diff1}{'strand'}=$cstrand;			
			    
			    if ($diff1>=500)
			      {
				$distHash{$diff1}{'type'}="5-distal";
			      }
			    else
			      {
				$distHash{$diff1}{'type'}="5-proximal";
			      }
			    
			    my $diff2=abs($stop-$mid);
			    $distHash{$diff2}{'gene'}=$gene;
			    $distHash{$diff2}{'strand'}=$cstrand;
			      
			    if ($diff2>=500)
			      {
				$distHash{$diff2}{'type'}="3-distal";
			      }
			    else
			      {
				$distHash{$diff2}{'type'}="3-proximal";
				}
			  }
		      }
		  }
	      }
	  }
	
	
	foreach my $key(sort {$a<=>$b} keys %distHash)
	  {
	    print "\t",$key,"\t",%distHash->{$key}{'gene'},"\t";
	    
	    my $pos2=%distHash->{$key}{'type'};
	    
	    print $pos2,"\t";
	    
	    if (%distHash->{$key}{'gene'}=~/$cgene/)
	      {
		print "Same\n";
		next LINE;
	      }
	    else
	      {
		print "Diff\n";
		next LINE;
	      }
	  }      
	print "\n";	
	}
    }
    
    sub checkSynToNewAtGene
      {    	
	my ($hash,$version,$cnsfile)=@_;
	
	open(IN,"$cnsfile");
	#<IN>;
    LINE: while (<IN>)
      {
	my $input=$_;
	chomp $input;
      	print $input,"\t";
	
	my @array=split(/\t/,$input);
	my ($cnsid,$aa,$aachr,$aastart,$aastop,$aastrand,$oldprime,$newat,$newprime)=($array[0],$array[6],$array[7],$array[8],$array[9],$array[10],$array[12],$array[17],$array[18]);
	
	if ($cnsid=~/cns_rna/i || $cnsid=~/cns_protein/i)
	  {
	    print "NA\n";
	    next LINE;
	  }

	if ($aa=~/19624b/)
	  {
	    print "NA\n";
	    next LINE;
	  }
	
	
	my ($chr,$start,$stop,$strand);
	
      LINE2: foreach my $feat_name ($coge->resultset('FeatureName')->search(
									    {name=>$aa}
									   ))
	  {
	    my $feat = $feat_name->feature;	
	    
	    next LINE2 unless $feat->type->name eq $type;	
	    #next LINE2 unless $feat->dataset->version eq $version;
	    $chr=$feat->chr;
	    $start=$feat->start;
	    $stop=$feat->stop;
	    $strand=$feat->strand;
	  }
	
	unless (defined $start)
	  {
	    print "NA\n";
	    next LINE;
	  }

	if ($start>$stop)
	  {
	    my $start2=$start;
	    $start=$stop;
	    $stop=$start2;
	  }
	my $prime;
	
	if ($strand=~/\-/)
	  {
	    if ($aastop<=$start)
	      {	
		$prime="3";		
	      }
	    elsif ($aastart>=$stop)
	      {
		$prime="5";
	      }
	    else
	      {
		$prime="Intron";
	      }
	  }
	else	
	  {
	    if ($aastop<=$start)
	      {	
		$prime="5";		
	      }
	    elsif ($aastart=>$stop)
	      {
		$prime="3";
	      }		
	    else
	      {
		$prime="Intron";
	      }
	  }
	
	if ($newprime=~/intron/i || $newprime=~/UTR/i)
	  {
	    print "Old","\n";
	    next LINE;
	  }
	
	elsif ($newprime=~/$prime/i)
	  {
	    print "New","\n";
	    next LINE;	    
	  }
	else
	  {
	    print "Old","\n";
	    next LINE;
	  }	
      }
  }


sub checkuORF
  {
    my ($cnsfile,$uorfile)=@_;
    
    open(IN,"$uorfile");
    <IN>;
    while(<IN>)
    {
      chomp;
      my @array=split(/\t/,$_);
      my $tag=$array[8];
      my $uorf=$1 if $tag=~/Parent\=(.*?);/;
      my $gene=$1 if $tag=~/major ORF .*?(AT.*?)\./i;
      my $grep=`grep -i $uorf "$cnsfile"`;
      if ($grep=~/At/i)
	{
	  print $uorf,"\t",$gene,"\n";
	}
    }
  }

sub sortLists
  {
    my ($file1,$file2)=@_;
    open(IN,"$file1");
    while(<IN>)
      {
	my $input=$_;
	chomp $input;
	my @array=split(/\t/,$input);
	my $id=$array[0];
	my $grep=`grep "$id" "$file2"`;
	if ($grep=~/At/i)
	  {
	    chomp $grep;
	    my @greparray=split(/\t/,$grep);
	    my $end=pop (@greparray);
	    my $diff=pop (@greparray);
	    my $pos=pop (@greparray);
	    my $gene=pop (@greparray);
	    my $dist=pop (@greparray);
	    
	    my $oldend=pop (@array);
	    my $olddiff=pop (@array);
	    my $oldpos=pop (@array);
	    my $oldgene=pop (@array);
	    my $olddist=pop (@array);
	    
	    print join("\t",@array,$dist,$gene,$pos,$diff,$end,"\n");
	  }
	else
	  {
	    print join("\t",@array,"\n");
	  }
      }
  }




sub cns_rna_fix_part2
  {
    my ($hash,$file)=@_;
    
    open(IN,"$file");
    while(<IN>)
      {
	my $input=$_;
	chomp $input;
	my @array=split(/\t/,$input);
	my $end=pop(@array);
	if ($end=~/to valid/)
	  {	 	 	    
	    my ($chr,$start,$stop,$cns,$rna)=split(/\_/,$array[1]);
	    	    
	    %$hash->{$chr}{$start}{'gene'}=$array[1];
	    %$hash->{$chr}{$start}{'stop'}=$stop;
	    %$hash->{$chr}{$start}{'strand'}="1";
	    %$hash->{$chr}{$start}{'uorf'}="No";
	    %$hash->{$chr}{$start}{'aa'}=$array[5];
	  }
      }
    return $hash;
  }

sub cns_rna_fix_part1
  { 
   my ($file,$file2)=@_;
   
   open(IN,"$file");
   while (<IN>)
     {
       #get gene name, grep into cns_rna file. are there genes in the cns_rna file that do not have a cns_rn notation in the cnsid?
       my $input=$_;       
       
       chomp $input;
              
       my @array=split(/\s+/,$input);
       my $cns=$array[0];
       
       if ($cns=~/cns_protein/i)
	 {
	   print $input,"\t";
	   	   
	   my @array2=split(/\|/,$cns);
	   my $id=$1 if $array2[0]=~/(.*?)_cns_protein/;
	   $id=~s/\_/\|/g;
	   
	   my $grep=`grep -i "$id" "$file2"`;
	   
	   if ($grep=~/\|/i)
	     {
	       print "Asssigned to valid rna gene.\n";
	     }
	   else
	     {
	       print "Assigned to invalid cns_protein gene.reassign\n";
	     }
	 }
       #else
       #{
       #print "NA\n";
       # }
     }
  }



#remove uORFs from the gene list
sub findUORFSandAddtoHash
  {    
    my $hashfile=shift;
    my $uorfile=shift;
    my $hash=do $hashfile;
    
    foreach my $chr (sort keys %$hash)
      {
	foreach my $start(sort {$a<=>$b} keys %{$hash->{$chr}})
	  {
	    my $gene=$hash->{$chr}{$start}{'gene'};	
	    my $grep=`grep "ID=$gene" "$uorfile"`;
	    if ($grep=~/At/i)
	      {
		$hash->{$chr}{$start}{'uorf'}="Yes";
	      }
	    else
	      {
		$hash->{$chr}{$start}{'uorf'}="No";
	      }
	  }
      }
    
    print Dumper $hash;
  }



  
sub findSyntenic
  {
    my ($hash,$cnsfile)=@_;

    foreach my $key (sort keys %$hash)
      {
	foreach my $key2(sort {$a<=>$b} keys %{$hash->{$key}})
	  {
	    my $gene=$hash->{$key}{$key2}{'gene'};
	    	    
	    my $grep=`grep "$gene" "$cnsfile"`;
	    if ($grep=~/At/i)
	      {
		foreach my $grepline(split(/\n/,$grep))
		  {
		    my @greparray=split(/\t/,$grepline);
		    
		    #for cns file
		    #$hash->{$key}{$key2}{'aa'}=$greparray[7]; 
		    
		    #for synmap file
		    next if exists $hash->{$key}{$key2}{'aa'};
		    my $aatag=$greparray[1];		    
		    my @tagarray=split(/\|\|/,$aatag);		
		    $hash->{$key}{$key2}{'aa'}=$tagarray[3];
		  }
	      }
	    else
	      {
		next;
	      }
	  }
      }
    return $hash;
  }

sub getOrthsIntoHash
  {
    my $hashfile=shift;
    my $file=shift;
    my $hash=do $hashfile;
    
    foreach my $key (sort keys %$hash)
      {
	foreach my $key2(sort {$a<=>$b} keys %{$hash->{$key}})
	  {
	    my $gene=$hash->{$key}{$key2}{'gene'};
	    my $grep=`grep "$gene" "$file"`;
	    chomp $grep;
	    my ($gene2,$b1,$b2,$b3)=split(/\t/,$grep);
	    $hash->{$key}{$key2}{'b1'}=$b1;
	    $hash->{$key}{$key2}{'b2'}=$b2;
	    $hash->{$key}{$key2}{'b3'}=$b3;
	    
	  }
      }
    print Dumper $hash;
    exit;
  }


sub getDups 
  {
    my $hashfile=shift;
    my $file=shift;
    my $hash=do $hashfile;
    
    foreach my $key (sort keys %$hash)
      {
	foreach my $key2(sort {$a<=>$b} keys %{$hash->{$key}})
	  {
	    my $gene=$hash->{$key}{$key2}{'gene'};
	    my $grep=`grep "$gene" "$file"`;
	    next unless  $grep=~/At/i;
	    foreach my $grepline(split(/\n/,$grep))
	      {
		chomp $grepline;
		my @array=split(/\t/,$grepline);
		if ($array[0]=~/$gene/)
		  {		    
		    if ($array[32]=~/mother/)
		      {
			$hash->{$key}{$key2}{'dup'}="Mother";
			$hash->{$key}{$key2}{'mother'}=$gene;
		      }
		    elsif($array[32]=~/At/i)
		      {
			$hash->{$key}{$key2}{'dup'}="Child";
			$hash->{$key}{$key2}{'mother'}=$array[32];		    
		      }
		    elsif($array[31]=~/I/)
		      {
			$hash->{$key}{$key2}{'dup'}="Interruptor";
		      }
		  }
	      }	
	  }
      }
    
    print Dumper $hash;
    exit;
  }


#open gene start stop file and make hash and print hash out with starts for + in start and start for minus is stop.

sub printhash
  {
    my %hash;
    my $file;
    open(IN,"$file");
    <IN>;
    while(<IN>)
      {
	chomp;
	my ($type,$chr,$start,$stop,$strand,$tag)=split(/\t/,$_);
	my @array=split(/\,/,$tag);
	my $gene=$array[0];
	if ($strand=~/^1$/)
	  {
	    $hash{$chr}{$start}{'gene'}=$gene;
	    $hash{$chr}{$start}{'strand'}=$strand;
	    $hash{$chr}{$start}{'stop'}=$stop;
	  }
	elsif ($strand=~/-1/)
	  {	
	    $hash{$chr}{$stop}{'stop'}=$start; 
	    $hash{$chr}{$stop}{'gene'}=$gene; 
	    $hash{$chr}{$stop}{'strand'}=$strand; 
	  }    
      }
    print Dumper (\%hash);
  }

