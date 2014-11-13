#!/usr/bin/perl
use warnings;
use strict;

my $temp = 'temp.txt';
open(my $FH,'>',$temp) or die "can't open file $temp $!\n";

#This script will take a CGH array input file and parse through to produce a hash of gene names and the corresponding patients and metrics (i.e. ABPOEC3G => chr22:..-.. => NA12891 => p-val=0.05 cnv=2)#
#This is a subroutine and will be transferred over to a module with other parsing subroutines.

    my $input = shift;
    my $filtered = shift;
    my $key = shift;
    my $cnv_cutoff = shift;

    unless ($cnv_cutoff){$cnv_cutoff =0}
open(IN, '<', $input) or die "can't open file $input $!\n";

    my @features;
    my %genes;
    my $sample;

#--------------------------------------------------------------------#
# Open the file and parse through each line.
#--------------------------------------------------------------------#

    while(my $line = <IN>){
	chomp $line;
	@features = split("\t",$line);

	my $id;
	my $chr;
	my $start;
	my $stop;
	my $probeset;
	my $cnv;
	my @multiple;

#--------------------------------------------------------------------#
# Look for the sample name i.e. NA12891.
#--------------------------------------------------------------------#

	if($features[0] =~ /US_/){
	    $sample = $features[0];
	}
    
	elsif(($features[8]<0.05) && ($features[9])){
	    $id = $features[9];
	    $id = uc($id);
	        
#--------------------------------------------------------------------#
# If the gene list has multiple genes separated by commas, this
# separates the criteria per each gene.
#--------------------------------------------------------------------#

	    if($id =~ /.+,\s.+/){
		    
		@multiple = split(/, /,$id);
		    
		foreach my $individualgene (@multiple){
		    $chr = $features[1];
		    $start = $features[3];
		    $stop = $features[4];
		    $probeset = $chr . ':' . $start . '-' . $stop;
		    
		    if($cnv_cutoff == 0){
			if($features[6]>$cnv_cutoff){
			    $cnv = $features[6];
			}
			elsif($features[7]<$cnv_cutoff){
			    $cnv = $features[7];
			}       
		    }

		    else{
			if ($cnv_cutoff > 0){ 
			    if($features[6]>$cnv_cutoff){
				$cnv = $features[6];
			    }
			    else{next;}
			}
			if ($cnv_cutoff < 0){   
			    if($features[7]<$cnv_cutoff){
				$cnv = $features[7];
			    }
			    else{next;}
			}}
		    $genes{$individualgene}{$probeset} = $sample;
#		    $genes{$individualgene}{$probeset}{$sample;
		    
		}
	    }
	        
#--------------------------------------------------------------------#
# If the gene list has only ONE gene, populate the hash with the 
# criteria from that line.
#--------------------------------------------------------------------# 

	    else{
		$chr = $features[1];
		$start = $features[3];
		$stop = $features[4];
		$probeset = $chr . ':' . $start . '-' . $stop;

		if($cnv_cutoff == 0){
		    if($features[6]>$cnv_cutoff){
			$cnv = $features[6];
		    }
		    elsif($features[7]<$cnv_cutoff){
			$cnv = $features[7];
		    }       
		}
		else{
		    if ($cnv_cutoff > 0){ 
			if($features[6]>$cnv_cutoff){
			    $cnv = $features[6];
			}
			else{next;}
		    }
		    if ($cnv_cutoff < 0){   
			if($features[7]<$cnv_cutoff){
			    $cnv = $features[7];
			}
			else{next;}
		    }
		}
		$genes{$id}{$probeset} = $sample;
	#	$genes{$id}{$probeset}{$sample}{'CNV'} = $cnv;
	    }
	}

    }
#    return %genes;

open(MATCH,'<',$filtered) or die "can't open file $filtered $!\n";

print $FH 'Gene',"\t",'Probeset',"\t",'Sample',"\t",'log2',"\t",'PVAL',"\t",'EXONTRUE',"\n";

my @feat_filt;
#my $sample;

#-----------------------------------------------------------------------------------------------------#
# This portion of the script takes the file based on Omer's stringencies (yet without sample names)   #
# and matches the probeset to that in the parsed unfiltered data occupying %genes. Result is to print #
# out a file with the header and the gene name with the correspnding sample that matches and the info #
# from the stringent file.                                                                            #
#-----------------------------------------------------------------------------------------------------#

while(my $filt_line = <MATCH>){
    chomp $filt_line;
    @feat_filt = split("\t",$filt_line);

    my $id_filt;
    my $chr_filt;
    my $start_filt;
    my $stop_filt;
    my $probeset_filt;
    my $log2;
    my $PVAL;
    my $EXONTRUE;	
    my $sample_match;

	if($feat_filt[6] !~ /NA/){
        
	    $chr_filt = $feat_filt[0];
	    $start_filt = $feat_filt[1];
	    $stop_filt = $feat_filt[2];
	    $id_filt = $feat_filt[6];
	    $log2 = $feat_filt[4];
	    $PVAL = $feat_filt[5];
	    $EXONTRUE = $feat_filt[9];
	    $probeset_filt = $chr_filt . ':' . $start_filt . '-' . $stop_filt;
        
#	    $sample_match = $genes{$id_filt}{$probeset_filt};
	    $sample_match = $genes{$key}{$probeset_filt};

	    if($sample_match){
	    print $FH "$key","\t","$probeset_filt","\t","$sample_match","\t","$log2","\t","$PVAL","\t","$EXONTRUE\n"; 
	    }
	    
	    else{
		next;
	    }
       }

    
        else{
	    next;
	}
}

#print "temp filename: $tmp_fh\n";
