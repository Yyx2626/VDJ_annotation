#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2021-12-21)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: cat <HTGTS_VDJ_annotated.tsv> | $0 <output_prefix> <ref_D_usage.tsv>
	[allowed_possible_D_num (default: 99)] [bait_IGHJ]
Input:
	<ref_D_usage.tsv>	two columns: D, usage(%)
Output:
	<output_prefix>.add_DinVDJ.tsv
	<output_prefix>.VH_DinVDJ.tsv
".$version;

if(@ARGV < 2){
	die $usage;
}
my $output_prefix = shift(@ARGV);

my $D_usage_pseudocount = 0.000;
my $ref_D_usage_filename = undef;
my @D_vec = ();
my %D2idx = ();
my %D_usage = ();
my ($now_D, $now_usage);
my $total_usage = 0;
my (@F);
if(@ARGV > 0){
	$ref_D_usage_filename = shift(@ARGV);
	open(IN, $ref_D_usage_filename) or die "Error: cannot open ref_D_usage $ref_D_usage_filename for input\n";
	while(<IN>){
		if(/^\s*$/){  next;  }
		s/[\r\n]+$//;
		@F = split/\t/;
		$now_usage = $F[1] + $D_usage_pseudocount;
		$D_usage{$F[0]} = $now_usage;
		$total_usage += $now_usage;
		if($F[0] ne "-"){
			$D2idx{$F[0]} = scalar(@D_vec);
			push(@D_vec, $F[0]);
		}
	}
	close(IN);
	
	foreach $now_D (keys %D_usage){
		$D_usage{$now_D} /= $total_usage;
	}
}
my $D_num = scalar(@D_vec);

my $allowed_possible_D_num = 99;
if(@ARGV > 0){
	$allowed_possible_D_num = shift(@ARGV);
}

my $bait_IGHJ = "IGHJ";
if(@ARGV > 0){
	$bait_IGHJ = shift(@ARGV);
}


my @V_vec = ();
my %V_D_mat = ();


my $output_filename = $output_prefix . ".add_DinVDJ.tsv";
open(OUT, ">" . $output_filename) or die "Error: cannot open $output_filename for output\n";
my $headline = <STDIN>;
$headline =~ s/[\r\n]+$//;
my @fields = split(/\t/, $headline);
print OUT join("\t", @fields, @D_vec)."\n";
my %field2idx = ();
my ($i, $k);
for($i=0; $i<@fields; $i++){
	$field2idx{$fields[$i]} = $i;
}
my ($mid_D_annotate_str, @possible_mid_Ds, @possible_D_prob_vec, $now_sum, @G, $now_idx, $now_V);
#my $total_reads = 0;
while(<STDIN>){
	s/[\r\n]+$//;
	@F = split/\t/;
	@G = (0) x $D_num;

	if($F[$field2idx{"bait_overlap_features"}] =~ /^$bait_IGHJ/i){   # bait = IGHJ
		$now_V = $F[$field2idx{"junction_overlap_features"}];
		if($now_V =~ /^IGHV/i){   # junction = IGHV
			if(!exists($V_D_mat{$now_V})){
				$V_D_mat{$now_V} = [(0) x $D_num];
				push(@V_vec, $now_V);
			}
			$mid_D_annotate_str = $F[$field2idx{"mid_D_annotate"}];
			@possible_mid_Ds = unique(map{  s/\s*\(rC\)//g; $_;  } split(/, /, $mid_D_annotate_str));
#			$total_reads++;
			if($mid_D_annotate_str eq "-" || @possible_mid_Ds > $allowed_possible_D_num){
#				$mid_D_annotate_str = "-";
#				$sum_D_reads{"-"} ++;
			}else{
				@possible_D_prob_vec = (1/scalar(@possible_mid_Ds)) x @possible_mid_Ds;   # evenly assigned, if no ref_D_usage is proviced
				if(defined($ref_D_usage_filename)){
					@possible_D_prob_vec = @D_usage{@possible_mid_Ds};
					$now_sum = sum(@possible_D_prob_vec);
					@possible_D_prob_vec = map {  $_ / $now_sum;  } @possible_D_prob_vec;
				}
				
#				print STDERR join("  ", @possible_mid_Ds)."\n";
#				print STDERR join("  ", @possible_D_prob_vec)."\n";
				
				for($k=0; $k<@possible_mid_Ds; $k++){
					$now_idx = $D2idx{$possible_mid_Ds[$k]};
					$G[$now_idx] += $possible_D_prob_vec[$k];
					$V_D_mat{$now_V}->[$now_idx] += $possible_D_prob_vec[$k];
				}
			}
		}
	}
	print OUT join("\t", @F, @G)."\n";
}
close(OUT);

#print STDERR $total_reads."\n";
#if($total_reads < 0.0001){
#	$total_reads = 0.0001
#}
$output_filename = $output_prefix . ".VH_DinVDJ.tsv";
open(OUT, ">" . $output_filename) or die "Error: cannot open $output_filename for output\n";
print OUT join("\t", "VH", @D_vec)."\n";
foreach $now_V (@V_vec){
	print OUT join("\t", $now_V, @{$V_D_mat{$now_V}})."\n";
}
close(OUT);


0;



sub unique{
	my %S = ();
	my @ans = ();
	foreach (@_){
		if(!exists($S{$_})){
			$S{$_} = 1;
			push(@ans, $_);
		}
	}
	return(@ans);
}

sub sum{
	my $ans = 0;
	foreach (@_){
		$ans += $_;
	}
	return($ans);
}

sub setdiff{
	my ($set1, $set2) = @_;
	my %S = ();
	my @ans = ();
	foreach (@$set1){
		$S{$_} = 1;
	}
	foreach (@$set2){
		$S{$_} = 0;
	}
	foreach (sort keys %S){
		if($S{$_} == 1){
			push(@ans, $_);
		}
	}
	return @ans;
}




