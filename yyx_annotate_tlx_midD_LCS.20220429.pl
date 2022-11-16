#!/usr/bin/perl

### ref: yyx_annotate_tlx_midD_LCS.20190116.py

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2022-04-29)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "cat <input.tlx> | perl $0 <D.fa> [score_cutoff (default:5)]
Input: STDIN   <input.tlx>
    should have at least these columns:
	B_Qstart, B_Qend, Qstart, Qend, Seq
    or
	mid
Output: STDOUT   append (5)+2 columns
	pre  bait  mid  prey  post
	mid_D_score  mid_D_annotate
".$version;


## 2018-12-23 turn back to use continuous_LCS_DP (LCS=longest continuous substring, not allow gap) instead of Bio.pairwise2.align.globalms
sub continuous_LCS_DP{
	my ($query, $subjt, $caseSensitive) = @_;
	if(!defined($caseSensitive)){
		 $caseSensitive = 0;
	}
	if(! $caseSensitive){
		$query =~ tr/a-z/A-Z/;
		$subjt =~ tr/a-z/A-Z/;
	}
	my @Q = split(//, $query);
	my @S = split(//, $subjt);
	my $len1 = scalar(@Q);
	my $len2 = scalar(@S);
	my $max_LCS_len = 0;
	my ($max_len_i, $max_len_j) = (-1, -1);
	my ($i, $j);
#	DPmat = np.zeros((len1+1, len2+1), dtype='int')
	my @DPmat = ([ (0) x ($len2+1) ]) x ($len1+1);
	for($i=1; $i<=$len1; $i++){
		for($j=1; $j<=$len2; $j++){
			if($Q[$i-1] eq $S[$j-1]){
				$DPmat[$i]->[$j] = $DPmat[$i-1]->[$j-1] + 1;
				if($DPmat[$i]->[$j] > $max_LCS_len){
					$max_LCS_len = $DPmat[$i]->[$j];
					($max_len_i, $max_len_j) = ($i, $j);
				}
			}
		}
	}
	return ($max_LCS_len, $max_len_i, $max_len_j);
	# max_match_shift = j - i
}


#trailing_branket_pattern = re.compile('\(.*$')


if(@ARGV < 1){
	die $usage;
}

my $D_ref_fa_filename = shift(@ARGV);

my $score_cutoff = 5;
if(@ARGV > 0){
	$score_cutoff = shift(@ARGV);
}


### read in D reference fasta file

sub read_fasta{
	my ($input_fa_filename) = @_;
	my $title = undef;
	my $seq = "";
	my $ans = {};
	open(IN, $input_fa_filename) or die "Error: cannot open fasta file $input_fa_filename for input\n";
	print STDERR "Now reading fasta file $input_fa_filename ...\n";
	while(<IN>){
		s/[\r\n]+$//;
		if(/^>(.*?)\s*$/){
			if(defined($title)){
				if(exists($ans->{$title})){
					print STDERR "Warning: duplicated sequence title $title in fasta file $input_fa_filename\n";
				}
				$ans->{$title} = $seq;
			}
			$title = $1;
			$seq = "";
		}else{
			$seq .= $_;
		}
	}
	if(defined($title)){
		if(exists($ans->{$title})){
			print STDERR "Warning: duplicated sequence title $title in fasta file $input_fa_filename\n";
		}
		$ans->{$title} = $seq;
	}
	close(IN);
	return $ans;
}

my @four_bases = qw/ A C G T /;
my $i;
my %complement_base_hash = ();
my ($b1, $b2);
for($i=0; $i<=2; $i++){
	$b1 = $four_bases[$i];
	$b2 = $four_bases[3-$i];
	if($i==2){
		$b1 = "N";
		$b2 = "N";
	}
	$complement_base_hash{$b1} = $b2;
	$complement_base_hash{$b2} = $b1;
	$b1 =~ tr/A-Z/a-z/;
	$b2 =~ tr/A-Z/a-z/;
	$complement_base_hash{$b1} = $b2;
	$complement_base_hash{$b2} = $b1;
}
sub hash_get_default{
	my ($hash, $query, $default, $error_message) = @_;
	if(exists($hash->{$query})){
		return $hash->{$query};
	}
	if(defined($error_message)){
		print STDERR $error_message;
	}
	return $default;
}
sub reverse_complement{
	my @S;
	my @ans = map{
		join("", map{
			hash_get_default(\%complement_base_hash, $_, $_);
		} reverse split// );
	} @_;
	if(@ans==1){
		return $ans[0];
	}
	return @ans;
}


my $D_raw_hash = read_fasta($D_ref_fa_filename);
my $D_ref_hash = {};
my $D_rC_hash = {};   # reverseComplement
my ($now_ID, $seq);
foreach $now_ID (keys %$D_raw_hash){
	$seq = $D_raw_hash->{$now_ID};
	$now_ID =~ s/lcl[|]//;
	$D_ref_hash->{$now_ID} = $seq;
	$D_rC_hash->{$now_ID} = reverse_complement($seq);
}

sub best_align{
	my ($query, $ref_hash) = @_;
	my $max_score = -1;
	my $max_score_hits = [];
	my ($k, $v, @alignment, $score);
	foreach $k (keys %$ref_hash){
		$v = $ref_hash->{$k};
		@alignment = continuous_LCS_DP($query, $v);
		$score = $alignment[0];
		if($score > $max_score){
			$max_score = $score;
			$max_score_hits = [$k];
		}elsif(abs($score-$max_score) < 1e-5){
			push(@$max_score_hits, $k);
		}
	}
	return ($max_score, $max_score_hits);
}

sub min{
	my $ans = undef;
	foreach (@_){
		if(!defined($ans)){
			$ans = $_;
		}elsif($_ < $ans){
			$ans = $_;
		}
	}
	return $ans;
}
sub max{
	my $ans = undef;
	foreach (@_){
		if(!defined($ans)){
			$ans = $_;
		}elsif($_ > $ans){
			$ans = $_;
		}
	}
	return $ans;
}

my %Dalign_storage = ();
sub Dalign{
	my ($query) = @_;
	if(exists($Dalign_storage{$query})){
		return @{$Dalign_storage{$query}};
	}
	my ($D_score, $D_hit) = best_align($query, $D_ref_hash);
	my ($D_rC_score, $D_rC_hit) = best_align($query, $D_rC_hash);
	my $max_score = max($D_score, $D_rC_score);
	my $max_score_hits = [];
	if(abs($D_score - $max_score) < 1e-5){
		push(@$max_score_hits, sort @$D_hit);
	}
	if(abs($D_rC_score - $max_score) < 1e-5){
		push(@$max_score_hits, map { $_." (rC)"; } sort @$D_hit);
	}
	$Dalign_storage{$query} = [$max_score, $max_score_hits];
	return ($max_score, $max_score_hits);
}

my $is_headline = 1;
my @fields = ();
my %f2i = ();
my $NR = 0;
my (@output_fields, @F, $now_colname);
my ($B_Qstart, $B_Qend, $P_Qstart, $P_Qend, $prey, $bait, $mid, $pre, $post);
my ($P_Rname, $P_Rstart, $P_Rend, $P_Strand, $B_Rname, $B_Rstart, $B_Rend, $B_Strand, $Rbait, $Rprey);
my ($mid_annotate, $score);
while(<STDIN>){
	$NR++;
	if($NR % 10000 == 0){
		print STDERR "Now process ".$NR."-th read ...\n";
	}
	s/[\r\n]+$//;
	@F = split/\t/;
	if($is_headline){
		$is_headline = 0;
		@fields = @F;
		for($i=0; $i<@F; $i++){
			$f2i{$F[$i]} = $i;
		}
		if(!exists($f2i{"mid"})){
			foreach $now_colname (qw/B_Qstart B_Qend Qstart Qend Seq/){
				if(!exists($f2i{$now_colname})){
					die "Error: cannot find colname '$now_colname'\n";
				}
			}
		}
		@output_fields = @fields;
		if(!exists($f2i{"mid"})){
			push(@output_fields, qw/pre bait mid prey post/);
		}
		push(@output_fields, qw/mid_D_score mid_D_annotate/);
		print join("\t", @output_fields)."\n";
		next;
	}else{
		if(!exists($f2i{"mid"})){
			$seq = $F[$f2i{"Seq"}];
			$B_Qstart = $F[$f2i{"B_Qstart"}];
			$B_Qend = $F[$f2i{"B_Qend"}];
			$P_Qstart = $F[$f2i{"Qstart"}];
			$P_Qend = $F[$f2i{"Qend"}];
			if($P_Qstart < $B_Qstart){
				print STDERR "Warning: prey start < bait start for Line $NR\n";
			}
			$bait = substr($seq, $B_Qstart-1, $B_Qend-($B_Qstart-1));
			$prey = substr($seq, $P_Qstart-1, $P_Qend-($P_Qstart-1));
			$pre = "-";
			if($B_Qstart-1 > 0){
				$pre = substr($seq, 0, $B_Qstart-1);
			}
			$mid = "-";
			if(($P_Qstart-1)-$B_Qend > 0){
				$mid = substr($seq, $B_Qend, ($P_Qstart-1)-$B_Qend);
			}
			$post = "-";
			$post = substr($seq, $P_Qend);
			$bait =~ tr/a-z/A-Z/;
			$mid =~ tr/a-z/A-Z/;
			$prey =~ tr/a-z/A-Z/;
			push(@F, $pre, $bait, $mid, $prey, $post);
		}else{
			$mid = $F[$f2i{"mid"}];
		}
		
		## alignment
#		score, mid_annotate = Dalign(Bio.Seq.Seq(mid, Bio.Alphabet.IUPAC.unambiguous_dna))
		## 2018-12-23 turn back to use continuous_LCS_DP (LCS=longest continuous substring, not allow gap) instead of Bio.pairwise2.align.globalms
		($score, $mid_annotate) = Dalign($mid);
		if($score < $score_cutoff){
			$mid_annotate = ["-"];
			$score = "-";
		}
		
		push(@F, $score, join(", ", @$mid_annotate));
		
		print join("\t", @F)."\n";
	}
}

0;

