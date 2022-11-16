#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2022-04-29)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: cat <input.tlx> | perl $0 <anno.bed>
Output: STDOUT   append several columns as follows
	pre  bait  mid  prey  post
	junction_overlap_bps  junction_overlap_features
	prey_overlap_bps  prey_overlap_features
	bait_overlap_bps  bait_overlap_features
".$version;

if(@ARGV < 1){
	die $usage;
}

#my $input_filename = shift(@ARGV);
my $anno_bed_filename = shift(@ARGV);


sub read_bed{
	my ($input_bed_filename) = @_;
	my %H = ();
	my @S = ();
	my $NR = 0;
	my (@F, $chr, $start, $end, $name);
	open(IN, $input_bed_filename) or die "Error: cannot open bed file $input_bed_filename for input\n";
	print STDERR "Now reading bed file $input_bed_filename ...\n";
	while(<IN>){
		$NR++;
		s/[\r\n]+$//;
		@F = split/\t/;
#		push(@S, [@F]);
#		if($NR != (@S - 1)){
#			die "Error: strange unalignment between NR=$NR and \@S in read_bed\n";
#		}
		if(@F < 4){
			print STDERR "Warning: Line $NR in bed file ($_) does not have >= 4 columns, so I skip it\n";
			next;
		}
		($chr, $start, $end, $name) = @F[0..3];
#		$name =~ s/ [(][^()]*[)]$//g;
		if(!exists($H{$name})){
			$H{$name} = [$chr, $start, $end];
			push(@S, $name);
		}else{
			print STDERR "Warning: duplicated records for $name in bed file\n";
		}
#		$NR++;
	}
	close(IN);
	return (\@S, \%H);
}
my ($anno_S, $anno_H) = read_bed($anno_bed_filename);



sub calculate_overlap_len{
	my ($s1, $e1, $s2, $e2) = @_;
	## assume sorted $s1 <= $e2, $s2 <= $e2
	if($s1 > $e2 or $s2 > $e1){
		if($s1==$e2+1 or $s2==$e1+1){
			return 0.1;
		}
		return 0;
	}else{
		if($s1 <= $s2 and $s2 <= $e1 and $e1 <= $e2){
			return ($e1-$s2+1);
		}elsif($s2 <= $s1 and $s1 <= $e2 and $e2 <= $e1){
			return ($e2-$s1+1);
		}elsif($s1 <= $s2 and $e2 <= $e1){
			return ($e2-$s2+1);
		}elsif($s2 <= $s1 and $e1 <= $e2){
			return ($e1-$s1+1);
		}else{
			return -1;
		}
	}
}

my ($guess_start_base, $guess_end_included, $guess_strand) = (0, 0, "+");
sub adjust_start_end{
	my ($chr, $start, $end, $format_key) = @_;
	my ($start_base, $end_included, $strand);
	my $is_format_guessed = 0;
	if(defined($format_key)){
		($start_base, $end_included, $strand) = split(//, $format_key);
	}else{
		($start_base, $end_included, $strand) = ($guess_start_base, $guess_end_included, $guess_strand);
		$is_format_guessed = 1;
	}
	## convert to 1-based, end-included
	$start += 1 - $start_base;
	$end += 1 - $start_base;
	$end -= 1 - $end_included;
	return ($chr, $start, $end, $strand, $is_format_guessed);
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

sub get_overlapping_features{
	### judge overlapping features (low efficiency)
	my ($Qchr, $Qstart, $Qend, $Qstrand, $Qjunction) = @_;
	my ($chr, $start, $end, $name, $format_key, $start_base, $end_included, $strand);
	my ($overlap_len, $Jdistance, $Jstart, $Jend);
	my @ans = ();
	foreach $name (@$anno_S){
#		($chr, $start, $end, $format_key) = @{$anno_H->{$name}};
		($chr, $start, $end, $strand) = adjust_start_end(@{$anno_H->{$name}});
		next if($chr ne $Qchr);
		$overlap_len = calculate_overlap_len($Qstart, $Qend, $start, $end);
		if($overlap_len >= 1){
			$Jdistance = 0;
			$Jstart = $Qjunction - $start;
			$Jend = $Qjunction - $end;
			if($Jstart * $Jend > 0){
				$Jdistance = min(abs($Jstart), abs($Jend));
			}
			push(@ans, [$Jdistance, $overlap_len, $name, $Qstrand eq $strand]);
		}
	}
	return sort { $a->[0] <=> $b->[0] } @ans;
}


sub str_left{
	my ($str, $len) = @_;
	return substr($str, 0, $len);
}
sub str_right{
	my ($str, $len) = @_;
	return substr($str, length($str)-$len, $len);
}
sub str_trim_left{
	my ($str, $len) = @_;
	return substr($str, $len);
}
sub str_trim_right{
	my ($str, $len) = @_;
	return substr($str, 0, length($str)-$len);
}


my @requried_fields = qw/B_Qstart B_Qend Qstart Qend Seq Rname Rstart Rend Strand B_Rname B_Rstart B_Rend B_Strand/;

my @fields;
my %f2i;
my (@F, $i, $NR);
my @output_fields;
my ($seq, $B_Qstart, $B_Qend, $P_Qstart, $P_Qend, $prey, $bait, $mid, $pre, $post);
my ($P_Rname, $P_Rstart, $P_Rend, $P_Strand, $B_Rname, $B_Rstart, $B_Rend, $B_Strand, $Rbait, $Rprey);
my ($P_Junction, $B_Junction, @P_overlapping_features, @B_overlapping_features, @J_overlapping_features);
my ($B_name, $P_name, $J_name, $B_overlap_len, $P_overlap_len, $J_overlap_len);
#open(IN, $input_filename) or die "Error: cannot open tlx file $input_filename for input\n";
print STDERR "Now parsing tlx file from STDIN ...\n";
$_ = <STDIN>;
s/[\r\n]+$//;
@F = split/\t/;
@fields = @F;
for($i=0; $i<@F; $i++){
	$f2i{$F[$i]} = $i;
}
foreach (@requried_fields){
	if(!exists($f2i{$_})){
		die "Error: cannot find colname $_ in STDIN\n";
	}
}
@output_fields = @fields;
push(@output_fields, qw/pre bait mid prey post 	junction_overlap_bps junction_overlap_features  prey_overlap_bps prey_overlap_features  bait_overlap_bps bait_overlap_features/);
print join("\t", @output_fields)."\n";
$NR = 1;
while(<STDIN>){
	$NR++;
	s/[\r\n]+$//;
	@F = split/\t/;
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
#	push(@F, $pre, $bait, $mid, $prey, $post);
#	$bait_to_prey = $bait.$mid.$prey;
	$bait =~ tr/a-z/A-Z/;
	$mid =~ tr/a-z/A-Z/;
	$prey =~ tr/a-z/A-Z/;
	push(@F, $pre, $bait, $mid, $prey, $post);

	$P_Rname  = $F[$f2i{"Rname"}];
	$P_Rstart = $F[$f2i{"Rstart"}];
	$P_Rend   = $F[$f2i{"Rend"}];
	$P_Strand = $F[$f2i{"Strand"}];
	$P_Strand = ($P_Strand > 0)? "+" : "-";
	$P_Junction = ($P_Strand eq "+")? $P_Rstart : $P_Rend;
	$B_Rname  = $F[$f2i{"B_Rname"}];
	$B_Rstart = $F[$f2i{"B_Rstart"}];
	$B_Rend   = $F[$f2i{"B_Rend"}];
	$B_Strand = $F[$f2i{"B_Strand"}];
	$B_Strand = ($B_Strand > 0)? "+" : "-";
	$B_Junction = ($B_Strand eq "+")? $B_Rend : $B_Rstart;
	
	### judge overlapping features (low efficiency)
	@P_overlapping_features = get_overlapping_features($P_Rname, $P_Rstart, $P_Rend, $P_Strand, $P_Junction);
	@J_overlapping_features = get_overlapping_features($P_Rname, $P_Junction, $P_Junction, $P_Strand, $P_Junction);
	@B_overlapping_features = get_overlapping_features($B_Rname, $B_Rstart, $B_Rend, $B_Strand, $B_Junction);
	
	$P_overlap_len = "-";
	$P_name = "-";
	if(@P_overlapping_features > 0){
#		$P_name = join(",", map { $_->[1].":".$_->[0].":".$_->[2]; } @P_overlapping_features);
		$P_overlap_len = join(",", map { $_->[1]; } @P_overlapping_features);
		$P_name = join(",", map { $_->[2]; } @P_overlapping_features);
	}
	$B_overlap_len = "-";
	$B_name = "-";
	if(@B_overlapping_features > 0){
#		$B_name = join(",", map { $_->[1].":".$_->[0].":".$_->[2]; } @B_overlapping_features);
		$B_overlap_len = join(",", map { $_->[1]; } @P_overlapping_features);
		$B_name = join(",", map { $_->[2]; } @B_overlapping_features);
	}
	$J_overlap_len = "-";
	$J_name = "-";
	if(@J_overlapping_features > 0){
#		$J_name = join(",", map { $_->[1].":".$_->[0].":".$_->[2]; } @J_overlapping_features);
		$J_overlap_len = join(",", map { $_->[1]; } @J_overlapping_features);
		$J_name = join(",", map { $_->[2]; } @J_overlapping_features);
	}
	push(@F, $J_overlap_len, $J_name, $P_overlap_len, $P_name, $B_overlap_len, $B_name);
	
	print join("\t", @F)."\n";
}
#close(IN);

0;


