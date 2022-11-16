#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2022-03-29)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <input.tlx> <ref.fa> <VDJ.bed> <V.fa> <J.fa> <J.aux>
Input:
	<J.aux>  the optional file in IgBLAST to annotate J frame, 4 columns:
		gene name  (e.g. JH1 , IGHJ1*01)
		first coding frame start position (0-based)  (e.g. 0 , 1 , 2)
		chain type  (e.g. JH , JK , JL ; I will ignore this column)
		CDR3 end (e.g. 18 , 13 , 6 , 7)
Output: STDOUT   append several columns as follows
	pre  bait  mid  prey  post   (segmented sequences on each read in <input.tlx>)
	Bfeature  Pfeature   (annotation of overlapping features in <VDJ.bed>)
	Vpart  midO  Jpart   (extended V(D)J sequence)
	InFrame  Stop  Productive
".$version;

if(@ARGV < 6){
	die $usage;
}

my $input_filename = shift(@ARGV);
my $ref_fa_filename = shift(@ARGV);
my $VDJ_bed_filename = shift(@ARGV);
my $V_fa_filename = shift(@ARGV);
my $J_fa_filename = shift(@ARGV);
my $J_aux_filename = shift(@ARGV);



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

my $ref_hash = read_fasta($ref_fa_filename);
my $V_hash = read_fasta($V_fa_filename);
my $J_hash = read_fasta($J_fa_filename);
my ($tmp_hash, $value);
$tmp_hash = $V_hash;
$V_hash = {};
foreach (keys %$tmp_hash){
	$value = $tmp_hash->{$_};
	s/^lcl[|]//;
	$V_hash->{$_} = $value;
}
$tmp_hash = $J_hash;
$J_hash = {};
foreach (keys %$tmp_hash){
	$value = $tmp_hash->{$_};
	s/^lcl[|]//;
	$J_hash->{$_} = $value;
}
my $VJ_hash = {};
foreach (keys %$V_hash){
	$VJ_hash->{$_} = $V_hash->{$_};
}
foreach (keys %$J_hash){
	$VJ_hash->{$_} = $J_hash->{$_};
}


sub read_aux_file{
	my ($input_aux_filename) = @_;
	my $Jaux_hash = {};
	my @F;
	open(IN, $input_aux_filename) or die "Error: cannot open aux file $input_aux_filename for input\n";
	print STDERR "Now reading aux file $input_aux_filename ...\n";
	while(<IN>){
		next if(/^#/);
		s/[\r\n]+$//;
		@F = split/\s+/;
		$Jaux_hash->{$F[0]} = [@F[1..3]];
	}
	close(IN);
	return $Jaux_hash;
}

my $Jaux_hash = read_aux_file($J_aux_filename);
#foreach (sort keys %$Jaux_hash){
#	print STDERR "[DEBUG] Jaux $_\n";
#}
foreach (sort keys %$J_hash){
	if(!exists($Jaux_hash->{$_})){
		print STDERR "Warning: J $_ does not exist in the aux file\n";
	}
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
#print join("\n", reverse_complement("ACCGT", "AAAA"))."\n";



### Editing distance
### ref: https://www.geeksforgeeks.org/edit-distance-dp-5/
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
sub editDistDP{
	my ($str1, $str2) = @_;
	my @S1 = split(//, $str1);
	my @S2 = split(//, $str2);
	my $m = @S1;
	my $n = @S2;
	my @dp = ();
	my ($i, $j);
	for($i=0; $i<=$m; $i++){
		push(@dp, [(0) x ($n+1)]);
	}
	for($i=0; $i<=$m; $i++){
		for($j=0; $j<=$n; $j++){
			if($i==0){
				$dp[$i]->[$j] = $j;
			}elsif($j==0){
				$dp[$i]->[$j] = $i;
			}elsif($S1[$i-1] eq $S2[$j-1]){
				$dp[$i]->[$j] = $dp[$i-1]->[$j-1];
			}else{
				$dp[$i]->[$j] = 1 + min($dp[$i]->[$j-1], $dp[$i-1]->[$j], $dp[$i-1]->[$j-1]);
			}
		}
	}
	return $dp[$m][$n];
}



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
my ($VDJ_S, $VDJ_H) = read_bed($VDJ_bed_filename);


### guess bed format (start 0/1-based, end-excluded/included, and strand ?)
sub best_match_and_guess_bed_format{
	my ($chr, $start, $end, $name, $chr_seq, $VJ_hash) = @_;
	my @VJ_name_vec = ();
	my ($name2, $name3, $name4);
	my @VJ_names = keys %$VJ_hash;
	$name2 = $name;
	if(exists($VJ_hash->{$name2})){
#		print STDERR "[DEBUG] find exact $name\n";
		push(@VJ_name_vec, $name2);
	}else{
		foreach $name3 (@VJ_names){
			$name4 = $name3;
			$name4 =~ s/[ (].*$//;
			$name4 =~ tr/a-z/A-Z/;
			$name4 =~ s/[*].*$//;
			$name4 =~ s/P$//;
			if($name2 eq $name4){
#				print STDERR "[DEBUG] find $name3 for $name\n";
				push(@VJ_name_vec, $name3);
			}
		}
	}
	if(@VJ_name_vec == 0){
		$name2 =~ s/[ (].*$//;
		$name2 =~ tr/a-z/A-Z/;
		$name2 =~ s/P$//;
		if(exists($VJ_hash->{$name2})){
#			print STDERR "[DEBUG] find $name2 for $name\n";
			push(@VJ_name_vec, $name2);
		}else{
			foreach $name3 (@VJ_names){
				$name4 = $name3;
				$name4 =~ s/[ (].*$//;
				$name4 =~ tr/a-z/A-Z/;
				$name4 =~ s/[*].*$//;
				$name4 =~ s/P$//;
				if($name2 eq $name4){
#					print STDERR "[DEBUG] find $name3 for $name\n";
					push(@VJ_name_vec, $name3);
				}
			}
		}
	}
	if(@VJ_name_vec == 0){
		$name2 =~ s/^IG(.)(.)/$2$1/;
		if(exists($VJ_hash->{$name2})){
#			print STDERR "[DEBUG] find $name2 for $name\n";
			push(@VJ_name_vec, $name2);
		}else{
			foreach $name3 (@VJ_names){
				$name4 = $name3;
				$name4 =~ s/[ (].*$//;
				$name4 =~ tr/a-z/A-Z/;
				$name4 =~ s/[*].*$//;
				$name4 =~ s/P$//;
				if($name2 eq $name4){
#					print STDERR "[DEBUG] find $name3 for $name\n";
					push(@VJ_name_vec, $name3);
				}
			}
		}
	}
	if(@VJ_name_vec == 0){
		$name2 =~ s/^([VJ])([HKL])/IG$2$1/;
		if(exists($VJ_hash->{$name2})){
#			print STDERR "[DEBUG] find $name2 for $name\n";
			push(@VJ_name_vec, $name2);
		}else{
			foreach $name3 (@VJ_names){
				$name4 = $name3;
				$name4 =~ s/[ (].*$//;
				$name4 =~ tr/a-z/A-Z/;
				$name4 =~ s/[*].*$//;
				$name4 =~ s/P$//;
				if($name2 eq $name4){
#					print STDERR "[DEBUG] find $name3 for $name\n";
					push(@VJ_name_vec, $name3);
				}
			}
		}
	}
	if(@VJ_name_vec == 0){
		print STDERR "Warning: cannot find name $name in V.fa nor J.fa, so I will check for all\n";
		push(@VJ_name_vec, sort keys %$VJ_hash);
	}
	my $best_VJ_name = undef;
	my $best_key = undef;
	my $best_ED = undef;
	my ($VJ_name, $VJ_seq, %target_seqs, $target_len, $has_exact_match, $start_base, $end_included, $strand);
	my ($now_best_key, $now_best_ED, $ref_seq, $key, $ED);
	foreach $VJ_name (@VJ_name_vec){
		$VJ_seq = $VJ_hash->{$VJ_name};
		$VJ_seq =~ tr/a-z/A-Z/;
		%target_seqs = (
			"+" => $VJ_seq,
			"-" => reverse_complement($VJ_seq)
		);
		$target_len = length($VJ_seq);
		
		$has_exact_match = 0;
		$now_best_key = undef;
		$now_best_ED = undef;
		foreach $end_included (0..1){
			next if(($end-$start+$end_included)!=$target_len);
			foreach $start_base (0..1){
				$ref_seq = substr($chr_seq, $start-$start_base, $target_len);
				$ref_seq =~ tr/a-z/A-Z/;
				foreach $strand (qw/+ -/){
					$key = $start_base.$end_included.$strand;
					$ED = editDistDP($ref_seq, $target_seqs{$strand});
					if(!defined($now_best_ED) || $ED < $now_best_ED){
						$now_best_ED = $ED;
						$now_best_key = $key;
					}
					if($ED==0){
						$has_exact_match = 1;
					}
					last if($has_exact_match);
				}
				last if($has_exact_match);
			}
			last if($has_exact_match);
		}
		if(defined($now_best_ED) && (!defined($best_ED) || $now_best_ED < $best_ED)){
			$best_ED = $now_best_ED;
			$best_key = $now_best_key;
			$best_VJ_name = $VJ_name;
			last if($has_exact_match);
		}
	}
#	print STDERR "[DEBUG] for $name, best match is $best_VJ_name, with distance $best_ED, format $best_key\n";
	return ($best_key, $best_VJ_name, $best_ED);
}

my ($chr_seq, $best_key, $best_name, $best_ED);
my ($chr, $start, $end, $name, $ref_seq);
my %bed_coord_guess_hash = ();
my $NR = 0;
foreach $name (@$VDJ_S){
	$NR++;
	($chr, $start, $end) = @{$VDJ_H->{$name}};
	$chr_seq = hash_get_default($ref_hash, $chr, "-", "Warning: ($NR/".scalar(@$VDJ_S).") for $name, chr $chr does not exist in ref.fa\n");
	if($chr_seq ne "-"){
		if($name =~ /[VJ]/i){
			($best_key, $best_name, $best_ED) = best_match_and_guess_bed_format($chr, $start, $end, $name, $chr_seq, $VJ_hash);
			if(defined($best_ED)){
				if($best_ED <= 0){
					push(@{$VDJ_H->{$name}}, $best_key, $best_name);
					$bed_coord_guess_hash{$best_key}++;
					print STDERR "[DEBUG] ($NR/".scalar(@$VDJ_S).")  with format $best_key ,  $name  can be exactly matched to  $best_name\n";
				}else{
					print STDERR "Warning: ($NR/".scalar(@$VDJ_S).") cannot match $name in V.fa nor J.fa, best match to $best_name with distance $best_ED\n";
				}
			}else{
				print STDERR "Warning: ($NR/".scalar(@$VDJ_S).") cannot match $name in V.fa nor J.fa, no best match\n";
			}
		}
	}
}
if(scalar(keys %bed_coord_guess_hash)==0){
	die "Error: no one in VDJ.bed can be matched in V.fa nor J.fa, so I cannot guess the coordincate format of VDJ.bed.\n";
}
my ($key);
$best_key = undef;
foreach $key (sort { $bed_coord_guess_hash{$b} <=> $bed_coord_guess_hash{$a} } keys %bed_coord_guess_hash){
	print STDERR join("\t", "[DEBUG] bed_format", $key, $bed_coord_guess_hash{$key})."\n";
	if(!defined($best_key)){
		$best_key = $key;
	}
}
my ($guess_start_base, $guess_end_included, $guess_strand) = split(//, $best_key);
#die;



sub get_ref_seq{
	my ($chr, $start, $end, $strand) = @_;
	my $ref = "-";
	$ref = hash_get_default($ref_hash, $chr, "-", "Warning: $chr does not exist in reference fasta\n");
	if($ref ne "-"){
		$ref = substr($ref, $start-1, $end-($start-1));
	}
	if($strand eq "-"){
		$ref = reverse_complement($ref);
	}
	return $ref;
}

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

sub get_overlapping_features{
	### judge overlapping features (low efficiency)
	my ($Qchr, $Qstart, $Qend, $Qstrand, $Qjunction) = @_;
	my ($chr, $start, $end, $name, $format_key, $start_base, $end_included, $strand);
	my ($overlap_len, $Jdistance, $Jstart, $Jend);
	my @ans = ();
	foreach $name (@$VDJ_S){
#		($chr, $start, $end, $format_key) = @{$VDJ_H->{$name}};
		($chr, $start, $end, $strand) = adjust_start_end(@{$VDJ_H->{$name}});
		next if($chr ne $Qchr);
		$overlap_len = calculate_overlap_len($Qstart, $Qend, $start, $end);
		if($overlap_len >= 1){
			$Jdistance = 0;
			$Jstart = $Qjunction - $start;
			$Jend = $Qjunction - $end;
			if($Jstart * $Jend > 0){
				$Jdistance = min(abs($Jstart), abs($Jend));
			}
			push(@ans, [$Jdistance, $name, $Qstrand eq $strand]);
		}
	}
	return sort { $a->[0] <=> $b->[0] } @ans;
}

sub judge_VDJname_category{
	my ($name) = @_;
	my $name2 = $VDJ_H->{$name}->[4];
	if(defined($name2)){
		if(exists($V_hash->{$name2})){
			return "V";
		}elsif(exists($J_hash->{$name2})){
			return "J";
		}
	}else{
		if($name =~ /^V/i || $name =~ /^IG.V/i){
			return "V";
		}elsif($name =~ /^J/i || $name =~ /^IG.J/i){
			return "J";
		}
	}
	return "other";
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
sub get_V_part_seq{
	my ($name, $Oseq, $Ostart, $Oend) = @_;
	## Oseq should be in V-D-J orientation
	my ($Rshort, $Oshort);
	my $short_bp = 3;
	my ($chr, $start, $end, $strand) = adjust_start_end(@{$VDJ_H->{$name}});
	my $V_seq = get_ref_seq($chr, $start, $end, $strand);
	$V_seq =~ tr/a-z/A-Z/;
	my $V_overhang_O_len = $Ostart - $start;
	if($strand eq "-"){
		$V_overhang_O_len = $end - $Oend;
	}
	if($V_overhang_O_len >= 0){
		$Rshort = str_right(str_left($V_seq, ($V_overhang_O_len + $short_bp)),$short_bp);
		$Oshort = str_left($Oseq, $short_bp);
		if($Rshort ne $Oshort){
			print STDERR "Warning: R=$Rshort ne O=$Oshort ?! in get_V_part_seq() for V=$name Oseq=$Oseq Ostart=$Ostart Oend=$Oend\n";
		}
		return str_left($V_seq, $V_overhang_O_len) . $Oseq;
	}else{
		$Rshort = str_left($V_seq, $short_bp);
		$Oshort = str_right(str_left($Oseq, (-$V_overhang_O_len + $short_bp)),$short_bp);
		if($Rshort ne $Oshort){
			print STDERR "Warning: R=$Rshort ne O=$Oshort ?! in get_V_part_seq() for V=$name Oseq=$Oseq Ostart=$Ostart Oend=$Oend\n";
		}
		return str_trim_left($Oseq, -$V_overhang_O_len);
	}
}
sub get_J_part_seq{
	my ($name, $Oseq, $Ostart, $Oend) = @_;
	## Oseq should be in V-D-J orientation
	my ($Rshort, $Oshort);
	my $short_bp = 3;
	my ($chr, $start, $end, $strand) = adjust_start_end(@{$VDJ_H->{$name}});
	my $J_seq = get_ref_seq($chr, $start, $end, $strand);
	$J_seq =~ tr/a-z/A-Z/;
	my $J_overhang_O_len = $end - $Oend;
	if($strand eq "-"){
		$J_overhang_O_len = $Ostart - $start;
	}
	if($J_overhang_O_len >= 0){
		$Rshort = str_left(str_right($J_seq, ($J_overhang_O_len + $short_bp)),$short_bp);
		$Oshort = str_right($Oseq, $short_bp);
		if($Rshort ne $Oshort){
			print STDERR "Warning: R=$Rshort ne O=$Oshort ?! in get_J_part_seq() for J=$name Oseq=$Oseq Ostart=$Ostart Oend=$Oend\n";
		}
#		print STDERR "[DEBUG] R=".str_left(str_right($J_seq, ($J_overhang_O_len + 3)),3)." O=".str_right($Oseq, 3)."\n";
		return $Oseq . str_right($J_seq, $J_overhang_O_len);
	}else{
		$Rshort = str_right($J_seq, $short_bp);
		$Oshort = str_left(str_right($Oseq, (-$J_overhang_O_len + $short_bp)),$short_bp);
		if($Rshort ne $Oshort){
			print STDERR "Warning: R=$Rshort ne O=$Oshort ?! in get_J_part_seq() for J=$name Oseq=$Oseq Ostart=$Ostart Oend=$Oend\n";
		}
#		print STDERR "[DEBUG] R=".str_right($J_seq, 3)." O=".str_left(str_right($Oseq, (-$J_overhang_O_len + 3)),3)."\n";
		return str_trim_right($Oseq, -$J_overhang_O_len);
	}
}
sub judge_in_frame{
	my ($VpartAddMid_len, $Jstart, $Jend, $Jname) = @_;
	my ($chr, $start, $end, $strand) = adjust_start_end(@{$VDJ_H->{$Jname}});
#	if(@{$VDJ_H->{$Jname}} < 5){
#		print STDERR "Warning: cannot judge in-frame for J=$Jname because it cannot match any J in fa file\n";
#		return "-";
#	}
	my $Jname2 = $VDJ_H->{$Jname}->[4];
	if(!defined($Jname2) || !exists($Jaux_hash->{$Jname2})){
		print STDERR "Warning: cannot judge in-frame for J=$Jname because it is not in aux file\n";
		return "-";
	}
	my $first_codon_start_pos = $Jaux_hash->{$Jname2}->[0];
	if($strand eq "+"){
		if(($VpartAddMid_len + $first_codon_start_pos - ($Jstart-$start)) % 3 == 0){
			return "T";
		}else{
			return "F";
		}
	}else{
		if(($VpartAddMid_len + $first_codon_start_pos - ($end-$Jend)) % 3 == 0){
			return "T";
		}else{
			return "F";
		}
	}
}

my %stop_codon_hash = (
	"TAA" => 1,
	"TAG" => 1,
	"TGA" => 1,
);
sub any_stop_codon{
	my ($seq) = @_;
	my $len = length($seq);
	my ($ii);
	for($ii=0; $ii<$len; $ii+=3){
		if(exists($stop_codon_hash{substr($seq, $ii, 3)})){
			return "T";
		}
	}
	return "F";
}


my @requried_fields = qw/B_Qstart B_Qend Qstart Qend Seq Rname Rstart Rend Strand B_Rname B_Rstart B_Rend B_Strand/;

my @fields;
my %f2i;
my (@F);
my @output_fields;
my ($seq, $B_Qstart, $B_Qend, $P_Qstart, $P_Qend, $prey, $bait, $mid, $pre, $post);
my ($P_Rname, $P_Rstart, $P_Rend, $P_Strand, $B_Rname, $B_Rstart, $B_Rend, $B_Strand, $Rbait, $Rprey);
my ($P_Junction, $B_Junction, @P_overlapping_features, @B_overlapping_features);
my ($P_Jdistance, $P_name, $P_isSameStrand, $B_Jdistance, $B_name, $B_isSameStrand);
my ($V_part_seq, $mid_O_seq, $J_part_seq, $is_in_frame, $has_any_stop_codon, $is_productive);
open(IN, $input_filename) or die "Error: cannot open tlx file $input_filename for input\n";
print STDERR "Now parsing tlx file $input_filename ...\n";
$_ = <IN>;
s/[\r\n]+$//;
@F = split/\t/;
@fields = @F;
for($i=0; $i<@F; $i++){
	$f2i{$F[$i]} = $i;
}
foreach (@requried_fields){
	if(!exists($f2i{$_})){
		die "Error: cannot find colname $_ in tlx file $input_filename\n";
	}
}
@output_fields = @fields;
push(@output_fields, qw/pre bait mid prey post Bfeature Pfeature Vpart midO Jpart InFrame Stop Productive/);
print join("\t", @output_fields)."\n";
$NR = 1;
while(<IN>){
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
	push(@F, $pre, $bait, $mid, $prey, $post);
#	$bait_to_prey = $bait.$mid.$prey;
	$bait =~ tr/a-z/A-Z/;
	$mid =~ tr/a-z/A-Z/;
	$prey =~ tr/a-z/A-Z/;
	
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
	@B_overlapping_features = get_overlapping_features($B_Rname, $B_Rstart, $B_Rend, $B_Strand, $B_Junction);
	
	$P_name = "-";
	if(@P_overlapping_features > 0){
		($P_Jdistance, $P_name, $P_isSameStrand) = @{$P_overlapping_features[0]};
	}
	$B_name = "-";
	if(@B_overlapping_features > 0){
		($B_Jdistance, $B_name, $B_isSameStrand) = @{$B_overlapping_features[0]};
	}
	push(@F, $B_name, $P_name);
	
	### extract V J information, then judge In-frame , Stop codon
	$V_part_seq = "-";
	$J_part_seq = "-";
	$mid_O_seq = "-";
	$is_in_frame = "-";
	$has_any_stop_codon = "-";
	$is_productive = "-";
#	print "[DEBUG] $B_name category = ".judge_VDJname_category($B_name)."\n";
#	print "[DEBUG] $P_name category = ".judge_VDJname_category($P_name)."\n";
	if(judge_VDJname_category($B_name) eq "J" && judge_VDJname_category($P_name) eq "V"){
#		print "[DEBUG] bait on J=$B_name, prey on V=$P_name\n";
		if($P_isSameStrand || $B_isSameStrand){
			print STDERR "Warning: bait on J=$B_name, prey on V=$P_name, read same strand as V-D-J (".($P_isSameStrand?"T":"F").($B_isSameStrand?"T":"F").") ?! for Line $NR\n";
		}else{   # not $P_isSameStrand && not $B_isSameStrand
			$V_part_seq = get_V_part_seq($P_name, $P_isSameStrand ? $prey : reverse_complement($prey), $P_Rstart, $P_Rend);
			$J_part_seq = get_J_part_seq($B_name, $B_isSameStrand ? $bait : reverse_complement($bait), $B_Rstart, $B_Rend);
			$mid_O_seq = reverse_complement($mid);
			$is_in_frame = judge_in_frame(length($V_part_seq)+length($mid), $B_Rstart, $B_Rend, $B_name);
			$has_any_stop_codon = any_stop_codon($V_part_seq.$mid_O_seq.$J_part_seq);
			$is_productive = ($is_in_frame eq "T" && $has_any_stop_codon eq "F") ? "T" : "F";
		}
	}elsif(judge_VDJname_category($P_name) eq "J" && judge_VDJname_category($B_name) eq "V"){
#		print "[DEBUG] bait on V=$B_name, prey on J=$P_name\n";
		if(!$P_isSameStrand || !$B_isSameStrand){
			print STDERR "Warning: bait on V=$B_name, prey on J=$P_name, read not same strand as V-D-J (".($P_isSameStrand?"T":"F").($B_isSameStrand?"T":"F").") ?! for Line $NR\n";
		}else{   # $P_isSameStrand && $B_isSameStrand
			$V_part_seq = get_V_part_seq($B_name, $B_isSameStrand ? $bait : reverse_complement($bait), $B_Rstart, $B_Rend);
			$J_part_seq = get_J_part_seq($P_name, $P_isSameStrand ? $prey : reverse_complement($prey), $P_Rstart, $P_Rend);
			$mid_O_seq = $mid;
			$is_in_frame = judge_in_frame(length($V_part_seq)+length($mid), $P_Rstart, $P_Rend, $P_name);
			$has_any_stop_codon = any_stop_codon($V_part_seq.$mid_O_seq.$J_part_seq);
			$is_productive = ($is_in_frame eq "T" && $has_any_stop_codon eq "F") ? "T" : "F";
		}
	}
	
	push(@F, $V_part_seq, $mid_O_seq, $J_part_seq, $is_in_frame, $has_any_stop_codon, $is_productive);
	print join("\t", @F)."\n";
}
close(IN);

0;


