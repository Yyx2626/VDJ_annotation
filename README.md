# VDJ_annotation
 scripts and pipeline for VDJ annotation of HTGTS tlx files, including iteratively estimate D usage in VDJ joins

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School


## Setup

The scripts are command-line perl and python scripts, with little dependence on other packages.
They should be able to run in any modern Linux or Mac platform;
or you might need to first install perl and python 3 for your platform.

To use the scripts, just download the scripts and put/unzipped in a folder.
You can rename the folder with any name you like, such as `VDJ_annotation_scripts`.

When executing the following two pipeline scripts, just replace the command-line argument `<scripts_dir>` with the path of the folder.


## Pipeline

### Step 1. annotate HTGTS tlx files

```
python3 yyx_annotate_HTGTS_VDJ_pipeline.20200219.py
	<VDJ_annotation.bed> <D.fa> <scripts_dir> <input.tlx> <output_prefix>
```

Input:

- `<input.tlx>`  generated by HTGTS pipeline (<https://robinmeyers.github.io/transloc_pipeline/>)

- `<scripts_dir>`  should contain the script files called by this pipeline

  - `Yyx_check_col_num.pl`
  - `yyx_tlx2bed.20200126.py`
  - `yyx_sequence_segment_tlx.20181221.py`
  - `yyx_annotate_tlx_with_intersectBedResults.20181223.py`
  - `yyx_annotate_tlx_midD_LCS.20190116.py`
  - `yyx_uniq_count_and_merge.20190111.py`
  - `yyx_show_or_skip_or_retrieve_columns.20190122.py`

- `<VDJ_annotation.bed>`  contains the genomic coordinate range for each V/D/J segment

  - will used by `yyx_annotate_tlx_with_intersectBedResults.20181223.py`

- `<D.fa>`  contains the sequences of D segments, for

  - will used by `yyx_annotate_tlx_midD_LCS.20190116.py`
  
Output:  (where `yyyyMMdd` is the executing date)

- `<output_prefix>.intersectBed_annotated.yyyyMMdd.tlx`

  will append columns of
  - pre, bait, mid, prey, post sequences
  - \*_overlap_bps and \*_overlap_features  for junction, prey and bait

- `<output_prefix>.annotate_tlx_midD_LCS.yyyyMMdd.tlx`

  will append columns of  pre, bait, mid, prey, post, mid_D_score, mid_D_annotate
  
  for mid_D_annotate column, '(rC)' means reverse complement (orientation of D in the read)

- `<output_prefix>.HTGTS_VDJ_annotated.yyyyMMdd.tsv`

  final output tsv (tab-separated values) format file
  
  compared to `<input.tlx>`, append 5+6+2=13 columns
  - pre, bait, mid, prey, post (sequences) : 5 columns
  - \*_overlap_bps and \*_overlap_features  for junction, prey and bait : 3\*2=6 columns
  - mid_D_score and mid_D_annotate : 2 columns

`yyx_annotate_HTGTS_VDJ_pipeline.20200219.py` will sequentially call:

- `Yyx_check_col_num.pl`  to check the number of columns of <input.tlx>

- `yyx_tlx2bed.20200126.py`  to convert tlx to bed format for junction, prey, and bait

  - Output: `<output_prefix>.*.yyyyMMdd.tmp.bed` and `<output_prefix>.*.yyyyMMdd.bed`  
    (these intermediate files will be automatically removed)

- `yyx_sequence_segment_tlx.20181221.py`  to retrieve the segmented sequence (pre, bait, mid, prey, post) for each read

  - Output: `<output_prefix>.sequence_segmented.yyyyMMdd.tlx`  
    (this intermediate file will be automatically removed)

- `yyx_annotate_tlx_with_intersectBedResults.20181223.py`  to annotate tlx junction/prey/bait overlapping with <VDJ_annotation.bed>

  - Output: `<output_prefix>.intersectBed_annotated.yyyyMMdd.tlx`

- `yyx_annotate_tlx_midD_LCS.20190116.py`  to annotate mid by align to <D.fa> by longest continuous substring (LCS) algorithm

  - Output: `<output_prefix>.annotate_tlx_midD_LCS.yyyyMMdd.tlx`

- `yyx_uniq_count_and_merge.20190111.py`  to merge the annotation results above

  - Output: `<output_prefix>.HTGTS_annotate_merged.yyyyMMdd.tsv`  
    (this intermediate file will be automatically removed)

- `yyx_show_or_skip_or_retrieve_columns.20190122.py`  to skip some useless columns and output the final output tsv file

  - Output: `<output_prefix>.HTGTS_VDJ_annotated.yyyyMMdd.tsv`


### Step 2. iteratively estimate D usage in VDJ joins

```
perl yyx_reannotate_VDJ_mid_D_usage_iteration_pipeline.20200319.pl
	<scripts_dir> <HTGTS_VDJ_annotated.tsv> <output_prefix> <mm9|mm9AJ>
	[iter_num (default:9)] [allowed_possible_D_num (default:99)]
	[init_ref_D_usage.tsv (default:even_probabilities]
```

Input: 

- `<HTGTS_VDJ_annotated.tsv>`  generated by yyx_annotate_HTGTS_VDJ_pipeline.20200219.py
- `[init_ref_D_usage.tsv]`  two columns: D, usage(%)

Output:

- `<output_prefix>.allow_*.iter_*.D_reannotated.tsv`
- `<output_prefix>.allow_*.iter_*.D_usage.tsv`

`yyx_reannotate_VDJ_mid_D_usage_iteration_pipeline.20200319.pl` will iteratively call 

- `yyx_reannotate_calculate_mid_D_usage.20200319.pl`  
  which will focus on VDJ joins (Column junction_overlap_features starts with IGHV), 
  assign ambiguous D segments according to  `[init_ref_D_usage.tsv]`,
  then sum up D assignment to estimate D usage in VDJ joins


Later, you may also use some bash loops like below, 
to rearrange the final `<output_prefix>.allow_*.iter_*.D_usage.tsv` files of several samples 
into one table  
(Suppose `[iter_num]` is 9 (by default), `[allowed_possible_D_num]` is 99 (by default), and `<output_prefix>` are step2/sample1, step2/sample2, ...)

```
for allow_num in 99; do
    echo -ne "D\titer_0" > tmp.head
    ls step2/*.allow_${allow_num}.iter_9.D_usage.tsv | head -n1 | while read f; do
        cat $f | perl -ne '@F=split/\t/; if($F[0] eq "-"){ print $F[0]."\t0\n"; }else{ print $F[0]."\t1\n"; }' >tmp.body
    done
    ls step2/*.allow_${allow_num}.iter_9.D_usage.tsv | while read f; do
        g=${f##*/}
        g=${g%%.*}
        echo $g
        echo -ne "\t$g.iter_9.sum\t$g.iter_9.prob"  >>tmp.head
        mv tmp.body tmp
        paste tmp <(cut -f2-3 $f)  >tmp.body
    done
    echo >>tmp.head
    cat tmp.head tmp.body >step2/allow_${allow_num}.iter_9.merged_D_usage_summary.tsv
done
```

## Scripts

### General utilities

- `Yyx_check_col_num.pl`

```
Usage: perl Yyx_check_col_num.pl <files> ...
Options:
	-d STR	set delimiter (default: \t)
	-n INT	the number of head lines that will be checked
		(default: 4; -1 for all lines)

Version: 0.1.0 (2012-12-02)
Author: Adam Yongxin Ye @ CBI
```

- `yyx_show_or_skip_or_retrieve_columns.20190122.py`
```
Usage: cat <input> | python3 yyx_show_or_skip_or_retrieve_columns.20190122.py <show|skip|retrieve> [column_pattern_1] [column_pattern_2] ...
```

### Process tlx files

- `yyx_tlx2bed.20200126.py`
```
Usage: cat <input.tlx> | python this.py <which_part>
Options:
   <which_part> can be: junction (default) | bait | prey
Output: STDOUT   bed format 6 columns
```

- `yyx_sequence_segment_tlx.20181221.py`
```
Usage: cat <input.tlx> | python3 yyx_sequence_segment_tlx.20181221.py
Input: STDIN   <input.tlx>
    should have at least these columns:
      B_Qstart, B_Qend, Qstart, Qend, Seq
Output: STDOUT   append 5 columns
    pre  bait  mid  prey  post
```

- `yyx_annotate_tlx_with_intersectBedResults.20181223.py`
```
Usage: python3 yyx_annotate_tlx_with_intersectBedResults.20181223.py <input.tlx> <anno1.bed> [anno2.bed] ...
```

- `yyx_annotate_tlx_midD_LCS.20190116.py`
```
Usage: cat <input.tlx> | python3 yyx_annotate_tlx_midD_LCS.20190116.py <D.fa> [score_cutoff (default:5)]
Input: STDIN   <input.tlx>
    should have at least these columns:
      B_Qstart, B_Qend, Qstart, Qend, Seq
Output: STDOUT   append 5+2 columns
    pre  bait  mid  prey  post   mid_D_score  mid_D_annotate
```

- `yyx_uniq_count_and_merge.20190111.py`
```
Usage: python this.py <empty_fill> <should_split_output> <output_prefix> <file1> <file2> [file3 ...]
  <file1> can be 'filename' or 'filename:key_col_idx' or 'filename:col1-col2,col3,...' or 'filename:col:fieldname_prefix'
  (should have headline, columns separated by '\t', key_col_idx is 1-based)
Output:
  <output_prefix>.log  =  STDERR
  <output_prefix>.tsv  (separated by '\t')
     shared_field(s) , count_in_file_1, append_fields_in_file_1, count_in_file_2, append_fields_in_file_2, ...
  <output_prefix>.1=???.tsv, <output_prefix>.1x2.tsv, <output_prefix>.1x2x3.tsv ...  if <should_split_output> = True
```


- `yyx_reannotate_calculate_mid_D_usage.20200319.pl`
```
Usage: cat <HTGTS_VDJ_annotated.tsv> | perl yyx_reannotate_calculate_mid_D_usage.20200319.pl <output_prefix> <mm9|mm9AJ>
	[allowed_possible_D_num (default: 0)] [ref_D_usage.tsv] [bait_IGHJ]
Input:
	[ref_D_usage.tsv]	two columns: D, usage(%)
Output:
	<output_prefix>.D_reannotated.tsv
	<output_prefix>.D_usage.tsv

Version: 0.1.2 (2020-03-19)
Author: Adam Yongxin Ye @ BCH
```

