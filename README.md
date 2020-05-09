# VDJ_annotation
 scripts and pipeline for VDJ annotation of HTGTS tlx files, including iteratively estimate D usage in VDJ joins

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School


## Pipeline

### Step 1. annotate HTGTS tlx files

```
python3 yyx_annotate_HTGTS_VDJ_pipeline.20200219.py  <VDJ_annotation.bed>  <D.fa>  <scripts_dir>  <input.tlx>  <output_prefix>
```

### Step 2. iteratively calculate D usage in VDJ joins

```
perl yyx_reannotate_VDJ_mid_D_usage_iteration_pipeline.20200319.pl
	<scripts_dir> <HTGTS_VDJ_annotated.tsv> <output_prefix> <mm9|mm9AJ>
	[iter_num (default:9)] [allowed_possible_D_num (default:99)]
	[init_ref_D_usage.tsv (default:even_probabilities]
```


## Scripts

