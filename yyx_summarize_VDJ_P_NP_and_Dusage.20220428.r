# Usage: Rscript this.r <input_annotated.tlx> <input_stats.txt> <D_usage.tsv> <V_anno.bed> <D_anno.bed> <output.tsv> 
#	[normalize_total_reads (default:100000) (NULL to skip normalization)]


options(stringsAsFactors=FALSE)

args = commandArgs(TRUE)
if(length(args) != 6 && length(args) != 7){
	stop("Error: there should be 5 or 6 command-line arguments\nUsage: Rscript this.r <input_annotated.tlx> <input_stats.txt> <D_usage.tsv> <V_anno.bed> <D_anno.bed> <output.tsv>\n\t[normalize_total_reads (default:100000) (NULL to skip normalization)]")
}



library(tidyverse)

`%.%` = function(x,y) paste0(x,y)

join = function(sep, vec, ...) paste(collapse=sep, c(vec, ...))
cat0 = function(...) cat(sep="", ...)

echo_str <- function(x, sep=" =\t", collapse=", ") deparse(substitute(x)) %.% sep %.% join(collapse, x)
echo <- function(x, sep=" =\t", collapse=", ") cat0(deparse(substitute(x)), sep, join(collapse, x), "\n");



args_strsplit_list = strsplit(args, ",")
input_filenames = args_strsplit_list[[1]]
stats_filenames = args_strsplit_list[[2]]
Dusage_tsv_filenames = args_strsplit_list[[3]]
VH_coord_bed_filename = args_strsplit_list[[4]]
DH_coord_bed_filename = args_strsplit_list[[5]]
output_filename = args_strsplit_list[[6]]

normalize_total = 100000
if(length(args)==7){
	if(grepl("null", args[7], ignore.case=TRUE)){
		normalize_total = NULL
	}else{
		normalize_total = args[7] %>% as.integer
	}
}

echo(input_filenames)
echo(stats_filenames)
echo(Dusage_tsv_filenames)
echo(VH_coord_bed_filename)
stopifnot(length(VH_coord_bed_filename)==1)
echo(DH_coord_bed_filename)
stopifnot(length(DH_coord_bed_filename)==1)
echo(output_filename)
stopifnot(length(output_filename)==1)
echo(normalize_total)

if(file.exists(output_filename)){
	stop("Error: output_file " %.% output_filename %.% " already exists. Please remove it if you want to rerun.")
}



VH_coord_DF = read_tsv(VH_coord_bed_filename, col_names=FALSE) %>% as.data.frame
colnames(VH_coord_DF) = c("chr", "start", "end", "name")

VH_idx = grepl("^IGHV|^VH", VH_coord_DF$name, ignore.case=TRUE)
VH_coord_DF = VH_coord_DF[VH_idx, ]
stopifnot(nrow(VH_coord_DF) > 0)
VH_coord_DF$x = 1:nrow(VH_coord_DF)
echo(dim(VH_coord_DF))   # 235  4


DH_coord_DF = read_tsv(DH_coord_bed_filename, col_names=FALSE) %>% as.data.frame
colnames(DH_coord_DF) = c("chr", "start", "end", "name")

DH_idx = grepl("^IGHD", DH_coord_DF$name, ignore.case=TRUE)
DH_coord_DF = DH_coord_DF[DH_idx, ]
stopifnot(nrow(DH_coord_DF) > 0)
DH_coord_DF$x = -(nrow(DH_coord_DF):1)
echo(dim(DH_coord_DF))   # 235  4



summarize_one_sample_input = function(input_filenames, stats_filenames, Dusage_tsv_filenames, normalize_total = NULL){
	stopifnot(length(input_filenames)==length(stats_filenames))
	germline = c()
	junction = c()
	for(stats_filename in stats_filenames){
		suppressWarnings({
		now_input = read_tsv(stats_filename, col_types = cols())
		})
		colnames(now_input)[1] = "type"
		germline = c(germline, now_input$reads[now_input$type=="uncut"])
		junction = c(junction, now_input$reads[now_input$type=="result"])
	}
	echo(stats_filenames)
	echo(germline)
	echo(junction)
	total = germline+junction
	echo(total)
	germline = sum(germline)
	junction = sum(junction)
	total = sum(total)
	head_DF = data.frame(name=c("Total", "Germline", "Junction"), count=c(total, germline, junction))
	head_DF$count_P = NA
	head_DF$count_NP = NA
	head_DF$count_D_in_VDJ = NA
	
	all_input = character(0)
	for(input_filename in input_filenames){
		suppressWarnings({
		now_input = read_tsv(input_filename, col_types = cols()) %>% select(Jfeatures, Productive) %>% as.data.frame
		})
		echo(dim(now_input))
		all_input = rbind(all_input, now_input)
	}
	echo(dim(all_input))
	
	Dusage_input = NULL
	for(Dusage_tsv_filename in Dusage_tsv_filenames){
		now_input = read_tsv(Dusage_tsv_filename, col_types = cols(), col_names=FALSE)
		colnames(now_input) = c("D", "count_in_VDJ", "pct_in_VDJ")
		if(is.null(Dusage_input)){
			Dusage_input = now_input %>% select(D, count_in_VDJ)
		}else{
			colnames(now_input) = c("D", "count_in_VDJ_2", "pct_in_VDJ_2")
			Dusage_input = Dusage_input %>% merge(now_input) %>% mutate(count_in_VDJ = count_in_VDJ + count_in_VDJ_2) %>% select(D, count_in_VDJ)
		}
	}
	echo(dim(Dusage_input))
	
	tmp_DF = (all_input %>% filter(Jfeatures %in% DH_coord_DF$name))$Jfeatures %>% table %>% as.data.frame
	if(ncol(tmp_DF) == 2){
		colnames(tmp_DF) = c("name", "count")
		DH_count_DF = DH_coord_DF %>% select(x, name) %>% merge(tmp_DF, all.x=TRUE) %>% arrange(x) %>% select(-x)
		DH_count_DF$count[is.na(DH_count_DF$count)] = 0
		DH_count_DF$count_P = NA
		DH_count_DF$count_NP = NA
	}else{
		DH_count_DF = DH_coord_DF %>% select(name)
		DH_count_DF$count = 0
		DH_count_DF$count_P = NA
		DH_count_DF$count_NP = NA
	}
	DH_count_DF$count_D_in_VDJ = NA
	for(ii in 1:nrow(Dusage_input)){
		now_DH_rowidx = which( grepl(Dusage_input$D[ii] %.% " ", DH_count_DF$name) | grepl(Dusage_input$D[ii] %.% "$", DH_count_DF$name) )
		if(length(now_DH_rowidx)==1){
			DH_count_DF$count_D_in_VDJ[now_DH_rowidx] = Dusage_input$count_in_VDJ[ii]
		}else{
			cat("Warning: find " %.% length(now_DH_rowidx) %.% " times for " %.% Dusage_input$D[ii] %.% " in DH table\n")
		}
	}
	
	tmp_DF = (all_input %>% filter(Jfeatures %in% VH_coord_DF$name))$Jfeatures %>% table %>% as.data.frame
	if(ncol(tmp_DF) == 2){
		colnames(tmp_DF) = c("name", "count")
		VH_count_DF = VH_coord_DF %>% select(x, name) %>% merge(tmp_DF, all.x=TRUE) %>% arrange(x) %>% select(-x)
		VH_count_DF$count[is.na(VH_count_DF$count)] = 0
	}else{
		VH_count_DF = VH_coord_DF %>% select(name)
		VH_count_DF$count = 0
	}
	
	tmp_DF = (all_input %>% filter(Jfeatures %in% VH_coord_DF$name, Productive=="T"))$Jfeatures %>% table %>% as.data.frame
	if(ncol(tmp_DF) == 2){
		colnames(tmp_DF) = c("name", "count")
		VH_P_count_DF = VH_coord_DF %>% select(x, name) %>% merge(tmp_DF, all.x=TRUE) %>% arrange(x) %>% select(-x)
		VH_P_count_DF$count[is.na(VH_P_count_DF$count)] = 0
	}else{
		VH_P_count_DF = VH_coord_DF %>% select(name)
		VH_P_count_DF$count = 0
	}
	stopifnot(all(VH_count_DF$name==VH_P_count_DF$name))
	VH_count_DF$count_P = VH_P_count_DF$count
	
	tmp_DF = (all_input %>% filter(Jfeatures %in% VH_coord_DF$name, Productive=="F"))$Jfeatures %>% table %>% as.data.frame
	if(ncol(tmp_DF) == 2){
		colnames(tmp_DF) = c("name", "count")
		VH_NP_count_DF = VH_coord_DF %>% select(x, name) %>% merge(tmp_DF, all.x=TRUE) %>% arrange(x) %>% select(-x)
		VH_NP_count_DF$count[is.na(VH_NP_count_DF$count)] = 0
	}else{
		VH_NP_count_DF = VH_coord_DF %>% select(name)
		VH_NP_count_DF$count = 0
	}
	stopifnot(all(VH_count_DF$name==VH_NP_count_DF$name))
	VH_count_DF$count_NP = VH_NP_count_DF$count
	VH_count_DF$count_D_in_VDJ = NA
	
	head_DF = rbind(head_DF, data.frame(name="DJ", count=sum(DH_count_DF$count), count_P=sum(DH_count_DF$count_P), count_NP=sum(DH_count_DF$count_NP), count_D_in_VDJ=NA))
	head_DF = rbind(head_DF, data.frame(name="VDJ", count=sum(VH_count_DF$count), count_P=sum(VH_count_DF$count_P), count_NP=sum(VH_count_DF$count_NP), count_D_in_VDJ=sum(DH_count_DF$count_D_in_VDJ)))
	
	echo(dim(head_DF))
	echo(dim(DH_count_DF))
	echo(dim(VH_count_DF))
	
	out_DF = rbind(head_DF, DH_count_DF, VH_count_DF)
	## normalize
	if(!is.null(normalize_total)){
		out_DF$count = out_DF$count / total * normalize_total
		out_DF$count_P = out_DF$count_P / total * normalize_total
		out_DF$count_NP = out_DF$count_NP / total * normalize_total
		out_DF$count_D_in_VDJ = out_DF$count_D_in_VDJ / total * normalize_total
	}
	echo(dim(out_DF))
	out_DF
}

summarize_one_sample_input(input_filenames, stats_filenames, Dusage_tsv_filenames, normalize_total) %>% write_tsv(output_filename)


