#!/usr/bin/env Rscript

source("/broad/IDP-Dx_work/nirmalya/research/gene_select_v3/utils.R")
suppressMessages(library(DESeq2))
suppressMessages(library(metap))
suppressMessages(library(docopt))


'Gene selection script using two way factor analysis

Usage: GeneSelect.R --stag <stag> -c <counts> -o <outdir> [--l_cand <l2fc_cand> --base_lim <baseMean % limit> --no_res] 

options:
  -c <count> --count <count>
  -o <outdir> --outdir <outdir>
  --stag <stag>
  --l_cand <l2fc_cand> [default: 1]
  --base_lim <baseMean % limit> [default: 50]' -> doc
# what are the options? Note that stripped versions of the parameters are added to the returned list

opts <- docopt(doc)
str(opts)  

count_file <- opts$count
outdir <- opts$outdir
no_res <- opts$no_res


# p_cand is hard coded since we are not using the adjusted p_value cutoff
p_cand <- 0.05

l_cand <- as.numeric(opts$l_cand)

stag <- opts$stag

base_lim_val <- as.numeric(opts$base_lim)

print(paste0("count_file: ", count_file))
print(paste0("l_cand: ", l_cand))
print(paste("base_lim: ", base_lim_val))
print(paste("stag: ", stag))
print(paste0("no_res", no_res))

dir.create(outdir, recursive = TRUE)
ma_dir <- paste0(outdir, "/MA_plots")
dir.create(ma_dir)
tbls_dir <- paste0(outdir, "/DESeq_tbls")
dir.create(tbls_dir)


logfile <- paste0(outdir, "/", stag, "_logfile.txt")
print(logfile)
file.create(logfile)

pval_logfile <- paste0(outdir, "/", stag, "_pval_logfile.txt")
print(pval_logfile)
file.create(pval_logfile)


print(paste0("count_file: " , count_file))
count_tbl <- get_count_tbl(count_file)
print(paste0("count_tbl: ", dim(count_tbl)[2]))

sample_ids <- colnames(count_tbl)
sample_groups <- get_sample_groups(sample_ids)
print("sample_groups")
print(sample_groups)
colData <- get_colData(sample_groups)
str(colData)
lsamples <- rownames(colData)
countData <- count_tbl[, lsamples] 
print("countData")
str(countData)
tp_parts <- sample_groups$exp_conds$tp_parts

print("Running the DESeq2 tests.")
res_lst <- list()
basemean_lst <- list()
gene_count <-  dim(countData)[1]
tp_len <- length(tp_parts)

for (k in 1:tp_len) {
    treated_term <- paste0('sus_treated_time', k)
    untreated_term <- paste0('sus_untreated_time', k) 
    print(treated_term)
    print(untreated_term)
    contrast_val <- get_contrast_str(sample_groups, treated_term, untreated_term)
    print(contrast_val)
    print("---------")
    lres <- deseq_condwise_part(countData, sample_groups, treated_term, 
        untreated_term, lfcth = l_cand, padjth = p_cand, 
        altH = "greaterAbs", use_beta_prior = FALSE, base_lim = base_lim_val) 
    name_term = paste0(treated_term, "__", untreated_term)
    res_lst[[name_term]] <- data.frame(lres$lres)
    basemean_lst[[name_term]] <- lres$baseMean_lim_val

    print_map_tbls(lres, stag, contrast_val, outdir)

}    

has_res <- !no_res
if (has_res) {
    for (k in 1:tp_len) {
        treated_term <- paste0('res_treated_time', k)
        untreated_term <- paste0('res_untreated_time', k) 
        print(treated_term)
        print(untreated_term)
        contrast_val <- get_contrast_str(sample_groups, treated_term, untreated_term)
        print(contrast_val)
        print("---------")
        lres <- deseq_condwise_part(countData, sample_groups, treated_term, 
            untreated_term, lfcth = l_cand, padjth = p_cand, 
            altH = "greaterAbs", use_beta_prior = FALSE, base_lim = base_lim_val) 
        name_term = paste0(treated_term, "__", untreated_term)
        res_lst[[name_term]] <- data.frame(lres$lres)
        basemean_lst[[name_term]] <- lres$baseMean_lim_val

        print_map_tbls(lres, stag, contrast_val, outdir)

    }
}

        

