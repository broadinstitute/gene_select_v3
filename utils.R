get_sign <- function(var) {
    if (var < 0) 
        return ("neg")
    else if (var >= 0)
        return ("nonneg")
    else {
        print("problem is getting sign")
    }
}


get_splitted_parts <- function(lstr, reg_str) {
    part_lst <- strsplit(lstr, reg_str)
    parts <- part_lst[[1]]
    return (parts)
}

combine_parts <- function(sample_ids, strains, treatments, times) {
    combined_res_pre <- apply(expand.grid(strains, treatments, times), 1, paste, collapse="_")
    combined_res_lst <- lapply(combined_res_pre,
        function(pre, sample_ids) {
            matched_ones <- sample_ids[startsWith(sample_ids, pre)]
            return (matched_ones)
         }, sample_ids)


    combined_res <- unlist(combined_res_lst)
    return (combined_res)

}

get_tp_parts <- function(sample_ids, time_pos) {
    vals <- lapply(sample_ids, function(x) {strsplit(x, '_')[[1]][time_pos]})
    tp_parts <- sort(unique(unlist(vals)))
    return (tp_parts)
}

get_samples <- function(sample_ids, pos, val) {

    retvals <- lapply(sample_ids, function(x, lpos, lval){
            lval1 <- strsplit(x, '_')[[1]][lpos]
            if(lval1 == lval)
                return (TRUE)
            else
                return (FALSE)
            }, pos, val)
    retvals1 <- unlist(retvals)
    ret_samples <- sample_ids[retvals1]
    return (ret_samples)
}

get_sample_groups <- function(sample_ids) {

    # A sample: Ebc_S_gent_+_010m_rep1
    # Strain types : S & R
    # treatment types : + & -

    print("sample_ids")
    print(sample_ids)

    straintype_pos <- 2
    treatment_pos <- 4
    time_pos <- 5
    sus_samples <- get_samples(sample_ids, straintype_pos, 'S')
    print("sus_samples: ")
    print(sus_samples)
    res_samples <- get_samples(sample_ids, straintype_pos, 'R')
    print("res_samples: ")
    print(res_samples)
    treated_samples <- get_samples(sample_ids, treatment_pos, '+')
    print("treated_samples: ")
    print(treated_samples)
    untreated_samples <- get_samples(sample_ids, treatment_pos, '-')
    print("untreated_samples")
    print(untreated_samples)
    tp_parts <- get_tp_parts(sample_ids, time_pos)
    print("tp_parts: ")
    print(tp_parts)


    retval <- list()

    tp_len <- length(tp_parts)
     for (j in 1:tp_len) {
        ltime_point <- tp_parts[j]
        ltime_samples <- get_samples(sample_ids, time_pos, ltime_point)

        treated_str <- paste0("sus_treated_time", j)
        lsamples_treated <- list(sus_samples, treated_samples, ltime_samples)
        retval[[treated_str]] <- Reduce(intersect, lsamples_treated)

        untreated_str <- paste0("sus_untreated_time", j)
        lsamples_untreated <- list(sus_samples, untreated_samples, ltime_samples)
        retval[[untreated_str]] <- Reduce(intersect, lsamples_untreated)

        res_treated_str <- paste0("res_treated_time", j)
        lsamples_res_treated <- list(res_samples, treated_samples, ltime_samples)
        retval[[res_treated_str]] <- Reduce(intersect, lsamples_res_treated)


        res_untreated_str <- paste0("res_untreated_time", j)
        lsamples_res_untreated <- list(res_samples, untreated_samples, ltime_samples)
        retval[[res_untreated_str]] <- Reduce(intersect, lsamples_res_untreated)
    }
    
    exp_conds = list("treated_parts" = c("+"), "untreated_parts" = c("-"), "tp_parts" = tp_parts)
    retval$"exp_conds" = exp_conds
    return (retval)
}


get_colData_small <- function(sample_groups, cond2, cond1) {
    lcond2_samples <- sample_groups[[cond2]]
    lcond1_samples <- sample_groups[[cond1]]
    lcond2_rep <- rep(cond2, length(lcond2_samples))
    lcond1_rep <- rep(cond1, length(lcond1_samples))
    condition <- c(lcond2_rep, lcond1_rep)
    all_samples <- c(lcond2_samples, lcond1_samples)

    colData <- data.frame(condition)
    colnames(colData) <- "condition"
    rownames(colData) <- all_samples
    return (colData)
}


get_repeated_vals <- function(repeat_tag, sample_groups) {
    tp_len <- length(sample_groups$"exp_conds"$"tp_parts")
    timepoints <- tp_len
    condition <- c()
    samples <- c()
    for (j in 1:timepoints) {
        rep_str <- paste0(repeat_tag, j)
        samples_arr <- unlist(sample_groups[[rep_str]])
        rep_arr <- rep(rep_str, length(samples_arr))
        condition <- c(condition, rep_arr)
        samples <- c(samples, samples_arr)
    }

    retval <- list("condition" = condition, "samples" = samples)
    return (retval)
}

get_colData <- function(sample_groups, hasRes = TRUE) {

    sus_treated_vals <- get_repeated_vals("sus_treated_time", sample_groups)
    sus_untreated_vals <- get_repeated_vals("sus_untreated_time", sample_groups)
    condition <- c(sus_treated_vals[["condition"]], sus_untreated_vals[["condition"]])
    all_samples <- c(sus_treated_vals[["samples"]], sus_untreated_vals[["samples"]])

    if (hasRes) {
        res_treated_vals <- get_repeated_vals("res_treated_time", sample_groups)
        res_untreated_vals <- get_repeated_vals("res_untreated_time", sample_groups)
        condition <- c(condition, res_treated_vals[["condition"]], res_untreated_vals[["condition"]])
        all_samples <- c(all_samples, res_treated_vals[["samples"]], res_untreated_vals[["samples"]])
    }

    colData <- data.frame(condition)
    colnames(colData) <- "condition"
    rownames(colData) <- all_samples
    return (colData)
}

get_colData_mixed <- function(sample_groups, hasRes = TRUE) {

    sus_untreated_vals <- get_repeated_vals("sus_untreated_time", sample_groups)
    condition <- c(sus_untreated_vals[["condition"]])
    all_samples <- c(sus_untreated_vals[["samples"]])

    if (hasRes) {
        res_untreated_vals <- get_repeated_vals("res_untreated_time", sample_groups)
        condition <- c(condition, res_untreated_vals[["condition"]])
        all_samples <- c(all_samples, res_untreated_vals[["samples"]])
    }

    colData <- data.frame(condition)
    colnames(colData) <- "condition"
    rownames(colData) <- all_samples
    return (colData)
}

get_count_tbl <- function(count_file) {
    count_tbl1 <- read.csv(count_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
    genes <- rownames(count_tbl1)
    cds_genes <- grepl("^CDS", genes)
    count_tbl <- count_tbl1[cds_genes, ]
    return (count_tbl)
}

get_baseMean_lim <- function(lres, base_lim) {
    lres_sorted <- lres[order(-lres$baseMean), ]
    rowlen <- dim(lres)[1]
    rowlen_lim <- floor(rowlen*(100-base_lim)/100)
    
    if (rowlen_lim > rowlen) {
        rowlen_lim = rowlen
    } else if (rowlen_lim < 1) {
        rowlen_lim = 1
    }
    baseMean_lim <- lres_sorted[rowlen_lim, "baseMean"]
    baseMean_100p <- lres_sorted[1, "baseMean"]
    print(paste0("baseMean_100p: ", baseMean_100p))
    return (baseMean_lim)
}


deseq_condwise <- function(dds, cond2, cond1, lfcth, padjth, altH = "greaterAbs", base_lim = 50) {
    lres <- results(dds, altHypothesis = altH, lfcThreshold = lfcth, contrast=c("condition", cond2, cond1))
    lres_top1 <- subset(lres, padj < padjth)
    baseMean_lim_val <- get_baseMean_lim(lres, base_lim)
    lres_top <- subset(lres_top1, baseMean >= baseMean_lim_val)
    retval <- list("lres_top" = lres_top, "lres" = lres, "baseMean_lim_val" = baseMean_lim_val)
    print_log(paste0("baseMean_", base_lim, "p: ", baseMean_lim_val))
    print_log(paste0("cond2: ", cond2))
    print_log(paste0("cond1: ", cond1))
    print_log(paste0("lfcth: ", lfcth))
    print_log(paste0("padjth: ", padjth))
    print_log(paste0("altH: ", altH))
    print_log(paste0("Gene_count: ", dim(lres)[1]))
    print_log(paste0("Gene_count after padj: ", dim(lres_top1)[1]))
    print_log(paste0("Gene_count_top: ", dim(lres_top)[1]))
    print_log("....................................")
    cat("\n")

    return (retval)
}

print_logf <- function(lstr) {
    write(lstr, file = logfile, append = TRUE)
}

print_log <- function(lstr) {
    print(lstr)
    write(lstr, file = logfile, append = TRUE)
}


deseq_condwise_part <- function(countData, sample_groups, cond2, cond1, lfcth, padjth, altH = "greaterAbs", use_beta_prior = FALSE, base_lim = 50) {
    colData_s <- get_colData_small(sample_groups, cond2, cond1)
    lsamples <- rownames(colData_s)
    countData_s <- countData[, lsamples]
    print(colData_s)
    print_log(paste0("Info: beta prior is set to ", use_beta_prior))
    dds <- NULL
    if (use_beta_prior) {
        dds <- DESeqDataSetFromMatrix(countData = countData_s, colData = colData_s, design = ~ condition)
    } else {
        dds <- DESeqDataSetFromMatrix(countData = countData_s, colData = colData_s, design = ~ 0 + condition)
    }
    dds <- DESeq(dds, betaPrior = use_beta_prior)
    retval <- deseq_condwise(dds, cond2, cond1, lfcth, padjth, altH, base_lim)
    return (retval)
}

print_data_ma_between_conds <- function(out_tbl, outdir, stag, exp_conds, timepos, beta_prior_str, data_usage_str) {
    prefix <- paste(stag, exp_conds$"treated_parts"[1], '_VS_', exp_conds$"untreated_parts"[1], exp_conds$"tp_parts"[timepos], beta_prior_str, data_usage_str, sep = "_")
    print_log(paste0("prefix_str: ", prefix))
    outfile <- paste0(outdir, "/", prefix, ".txt")
    write.table(as.data.frame(out_tbl$"lres_top"), outfile, sep = "\t")
    prefix <- paste(stag, exp_conds$"treated_parts"[1], '_VS_', exp_conds$"untreated_parts"[1], exp_conds$"tp_parts"[timepos], beta_prior_str, data_usage_str, "allgenes", sep = "_")
    outfile <- paste0(outdir, "/", prefix, ".txt")
    write.table(as.data.frame(out_tbl$"lres"), outfile, sep = "\t")
    outfile_MA <- paste0(outdir, "/", prefix, "_MA.pdf")
    pdf(outfile_MA)
    plotMA(out_tbl$"lres", main=prefix, ylim=c(-5,5))
    dev.off()

}

print_data_ma_between_tps <- function(out_tbl, outdir, stag, exp_conds, cond, timepos2, timepos1, beta_prior_str, data_usage_str) {
    cond_str <- NULL
    if (cond == "treated") {
        cond_str <- exp_conds$"treated_parts"[1]
    } else if (cond == "untreated") {
        cond_str <- exp_conds$"untreated_parts"[1]
    }
    prefix <- paste(stag, cond_str, exp_conds$"tp_parts"[timepos2], "_VS_", exp_conds$"tp_parts"[timepos1], beta_prior_str, data_usage_str, sep = "_")
    print_log(paste0("prefix_str: ", prefix))
    outfile <- paste0(outdir, "/", prefix, ".txt")
    write.table(as.data.frame(out_tbl$"lres_top"), outfile, sep = "\t")
    prefix <- paste(stag, cond_str, exp_conds$"tp_parts"[timepos2], "_VS_", exp_conds$"tp_parts"[timepos1], beta_prior_str, data_usage_str, "allgenes", sep = "_")
    outfile <- paste0(outdir, "/", prefix, ".txt")
    write.table(as.data.frame(out_tbl$"lres"), outfile, sep = "\t")
    outfile_MA <- paste0(outdir, "/", prefix, "_MA.pdf")
    pdf(outfile_MA)
    plotMA(out_tbl$"lres", main=prefix, ylim=c(-5,5))
    dev.off()
}



get_lst <- function(lst, end_str) {
    lnames <- names(lst)
    reg_str <- paste0(end_str, "$")
    lvals <- grepl(reg_str, lnames)
    ret_lst <- lst[lvals]
    return(ret_lst)
}

select_genes <- function(gene_lst, tp_len, res_lst, pval_lim_raw) {

    count <- 1
    gene_cond_lst <- list()
    gene_cond_raw_lst <- list()
    pval_term <- "padj"
    pval_term_raw <- "pvalue"

    gene_lst_len <- length(gene_lst)
    total_count <- 0
    last_print <- -1
    for (lgene in gene_lst) {
        for (j in 1:tp_len) {
            logFC_sign_j <- NA
            success_j <- FALSE
            treated_term <- paste0('sus_treated_time', j)
            lgene_cond <- paste0(lgene, "__", treated_term)
            pval_arr <- c()
            pval_raw_arr <- c()
            for (k in 1:tp_len) {
                untreated_term <- paste0('sus_untreated_time', k)
                name_term = paste0(treated_term, "__", untreated_term)
                ldata <- res_lst[[name_term]]
                basemean_lim_val <- basemean_lst[[name_term]]
                basemean_k <- ldata[lgene, "baseMean"]
                PValue_k <- ldata[lgene, pval_term]
                PValue_raw_k <- ldata[lgene, pval_term_raw]
                logFC_k <- ldata[lgene, "log2FoldChange"]

                drop_tag <- paste0("lgene: ", lgene, " treated: ", treated_term, " untreated: ", untreated_term)
                if (is.na(logFC_k)) {
                    success_j <- FALSE
                    err_str <- paste0("Dropped: NA val of logFC, ", drop_tag)
                    print_logf(err_str)
                    break
                }

                if (is.na(PValue_k)) {
                    success_j <- FALSE
                    err_str <- paste0("Dropped: NA val of adjusted PValue, ", drop_tag)
                    print_logf(err_str)
                    break
                }
                
                if (is.na(basemean_k)) {
                    success_j <- FALSE
                    err_str <- paste0("Dropped: NA val of basemean, ", drop_tag)
                    print_logf(err_str)
                    break
                }

                if (is.na(PValue_raw_k)) {
                    success_j <- FALSE
                    err_str <- paste0("Dropped: NA val of unadjusted p value, ", drop_tag)
                    print_logf(err_str)
                    break
                }
                 
                if (basemean_k < basemean_lim_val) {
                    success_j <- FALSE
                    err_str <- paste0("Dropped: basemean ", basemean_k, " < basemean_lim_val: ", basemean_lim_val, " , ", drop_tag)
                    print_logf(err_str)
                    break
                }

                if (PValue_raw_k > pval_lim_raw) {
                    err_str <- paste0("Dropped: raw p-value ", PValue_raw_k, " >  raw p-value limit ", pval_lim_raw, " , ", drop_tag)
                    print_logf(err_str)
                    success_j <- FALSE
                    break
                }

                logFC_sign_k <- get_sign(logFC_k)

                if (k == 1) {
                    logFC_sign_j <- logFC_sign_k
                    success_j <- TRUE
                    pval_arr <- c(pval_arr, PValue_k)
                    pval_raw_arr <- c(pval_raw_arr, PValue_raw_k)

                } else {
                    if (logFC_sign_j != logFC_sign_k) {
                        errstr <- paste0("Dropped: change of sign og logFC from ", logFC_sign_j, " to ", logFC_sign_k, " , ", drop_tag)
                        print_logf(errstr)
                        success_j <- FALSE
                        break
                    } else {
                        pval_arr <- c(pval_arr, PValue_k)
                        pval_raw_arr <- c(pval_raw_arr, PValue_raw_k)
                    }
                }
            }
            if (success_j) {
                success_str <- paste0("count: ", count, " gene: ", lgene, " cond: ", treated_term)
                print_logf(success_str)
                count <- count + 1
                lgene_cond <- paste0(lgene, "__", treated_term)
                gene_cond_lst[[lgene_cond]] <- pval_arr
                gene_cond_raw_lst[[lgene_cond]] <- pval_raw_arr
            }
        }
        total_count <- total_count + 1
        total_comp <- (total_count / gene_lst_len) *100
        total_comp_floor <- floor(total_comp)
        if (total_comp_floor %% 5 == 0) {
            if (total_comp_floor != last_print) {
                out_str <- paste0("Gene selection complete: ", total_comp_floor, "%")
                print_log(out_str)
                last_print <- total_comp_floor
            }
        }
    }

    ret_lst <- list("gene_cond_lst" = gene_cond_lst, "gene_cond_raw_lst" = gene_cond_raw_lst)
    return (ret_lst)
}

get_meta_pval <- function(lpvals, method) { 
   val <- NA 
   if (method == "fisher") { 
       val <- sumlog(lpvals) 
   } else if (method == "stouffer") { 
       val <- sumz(lpvals) 
   } 
   return (val$p)
}

