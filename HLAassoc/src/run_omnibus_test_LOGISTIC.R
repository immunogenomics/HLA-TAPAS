library(argparse)
library(stringr)
library(purrr)
library(dplyr)
library(multidplyr)
library(tidyr)
library(data.table)
library(parallel)
library(rcompanion)

options(datatable.fread.datatable = FALSE)
options(stringsAsFactors = FALSE)

read_data <- function(phased_fname, fam_fname, pheno_fname, covs_fname = NULL) {
  x <- fread(phased_fname, header = FALSE)
  phased <- x[6:nrow(x), 2:ncol(x)]
  colnames(phased) <- c("SNP", x[2, 3:ncol(x)])

  fam <- read.table(fam_fname, as.is = T, na.strings = "-9")
  colnames(fam) <- c("FID", "IID", "FAT", "MOT", "SEX", "PHENO")
  pheno <- read.table(pheno_fname, header = T, stringsAsFactors=F, na.strings = c("-9", "NA"))
  colnames(pheno) <- c("FID", "IID", "PHENO")
  pheno <- subset(pheno, IID %in% fam$IID)

  #sex <- read.table(sex_fname, as.is = T, na.strings = c("0", "-9"))
  #colnames(sex) <- c("FID", "IID", "SEX")
  #sex <- subset(sex, IID %in% fam$IID)

  if (!is.null(covs_fname) & file.exists(covs_fname)) {
    covs <- read.table(covs_fname, header = TRUE, stringsAsFactors = T, na.strings = c("-9","NA"))
    covs <- subset(covs, IID %in% fam$IID)
  } else {
    covs <- NULL
  }

  # (2020.05.12. by Wanson Choi)
  # pop = read.table(pop, as.is = T, na.strings = c("0", "-9"))
  # colnames(pop) = c("FID", "IID", "POP")
  # pop = subset(pop, IID %in% fam$IID)
  # print(head(pop))

  # bim <- read.table(bim_fname, as.is = T)

  # update pheno
  #fam$SEX[match(sex$IID, fam$IID)] <- as.factor(sex$SEX)
  fam$SEX <- fam[,5]
  fam$PHENO[match(pheno$IID, fam$IID)] <-  as.numeric(pheno$PHENO)
  fam$PHENO <- as.numeric(fam$PHENO)
  # fam$POP <- pop
  # fam$POP <- pop$POP

  # remove NA samples
  if (any(is.na(fam$PHENO) | is.na(fam$SEX))) {
    phased <- phased[, -(2 * rep(which(is.na(fam$PHENO) | is.na(fam$SEX)), each = 2) + (0:1))]
    fam <- na.omit(fam)
  } # added small exeption handling by Wanson Choi (2020.05.14.)
  n_sample <- nrow(fam)

  return(list(
    phased = phased,
    fam = fam,
    n_sample = n_sample,
    #bim = bim,
    covs = covs
  ))
}


################################################################################
# Functions -- general subfunctions
################################################################################
get_haplo <- function(phased) {
  n_sample <- nrow(phased) / 2
  uphased <- unique(data.table(phased))
  uniq_haplos <- apply(uphased, 1, function(x) {
    str_c(x, collapse = "")
  })
  haplo <- matrix(0, nrow = n_sample, ncol = length(uniq_haplos))
  for (i in 1:n_sample) {
    j <- match(str_c(phased[2 * i - 1, ], collapse = ""), uniq_haplos)
    k <- match(str_c(phased[2 * i, ], collapse = ""), uniq_haplos)
    haplo[i, j] <- haplo[i, j] + 1
    haplo[i, k] <- haplo[i, k] + 1
  }
  colnames(haplo) <- uniq_haplos
  return(haplo)
}

logit <- function(...) {
  glm(..., family = "binomial")
}

glm_obj_base_haplo <- function(pheno, covars = NULL, covars_aa = NULL, regfun = logit) {
  if (is.null(covars)) {
    obj_base <- regfun(pheno ~ 1)
  } else if (is.null(covars_aa)) {
    obj_base <- regfun(pheno ~ covars)
  } else {
    haplo <- get_haplo(covars_aa)
    obj_base <- regfun(pheno ~ haplo + covars)
  }
  return(obj_base)
}


next_covars_haplo <- function(aa_phased, AA_ID, n_sample, prev_covar_aa = NULL, subset = NULL) {
  if (is.null(subset)) {
    phased <- t(aa_phased[aa_phased$AA_ID %in% AA_ID, 2:(n_sample * 2 + 1)])
  } else {
    phased <- t(aa_phased[aa_phased$AA_ID %in% AA_ID, 2 * subset + 1])
  }
  colnames(phased) <- aa_phased[aa_phased$AA_ID %in% AA_ID, ]$SNP
  return(cbind(prev_covar_aa, phased))
}

haplo_to_aa <- function(haplo, SNP, output_composites=FALSE) {
  # SNP = str_split_fixed(SNP, "_", 5)
  # SNP = apply(SNP, 1, function(x){str_c(x[2], x[3], x[5])})
  aa <- do.call(c, lapply(str_locate_all(haplo, "T"), function(x) {
    alleles = SNP[x[, 1]]
    if (!output_composites) {
      alleles = alleles[!str_detect(alleles, "^AA_.*_[A-Zx]{2,}$")]
    }
    if (length(alleles) > 0) {
      str_c(alleles, collapse = ",")
    } else {
      NA_character_
    }
  }))
  return(aa)
}

get_covs = function(fam, covs) {
  as.matrix(do.call(cbind, lapply(covs, function(df) {
    ret = data.frame(
      FID = fam$FID,
      IID = fam$IID,
      matrix(0, nrow = nrow(fam), ncol = ncol(df) - 2)
    )
    df = subset(df, IID %in% fam$IID)
    ret[match(df$IID, ret$IID), 3:ncol(ret)] = df[,3:ncol(df)]
    colnames(ret)[3:ncol(ret)] = colnames(df)[3:ncol(df)]
    return(ret[,3:ncol(ret),drop=FALSE])
  })))
}

# cf. https://github.com/tidyverse/multidplyr/blob/cae3e5127053492bd6f6be51f69122d3a5796015/R/cluster-utils.R#L94
cluster_ls <- function(cluster) {
  cluster_call(cluster, ls(envir = globalenv()))
}

################################################################################
# Functions -- main omnibus test functions
################################################################################
omnibus_test_haplo <- function(x, pheno, obj_base, covars = NULL, covars_aa = NULL, subset = NULL, maf_threshold = NULL, output_composites=FALSE, regfun = logit) {
  n_sample <- length(pheno)
  if (is.null(covars_aa)) {
    phased <- t(x[, 2:(n_sample * 2 + 1)])
    colnames(phased) <- x$SNP
    n_haplo_prev <- 0
  } else {
    phased <- cbind(t(x[, 2:(n_sample * 2 + 1)]), covars_aa)
    colnames(phased) <- c(x$SNP, colnames(covars_aa))
    n_haplo_prev <- nrow(unique(covars_aa))
  }
  haplo <- get_haplo(phased)
  if (!is.null(subset)) {
    haplo <- haplo[subset, , drop = F]
    pheno <- pheno[subset]
    covars <- covars[subset, , drop = F]
  }
  freq <- colMeans(haplo) / 2

  if (ncol(haplo) == 1 | !all(rowSums(haplo) == 2) | (!is.null(maf_threshold) & sum(freq >= maf_threshold) <= 1)) {
    return(tibble(
      AA_ID = x$AA_ID[1],
      GENE = x$GENE[1],
      AA_POS = x$AA_POS[1],
      POS = x$POS[1],
      DOF = NA,
      SS = NA,
      PVALUE = NA,
      AA = str_c(colnames(phased), collapse = ","),
      N_AA = nrow(x),
      N_HAPLO = NA,
      N_HAPLO_ADD = NA,
      PctExp = NA,
      perHaplo = list(tibble(
        HAPLO = haplo_to_aa(colnames(haplo), colnames(phased), output_composites),
        FREQ = freq,
        BETA = NA,
        SE = NA,
        PVALUE = NA
      ))
    ))
  }

  if (!is.null(maf_threshold)) {
    haplo <- haplo[, freq >= maf_threshold]
    freq <- freq[freq >= maf_threshold]
    pheno[rowSums(haplo, na.rm = T) != 2] <- NA
  }
  obj_base <- glm_obj_base_haplo(pheno, covars, covars_aa, regfun = regfun)

  n_haplo <- ncol(haplo)
  # if all A haplotype is present, then make it reference, otherwise the most freq haplo is the reference
  if (any(str_detect(colnames(haplo), "^A+$"))) {
    most_freq <- which(str_detect(colnames(haplo), "^A+$"))
  } else {
    most_freq <- which.max(freq)
  }
  haplo <- haplo[, -most_freq, drop = FALSE]
  freq <- freq[-most_freq]
  if (is.null(covars)) {
    obj <- regfun(pheno ~ haplo)
  } else {
    obj <- regfun(pheno ~ haplo + covars)
  }

  ao <- anova(obj_base, obj)
  dof <- ao[2, "Df"]
  ss <- ao[2, "Sum of Sq"]
  ao2 <- anova(obj)
  if (identical(regfun, lm)) {
    pval <- ao[2, "Pr(>F)"]
    ao2$PctExp <- ao2[["Sum Sq"]] / sum(ao2[["Sum Sq"]])
  } else {
    ao <- anova(obj_base,obj, test="Chisq")
    ss <- ao[2, "Deviance"]
    pval <- ao[2, "Pr(>Chi)"]
    ao2$PctExp <- as.numeric(nagelkerke(obj, null = obj_base)$Pseudo.R.squared.for.model.vs.null[3])
  }
  coef <- as.data.frame(summary(obj)$coef)
  coef <- coef[2:(1 + ncol(haplo)), c(1, 2, 4)]
  rownames(coef) <- colnames(haplo)
  colnames(coef) <- c("BETA", "SE", "P")
  coef$haplo <- haplo_to_aa(rownames(coef), colnames(phased), output_composites)

  return(tibble(
    AA_ID = x$AA_ID[1],
    GENE = x$GENE[1],
    AA_POS = x$AA_POS[1],
    POS = x$POS[1],
    DOF = dof,
    SS = ss,
    PVALUE = pval,
    AA = str_c(colnames(phased), collapse = ","),
    N_AA = nrow(x),
    N_HAPLO = n_haplo,
    N_HAPLO_ADD = n_haplo - n_haplo_prev,
    PctExp = ao2$PctExp[1],
    perHaplo = list(tibble(
      HAPLO = coef$haplo,
      FREQ = freq,
      BETA = coef$BETA,
      SE = coef$SE,
      PVALUE = coef$P
    ))
  ))
}

# for parallel
omnibus_test_haplo_wrapper <- function(aa_phased, n_sample, pheno, obj_base, regfun = logit, covars = NULL, covars_aa = NULL, subset = NULL, maf_threshold = NULL, output_composites=FALSE, cluster = NULL) {
  obj_base <- glm_obj_base_haplo(pheno, covars, covars_aa)
  f <- function(x) {
    omnibus_test_haplo(
      x, pheno, obj_base,
      regfun = regfun,
      covars = covars,
      covars_aa = covars_aa,
      subset = subset,
      maf_threshold = maf_threshold,
      output_composites = output_composites
    )
  }

  if (is.null(cluster)) {
    ret <- aa_phased %>%
      group_by(AA_ID) %>%
      do(f(.)) %>%
      arrange(PVALUE, AA_ID)
  } else {
    # clean up cluster
    cluster_rm(cluster, cluster_ls(cluster)[[1]])

    cluster_copy(cluster, names = c(
      "f", "omnibus_test_haplo", "glm_obj_base_haplo", "get_haplo", "haplo_to_aa", "logit",
      "data.table", "tibble", "str_c", "str_split_fixed", "str_locate_all", "str_detect",
      "n_sample", "pheno", "obj_base", "regfun", "covars", "covars_aa", "subset", "maf_threshold"
    ))

    ret <- aa_phased %>%
      group_by(AA_ID) %>%
      partition(cluster = cluster) %>%
      do(f(.)) %>%
      collect() %>%
      arrange(PVALUE, AA_ID)
  }

  return(list(
    omnibus = ret %>% select(-starts_with("perHaplo")),
    perHaplo = ret %>% select(AA_ID, GENE, AA_POS, POS, perHaplo) %>% unnest(perHaplo)
  ))
}


################################################################################
# Main
################################################################################

main <- function(args) {
  n_files <- length(args$phased)

  data <- map(1:n_files, function(i) {
    print(sprintf("Reading files... %s", args$file[i]))
    read_data( args$phased[i], args$fam[i], args$pheno[i], args$covars[i])
  })

  AA_ID_splitted <- data.frame(str_split_fixed(data[[1]]$phased$SNP, "_", 5)) %>%
    mutate(X1 = if_else(X1 != "INDEL", X1, str_c(X1, X2, sep="_"))) %>%
    mutate(
      X2 = if_else(!str_detect(X1, "^INDEL_"), X2, X3),
      X3 = if_else(!str_detect(X1, "^INDEL_"), X3, X4),
      X4 = if_else(!str_detect(X1, "^INDEL_"), X4, X5),
      X5 = if_else(!str_detect(X1, "^INDEL_"), X5, "")
    )

  phased <- cbind(data.frame(SNP = data[[1]]$phased$SNP), map_dfc(1:n_files, function(i) {
    data[[i]]$phased[, 2:(2 * data[[i]]$n_sample + 1)]
  }))
  fam <- map_dfr(1:n_files, function(i) {
    data[[i]]$fam
  })
  n_sample <- sum(map_int(1:n_files, function(i) {
    data[[i]]$n_sample
  }))

  # sex <- as.factor(fam$SEX)
  #pop <- as.factor(fam$POP)
  covs <- get_covs(fam, map(data, "covs") %>% discard(is.null))
  print(sprintf("There are %d covariates in total.",ncol(covs)))

  phased <-
    phased %>%
    mutate(
      AA_ID = apply(AA_ID_splitted[, 1:4], 1, function(x) {
        str_c(x, collapse = "_")
      }),
      TYPE = case_when(
        str_detect(AA_ID_splitted[, 1], "^rs") | AA_ID_splitted[, 1] == "SNPS" ~ "SNPS",
        TRUE ~ AA_ID_splitted[, 1]
      ),
      GENE = AA_ID_splitted[, 2],
      AA_POS = AA_ID_splitted[, 3],
      POS = AA_ID_splitted[, 4]
    )

  if (args$aa_only) {
    phased <- phased %>% filter(TYPE == "AA")
  }

  if (args$exclude_composites) {
    phased <- phased %>% filter(!str_detect(SNP, "^AA_.*_[A-Zx]{2,}$"))
  }

  # remove samples by haplotype constructed using specific AAs (e.g., all HLA-B AAs)
  if (args$remove_samples_by_haplo) {
    aa_ids <- phased %>% filter(str_detect(AA_ID, args$remove_samples_aa_pattern)) %>% .$AA_ID
    haplo <- get_haplo(next_covars_haplo(phased, aa_ids, n_sample))
    haplo <- haplo[,colSums(haplo) < args$min_haplo_count, drop=F]

    if (ncol(haplo) > 0) {
      remove_idx <- rowSums(haplo) > 0
      print(sprintf("--remove-samples-by-haplo removes %d samples", sum(remove_idx)))
      phased <- phased[, -(2 * rep(which(remove_idx), each = 2) + c(0, 1)), drop=F]
      fam <- fam[!remove_idx,]
      #sex <- sex[!remove_idx]
      # pop <- pop[!remove_idx]
      covs <- covs[!remove_idx,]
      n_sample <- nrow(fam)
    }
  }

  if (args$n_threads > 1) {
    cluster <- new_cluster(args$n_threads)
  } else {
    cluster <- NULL
  }

  print("Starting regression...")
  # regression
  aa_phased <- phased
  covars <- covs

  if (args$omnibus) {
    print("Running omnibus test...")
    print(sprintf("Total inds %d",n_sample))
    if (!is.null(args$condition)) {
      covars_aa <- next_covars_haplo(aa_phased, args$condition, n_sample)
    } else {
      covars_aa <- NULL
    }

    # remove AA with only single pos from omnibus test
    aa_phased <-
      aa_phased %>%
      group_by(AA_ID) %>%
      filter(n() > 1) %>%
      ungroup()

    omnibus <-
      omnibus_test_haplo_wrapper(
        aa_phased, n_sample, fam$PHENO,
        covars = covars,
        covars_aa = covars_aa,
        maf_threshold = args$maf_threshold,
        output_composites = args$output_composites,
        cluster = cluster
      )

    write.table(omnibus[[1]], paste0(args$out, ".txt"), row.names = F, quote = F, sep = "\t")
    write.table(omnibus[[2]], paste0(args$out, ".haplo.txt"), row.names = F, quote = F, sep = "\t")
  }

  if (!is.null(args$exhaustive)) {
    print(sprintf("Running exhaustive testing for HLA-%s...", args$exhaustive))
    exhaustive_aa <- aa_phased %>%
      filter(TYPE == "AA" & GENE == args$exhaustive) %>%
      select(AA_ID, AA_POS) %>%
      mutate(AA_POS = as.numeric(AA_POS))

    uniq_aa_pos <- exhaustive_aa %>%
      .$AA_POS %>%
      unique() %>%
      sort()

    n_uniq_aa_pos <- length(uniq_aa_pos)
    print(sprintf("%d unique AA positions found in HLA-%s.", n_uniq_aa_pos, args$exhaustive))
    print(uniq_aa_pos)

    pheno <- as.numeric(as.character(fam$PHENO))
    
    obj_base <- glm_obj_base_haplo(pheno, covars, covars_aa=NULL)

    for (k in args$exhaustive_min_aa:args$exhaustive_max_aa) {
      outname <- paste0(args$out, ".exhaustive.", k, ".AA_", args$exhaustive)
      if (is.null(args$exhaustive_aa_pos)) {
        idx_comb <- data.frame(combn(n_uniq_aa_pos, k))
      } else {
        idx_comb <- combn(n_uniq_aa_pos, k - length(args$exhaustive_aa_pos))
        idx_comb <- rbind(matrix(rep(which(uniq_aa_pos %in% args$exhaustive_aa_pos), ncol(idx_comb)), ncol = ncol(idx_comb)), idx_comb)
        if (!args$exhaustive_no_filter) {
          for (i in (length(args$exhaustive_aa_pos) + 1):k) {
            idx_comb <- idx_comb[, idx_comb[i,] > max(which(uniq_aa_pos %in% args$exhaustive_aa_pos))]
          }
        }
        idx_comb <- data.frame(idx_comb)
        outname <- str_c(c(outname, args$exhaustive_aa_pos), collapse="_")
      }
      print(sprintf("%d combinations to test", ncol(idx_comb)))
      ret = mclapply(idx_comb, function(i) {
        comb_aa_id <- exhaustive_aa %>% filter(AA_POS %in% uniq_aa_pos[i]) %>% .$AA_ID %>% unique() %>% str_sort(numeric=TRUE)
        print(comb_aa_id)
        x <- aa_phased %>% filter(AA_ID %in% comb_aa_id) %>% mutate(AA_ID = str_c(comb_aa_id, collapse=','))
        omnibus <-
          omnibus_test_haplo(
            x, pheno, obj_base,
            regfun = lm,
            covars = covars,
            maf_threshold = args$maf_threshold,
            output_composites = args$output_composites
          ) %>%
          select(-starts_with("perHaplo")) %>%
          mutate(AA_POS = str_c(uniq_aa_pos[i], collapse=',')) %>%
          select(-POS)

        return(omnibus)
      }, mc.cores = args$n_threads)
      ret = do.call(rbind, ret) %>%
        arrange(-PctExp)

      write.table(ret, paste0(outname, ".txt"), row.names = F, quote = F, sep = "\t")
    }
  }
}


parser <- ArgumentParser()

parser$add_argument("--file", type = "character", nargs = "+")
parser$add_argument("--pop", type = "character", nargs = "+")
parser$add_argument("--phased", type = "character", nargs = "+")
parser$add_argument("--fam", type = "character", nargs = "+")
parser$add_argument("--bim", type = "character", nargs = "+")
parser$add_argument("--pheno", type = "character", nargs = "+")
parser$add_argument("--covars", type = "character", nargs = "+")
parser$add_argument("--maf-threshold", type = "double", default = 0.005)
parser$add_argument("--aa-only", action='store_true', help="Run association test only for AA changes")
parser$add_argument("--n-threads", "-n", type = "integer", default = 8)
parser$add_argument("--out", type = "character", required = TRUE)

parser$add_argument("--remove-samples-by-haplo", action='store_true', help="Remove samples based on haplotypes constructed using a subset of AAs")
parser$add_argument("--remove-samples-aa-pattern", type = "character")
parser$add_argument("--min-haplo-count", type = "integer", default = 10)

# omnibus test options
parser$add_argument("--omnibus", action='store_true', help="Run omnibus test")
parser$add_argument("--condition", type = "character", nargs = "+")
parser$add_argument("--condition-gene", type = "character", nargs = "+")
parser$add_argument("--exclude-composites", action='store_true', help="Exclude composite amino acids from omnibus test")
parser$add_argument("--output-composites", action='store_true', help="Output composite amino acids in .haplo.txt")

parser$add_argument("--exhaustive", type = "character", help="Run exhaustive testing of AA combinations in a single HLA gene")
parser$add_argument("--exhaustive-aa-pos", type = "integer", nargs = "+", help="Specify a AA pos of the first element(s) of combinations")
parser$add_argument("--exhaustive-min-aa", type = "integer", default=2, help="Minimum number of AA positions to form a combination")
parser$add_argument("--exhaustive-max-aa", type = "integer", default=2, help="Maximum number of AA positions to form a combination")
parser$add_argument("--exhaustive-no-filter", action='store_true', help="Don't filter non-increasing AA combinations")

args <- parser$parse_args()
print(args)

if (!args$omnibus & is.null(args$exhaustive)) {
  stop("At least --omnibus or --exhaustive should be specified")
}

if (args$exclude_composites & args$output_composites) {
  stop("Both --exclude-composites and --output-composites cannot be specified")
}

if (args$exhaustive_min_aa > args$exhaustive_max_aa) {
  stop("--exhaustive-max-aa should be greater than or equal to --exhaustive-min-aa")
}

if (length(args$exhaustive_aa_pos) >= args$exhaustive_min_aa) {
  stop("--exhaustive-min-aa should be greater than # AAs in --exhaustive-aa-pos")
}

if (is.null(args$file) & any(sapply(list(args$phased, args$fam, args$pheno), is.null))) {
  stop("Either --file or all of --phased, --fam, and --pheno should be specified.")
}

if (is.null(args$phased)) {
  args$phased <- paste0(args$file, ".bgl.phased")
}
if (is.null(args$fam)) {
  args$fam <- paste0(args$file, ".fam")
}
#if (is.null(args$bim)) {
#  args$bim <- paste0(args$file, ".bim")
#}
if (is.null(args$pheno)) {
  args$pheno <- paste0(args$file, ".pheno")
}
#if (is.null(args$sex)) {
#  args$sex <- paste0(args$file, ".sex")
#}
if (is.null(args$covars)) {
  args$covars <- paste0(args$file, ".covs")
}

if (!is.null(args$exhaustive_aa_pos)) {
  args$exhaustive_aa_pos <- sort(args$exhaustive_aa_pos)
}

main(args)
