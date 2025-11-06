#!/usr/bin/Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

impacts_file <- "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/snpEff/varscan/varscan_variant_impacts.txt"
matrix_file  <- "repeatfilt_sc150_splitgenome_VERYrelaxedAnc_strictercall_readper_table_orderedbypos.txt"
out_plot <- "stackedbar_variant_counts_by_impact.pdf"
out_plot2 <- "allsamples_stackedbar_variant_counts_by_impact.pdf"
thresh <- 10
impact_levels <- c("HIGH", "MODERATE", "LOW", "MODIFIER", "INTERGENIC")


impacts <- read.table(impacts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
#colnames(impacts) <- c("variant", "impact")
#impacts$variant <- trimws(impacts$variant)
#impacts$impact  <- trimws(impacts$impact)
#impacts$impact  <- factor(impacts$impact, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))

mat <- read.table(matrix_file, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
blacklists <- c("scaffold_3_pilon:5389696","scaffold_6_pilon:84121","scaffold_3_pilon:5788746","scaffold_3_pilon:4067561","scaffold_5_pilon:3302727","scaffold_6_pilon:1734032","scaffold_3_pilon:5608240")
mat <- mat[!mat$variant %in% blacklists,]
if (!"variant" %in% colnames(mat)) {
    colnames(mat)[1] <- "variant"
}
all_cols <- setdiff(colnames(mat), "variant")

is_prefix <- function(x, pref) startsWith(x, pref)
ath_cols <- all_cols[is_prefix(all_cols, "Ath_")]
hvu_cols <- all_cols[is_prefix(all_cols, "Hvu_")]
msly_cols <- all_cols[is_prefix(all_cols, "mSly_")]
sys_cols <- all_cols[is_prefix(all_cols, "system_") | is_prefix(all_cols, "Sys_")]  # accept both

to_num_cols <- function(df, cols) {
    for (nm in cols) {
	suppressWarnings(df[[nm]] <- as.numeric(df[[nm]]))
    }
    df
}
mat <- to_num_cols(mat, all_cols)

present_from_cols <- function(df, cols, thr) {
    if (length(cols) == 0) return(rep(FALSE, nrow(df)))
    m <- as.matrix(df[, cols, drop = FALSE])
    rs <- rowSums(m > thr, na.rm = TRUE)
    rs > 0
}

pres_Ath  <- present_from_cols(mat,ath_cols,thresh)
pres_Hvu  <- present_from_cols(mat,hvu_cols,thresh)
pres_mSly <- present_from_cols(mat,msly_cols,thresh)
pres_Sys  <- present_from_cols(mat,sys_cols,thresh)

presence_df <- tibble(
		      variant = mat$variant,
		      Ath  = as.integer(pres_Ath),
		      Hvu  = as.integer(pres_Hvu),
		      mSly = as.integer(pres_mSly),
		      Sys  = as.integer(pres_Sys)
		      ) |> 
pivot_longer(cols = -variant, names_to = "condition", values_to = "present") |>
filter(present == 1L) |>
select(-present)

joined <- presence_df |>
    left_join(impacts, by = "variant") |>
    mutate(
	   impact = coalesce(impact, "INTERGENIC"),
	   impact = toupper(trimws(impact)),
	   impact = factor(impact, levels = impact_levels)
    )


counts <- joined |>
    count(condition, impact, name = "count") |>
    complete(
	     condition = c("Ath","Hvu","mSly","Sys"),
	     impact = factor(impact_levels, levels = impact_levels),
	     fill = list(count = 0L)
    )

desired_order <- c("Ath","mSly","Hvu","Sys")
counts$condition <- factor(counts$condition, levels = rev(desired_order))

#plotting
p <- ggplot(counts, aes(y = condition, x = count, fill = impact)) +
    geom_col(color = "grey20", width = 0.7) +
    scale_fill_manual(
		      values = c(HIGH = "#ff1493",
				 MODERATE = "#ffa500",
				 LOW = "#00ff00",
				 MODIFIER = "#1e90ff",
				 INTERGENIC="#c3b091"),
		      drop = FALSE) +
    labs(x = "Number of variants", y = "Condition", fill = "Impact",
	 title = paste0("Variant counts by impact per condition (threshold > ", thresh, ")")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right")

#ggsave(out_plot, p, width = 8, height = 4.5, dpi = 300)


mat_long <- mat |>
    tidyr::pivot_longer(cols = all_of(all_cols), names_to = "sample", values_to = "value") |>
    dplyr::mutate(present = value > thresh) |>
    dplyr::filter(present) |>
    dplyr::select(variant, sample)

joined <- mat_long |>
    dplyr::left_join(impacts, by = "variant") |>
    dplyr::mutate(
		  impact = dplyr::coalesce(impact, "INTERGENIC"),
		  impact = toupper(trimws(impact))
    )

counts <- joined |>
    dplyr::count(sample, impact, name = "count")

all_samples <- tibble::tibble(sample = all_cols) |>
    dplyr::mutate(
		  condition = dplyr::case_when(
					       startsWith(sample, "Ath_") ~ "Ath",
					       startsWith(sample, "mSly_") ~ "mSly",
					       startsWith(sample, "Hvu_") ~ "Hvu",
					       startsWith(sample, "system_") | startsWith(sample, "Sys_") ~ "system", TRUE ~ "other"
		  )
    )

grid <- tidyr::expand_grid(sample = all_samples$sample,
                           impact = impact_levels)

counts_full <- grid |>
    dplyr::left_join(counts, by = c("sample","impact")) |>
    dplyr::mutate(count = dplyr::coalesce(count, 0L)) |>
    dplyr::left_join(all_samples, by = "sample") |>
    dplyr::mutate(impact = factor(impact, levels = impact_levels))

counts_full <- counts_full |>
    dplyr::mutate(
		  traj = stringr::str_extract(sample, "T[1-5]$"),
		  generation = stringr::str_extract(sample, "G\\d{2}")
    )
counts_full$traj <- factor(counts_full$traj, levels = c("T1","T2","T3","T4","T5"))
counts_full$generation <- factor(counts_full$generation, levels = c("G02","G04","G06","G08","G10","G12"))
counts_full$condition <- factor(counts_full$condition, levels = c("Ath","mSly","Hvu","system"))
# Order samples: Ath, mSly, Hvu, system; preserve original within-group order
#order_samples <- c("Ath_G02_T1","Ath_G04_T1","Ath_G06_T1","Ath_G08_T1","Ath_G10_T1","Ath_G12_T1","Ath_G02_T2","Ath_G04_T2","Ath_G06_T2","Ath_G08_T2","Ath_G10_T2","Ath_G12_T2","Ath_G02_T3","Ath_G04_T3","Ath_G06_T3","Ath_G08_T3","Ath_G10_T3","Ath_G12_T3","Ath_G02_T4","Ath_G04_T4","Ath_G06_T4","Ath_G08_T4","Ath_G10_T4","Ath_G12_T4","Ath_G02_T5","Ath_G04_T5","Ath_G06_T5","Ath_G08_T5","Ath_G10_T5","Ath_G12_T5","mSly_G02_T1","mSly_G04_T1","mSly_G06_T1","mSly_G08_T1","mSly_G10_T1","mSly_G12_T1","mSly_G02_T2","mSly_G04_T2","mSly_G06_T2","mSly_G08_T2","mSly_G10_T2","mSly_G12_T2","mSly_G02_T3","mSly_G04_T3","mSly_G06_T3","mSly_G08_T3","mSly_G10_T3","mSly_G12_T3","mSly_G02_T4","mSly_G04_T4","mSly_G06_T4","mSly_G08_T4","mSly_G10_T4","mSly_G12_T4","mSly_G02_T5","mSly_G04_T5","mSly_G06_T5","mSly_G08_T5","mSly_G10_T5","mSly_G12_T5","Hvu_G02_T1","Hvu_G04_T1","Hvu_G06_T1","Hvu_G08_T1","Hvu_G10_T1","Hvu_G12_T1","Hvu_G02_T2","Hvu_G04_T2","Hvu_G06_T2","Hvu_G08_T2","Hvu_G10_T2","Hvu_G12_T2","Hvu_G02_T3","Hvu_G04_T3","Hvu_G06_T3","Hvu_G08_T3","Hvu_G10_T3","Hvu_G12_T3","Hvu_G02_T4","Hvu_G04_T4","Hvu_G06_T4","Hvu_G08_T4","Hvu_G10_T4","Hvu_G12_T4","Hvu_G02_T5","Hvu_G04_T5","Hvu_G06_T5","Hvu_G08_T5","Hvu_G10_T5","Hvu_G12_T5","system_G02_T1","system_G04_T1","system_G06_T1","system_G08_T1","system_G10_T1","system_G12_T1","system_G02_T2","system_G04_T2","system_G06_T2","system_G08_T2","system_G10_T2","system_G12_T2","system_G02_T3","system_G04_T3","system_G06_T3","system_G08_T3","system_G10_T3","system_G12_T3","system_G02_T4","system_G04_T4","system_G06_T4","system_G08_T4","system_G10_T4","system_G12_T4","system_G02_T5","system_G04_T5","system_G06_T5","system_G08_T5","system_G10_T5","system_G12_T5")
#sample_levels <- c(ath_cols, msly_cols, hvu_cols, sys_cols)
#counts_full$sample <- factor(counts_full$sample, levels = order_samples)
counts_gen <- counts_full |>
    dplyr::group_by(condition, traj, generation, impact) |>
    dplyr::summarise(count = sum(count), .groups = "drop")



# plotting, in a grid setting, faceted by condition
p <- ggplot(counts_gen, aes(x = generation, y = count, fill = impact)) +
    geom_col(width = 0.8, color = "grey20") +
    facet_grid(rows = vars(condition), cols = vars(traj), scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_manual(
		      values = c(HIGH = "#ff1493",
				 MODERATE = "#ffa500",
				 LOW = "#00ff00",
				 MODIFIER = "#1e90ff",
				 INTERGENIC="#c3b091"),
		      drop = FALSE
		      ) +
    labs(x = "Sample", y = "Number of variants", fill = "Impact",
	 title = paste0("Variant counts by impact per sample (threshold > ", thresh, ")")) +
   theme_classic(base_size = 12) +
   theme(
      axis.text.x = element_text(angle = 0, vjust = 0.6),
      panel.spacing.x = grid::unit(0.2, "lines"),
      panel.spacing.y = grid::unit(1.2, "lines"),
      strip.placement = "outside",
      strip.background = element_blank(),
      legend.position = "right"
      )+
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

ggsave(out_plot2, p, width = 10, height = 10, dpi = 300)
message("Saved plot to: ", out_plot2)
