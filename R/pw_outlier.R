#' @title pw_outlier
#' @name pw_outlier
#' @description Pairwise outlier removal of replicate samples within analytical batches. \cr
#' This is useful for example to identify technical errors, particularly when there is not extensive replication among samples to conduct more conventional outlier detection.
#' @param df The dataframe containing peak data. At minimum it should contain columns labeled: name, sample, batch, compound, area, rep, rep_tech \cr
#' Additional columns are OK, but will not be used.
#' @param qc.label \code{character()} Label designating the QC sample in the sample column of df. These will not be subjected to pairwise outlier detection since they can be useful downstream, for example with signal drift correction.
#' @param n.cores \code{numeric()} The number of cores to be used for processing if being run on a multi-core machine
#' @param mad_threshold \code{numeric()} The median absolute deviation (MAD) threshold to be used when identifying potential outliers
#' @param pw_threshold \code{numeric()} The pairwise outlier difference threshold. Should be between 0 and 1. Default is set to 0.95, meaning 5% of the data will be identified as a potential outlier.
#' @param peak_shrinkage \code{logical()} For samples surpassing the MAD threshold, should they be shrunk towards the distribution (TRUE) or completely removed from the analysis (FALSE). \cr
#' If set to TRUE (default), values are shrunk towards the distribution of samples and will maintain the same rankings of samples while reducing the overall distribution skew during the pairwise elimination calculation.
#' @param grouping_factor The column label containing the grouping factor from which pairwise differences will be calculated. Default is batch column.
#' @param return_plot \code{logical()} Should a density plot be returned showing the pairwise difference threshold distributions per compound. Default is FALSE.
#' @param plot_name \code{character()} If return_plot is TRUE, what is the name of the output .pdf file. This will be saved to the current working directory.
#' @return \code{list()} containing:
#' \itemize{
#' \item df (original input data)
#' \item df_cleaned (pairwise outlier cleaned data).
#' \item df_rm (samples removed by the pairwise outlier elimination).
#' }
#' @import dplyr tidyr
#' @importFrom data.table fread
#' @importFrom parallel mclapply
#' @export
#' @examples
#' pw_out = pw_outlier(
#' df = dat,
#' qc.label = "Pool",
#' n.cores = 2,
#' mad_threshold = 3,
#' pw_threshold = 0.95,
#' peak_shrinkage = TRUE,
#' grouping_factor = "batch",
#' return_plot = FALSE,
#' plot_name = "pw_outlier_plot")
#'
#' list2env(pw_out ,.GlobalEnv)
pw_outlier <- function(
  df = NULL,
  n.cores = 1,
  mad_threshold = 3,
  pw_threshold = 0.95,
  peak_shrinkage = TRUE,
  grouping_factor = "batch",
  return_plot = FALSE,
  plot_name = "pw_outlier_plot",
  qc.label = "Pool"){
  # check input df colnames -------------------------------------------------
  c_names = colnames(df)
  c_names_need = c("name","sample","batch","compound","area","rep","rep_tech")
  c_names_check = c_names_need[!c_names_need%in%c_names]
  if (length(c_names_check)>=1) {
    print(paste0("Please check your input df to make sure it contains all necessary columns... ",c_names_check, " missing"))
    stop()
  }
  m = df
  gp = grouping_factor
  gplab = m[,gp] %>%
    pull(all_of(gp))
  gplab = factor(gplab, levels = unique(gplab))
  m = m %>%
    mutate(g = gplab)
  # make lists to hold output -----------------------------------------------
  samps_rm = list()
  samps_plot = list()
  # additional df calculations by compound ----------------------------------
  m = m %>%
    group_by(sample, compound, rep_tech, g) %>%
    mutate(geno_tmp = paste0("TR",rep_tech,"_", sample), .before = 1) %>%
    ungroup() %>%
    group_by(geno_tmp, compound, rep_tech, g) %>%
    mutate(tmp_rep = 1:n(),
           uid = paste0(geno_tmp, "_", tmp_rep, "_", g, "_",compound)) %>%
    ungroup() %>%
    group_by(g, compound) %>%
    mutate(my_rank = rank(area, na.last = "keep", ties.method = "average"),
           my_quant = qnorm(my_rank/(n()+1)),
           my_mad = mad(area, na.rm = TRUE),
           my_thresh = all_of(mad_threshold)*my_mad + median(area, na.rm = TRUE),
           area_tmp = ifelse(area>=my_thresh, NA,area),
           area_shrink = ifelse(area>=my_thresh, max(area_tmp, na.rm = TRUE)+abs(my_quant),area))
  if (peak_shrinkage) {
    m = m %>%
      mutate(area_og = area,
             area = area_shrink)
  }else{
    m = m %>%
      mutate(area_og = area,
             area = area_tmp)
  }
  mog = m
  m = m %>%
    filter(!sample==all_of(qc.label))
  m_na = m %>%
    filter(is.na(area))
  m = m %>%
    drop_na(area)
  ## Change option to supress warning message when summarizing
  options(dplyr.summarise.inform = FALSE)
  ## Compounds in data
  integrated_compounds = sort(unique(m$compound))
  ## Calculate pairwise difference between genotypes separated by batch or other grouping factor.
  pw_fn = function(ii){
    all_compounds = integrated_compounds[ii]
    for (i in seq_along(all_compounds)) {
      ## Make an outlier shrikage while maintaining rankings
      my_df = m %>%
        filter(compound==all_of(all_compounds[i]))
      my_batches = unique(my_df$g)
      diff_vect = list()
      most_diff = list()
      plot_dat = list()
      for (j in seq_along(my_batches)) {
        current_batch = my_df %>%
          filter(g == all_of(my_batches[j]))
        my_genos = unique(current_batch$geno_tmp)
        diff_vect1 = list()
        most_diff1 = list()
        for (k in seq_along(my_genos)) {
          current_geno = current_batch %>%
            filter(geno_tmp==all_of(my_genos[k]))
          current_geno_name = unique(current_geno$geno_tmp)
          n_obs = nrow(current_geno)
          if (n_obs<3) {
            most_diff1[[k]] = NULL
            diff_vect1[[k]] = NULL
            # print(paste0("Genotype: ",my_genos[k], " in ", my_batches[j], " has less than three obs for compound: ", all_compounds[i]))
          }else{
            geno_diffs = as_tibble(as.matrix(dist(current_geno$area, method = "manhattan",upper = TRUE))) %>%
              rowwise() %>%
              mutate(sum_diffs = sum(c_across(cols = everything())))
            most_diff1[[k]] = as_tibble(geno_diffs) %>%
              mutate(tmp_rep = 1:n()) %>%
              mutate(geno_tmp = all_of(current_geno_name)) %>%
              slice_max(., order_by = sum_diffs, prop = 0.34) %>%
              mutate(compound = all_of(all_compounds[i])) %>%
              mutate(g = all_of(my_batches[j])) %>%
              mutate(uid = paste0(geno_tmp, "_", tmp_rep, "_", g))
            diff_vect1[[k]] = as_tibble(as.numeric(dist(current_geno$area, method = "manhattan"))) %>%
              mutate(compound = all_of(all_compounds[i])) %>%
              mutate(g = all_of(my_batches[j]))
          }
        }
        diff_vect[[j]] = bind_rows(diff_vect1) %>%
          group_by(compound, g) %>%
          summarise(thresh = quantile(value, pw_threshold, na.rm = TRUE))
        plot_dat[[j]] = bind_rows(diff_vect1) %>%
          left_join(., diff_vect[[j]], by = c("compound", "g"))
        most_diff[[j]] = bind_rows(most_diff1) %>%
          select(-c(geno_tmp, sum_diffs, tmp_rep)) %>%
          pivot_longer(!c("uid", "compound", "g"),names_to = "variable") %>%
          left_join(., diff_vect[[j]], by = c("compound", "g")) %>%
          filter(value>=thresh)
      }
      samps_rm[[i]] = bind_rows(most_diff)
      samps_plot[[i]] = bind_rows(plot_dat)
    }
    return(list(samps_rm = samps_rm,samps_plot = samps_plot))
  }
  pw_results = parallel::mclapply(seq_along(integrated_compounds), pw_fn, mc.cores = n.cores)
  samps_rm1 = lapply(pw_results,"[[", 1)
  samps_rm1 = bind_rows(samps_rm1) %>%
    select(-c(variable, value)) %>%
    distinct() %>%
    mutate(uid = paste0(uid, "_", compound))
  samps_plot1 = lapply(pw_results,"[[", 2)
  samps_plot1 = bind_rows(samps_plot1)
  if(return_plot){
    pdf(paste0(plot_name,".pdf"),width = 10,height = 10)
    for (i in seq_along(integrated_compounds)) {
      comp_thresh = samps_plot1 %>%
        filter(compound==all_of(integrated_compounds[i])) %>%
        select(-c(value)) %>%
        distinct()
      p = samps_plot1 %>%
        filter(compound==all_of(integrated_compounds[i])) %>%
        ggplot(., aes(x=value)) +
        geom_density(alpha=0.2) +
        facet_wrap(~ g) +
        geom_vline(aes(xintercept=thresh), colour = "#2c7fb8")+
        geom_label(data = comp_thresh, mapping = aes(x = -Inf, y = -Inf, label = paste0(pw_threshold, " threshold: ", round(thresh, 1) )), hjust = -0.79, vjust = -1.8, colour = "black")+
        theme_classic()+
        theme(
          axis.text.y = element_text(size = 12,face = "bold", color = "black"),
          axis.title.y = element_text(size = 12, color = "black", face = "bold"),
          axis.title.x = element_text(size = 12, color = "black", face = "bold"),
          strip.text.x = element_text(size = 9, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5,color = "black", face = "bold"),
          legend.key.width = unit(2,"cm"),
          legend.position="top",
          legend.box="vertical",
          legend.margin=margin(),
          legend.justification='center',
          # legend.title = element_blank(),
          legend.direction='horizontal'
        )+
        labs(x = "Replication pairwise differences", y = "Density",title = integrated_compounds[i])
      print(p)
    }
    dev.off()
  }
  # outputs to return -------------------------------------------------------
  dat_removed = m %>%
    ungroup() %>%
    filter(uid%in%samps_rm1$uid) %>%
    bind_rows(m_na) %>%
    mutate(area = area_og)
  m1 = mog %>%
    ungroup() %>%
    filter(!uid%in%dat_removed$uid) %>%
    mutate(area = area_og) %>%
    select(all_of(c_names))
  dat_removed = dat_removed %>%
    select(all_of(c_names))
  dd <- c(" ", " -------------- ", " Done! ", " --------------")
  cat(dd, sep = "\n")
  return(list(df = df, df_cleaned = m1, df_rm = dat_removed))
}
