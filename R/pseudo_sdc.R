#' @title pseudo_sdc
#' @name pseudo_sdc
#' @description Signal drift correction using QC samples present in some batches but absent in others.
#' @param df The dataframe containing peak data. At minimum should contain columns labeled: name, sample, batch, compound, area, experiment_index, batch_index. \cr
#' @param n.cores \code{numeric()} The number of cores to be used for processing if being run on a multi-core machine
#' @param train.batch \code{character()} The batch name in the df which contains QC samples (should only be one batch which will be used for training regression spline model).
#' @param test.breaks \code{numeric()} A numeric vector indicating the number of equal sized sub-batches for the train.batch to be divided into and tested.
#' @param test.window \code{numeric()} A numeric vector containing the sizes of sliding windows to test when performing the sliding window median calculation.
#' @param test.index \code{numeric()} A numeric vector containing the injection position offset for pseudo QC inclusion in the peak data matrix.
#' @param criteria \code{character()} What criteria should be minimized when determining the optimal set of parameters. Should be one of:
#' \itemize{
#' \item "RSD" relative standard deviation assuming a Gaussian distribution of errors
#' \item "RSD_robust" relative standard deviation assuming a non-Gaussian distribution of errors
#' \item "MSE" mean squared error
#' \item "TSS" total sum of squares
#' }
#' @param qc.label \code{character()} Label designating the QC sample in the sample column of df.
#' @param min.qc \code{numeric()} The minimum number of pseudo-QC samples to consider during model training. Should be a value greater than 2.
#' @param quantile.increment \code{numeric()} Incremental step for qunatiles of peak areas to retain in training model.
#' @param log_transform \code{logical()} TRUE(default)/FALSE should data be log transformed
#' @param mad_outlier \code{logical()} TRUE(default)/FALSE should median absolute deviation (MAD) based outliers be excluded from signal drift calculation
#' @param mad_threshold \code{numeric()} How many MAD from the median a value can be before considered an outlier. Default is 3.
#' @return \code{list()} containing:
#' \itemize{
#' \item df (original input data)
#' \item df_pseudoQC (data with pseudoQC calculated samples included). Includes an additional column labeled 'class' which categorizes true QC, Sample, Pseudo_QC samples.
#' \item df_pseudoQC_corrected (signal drift corrected data using pseudoQC samples). Same columns as df_pseudoQC returned, with an aditional 'area_corrected' column designating the signal drift corrected data.
#' \item criteria_table (table with results for criteria applied along with the others not-used).
#' }
#' @import dplyr pmp zoo
#' @importFrom data.table fread
#' @importFrom parallel mclapply
#' @export
#' @examples
#' sim_dat = simulate_data(db_ids = "FIO00738",
#'                         nsamps_per_batch = 100,
#'                         xls_file_name = system.file("extdata", "Index.xls", package = "pseudoDrift"),
#'                         valid_sdf_file = system.file("extdata", "valid-test.sdf", package = "pseudoDrift"))
#'
#' df = sim_dat[["t4_sim_mat"]][[1]]
#'
#' sdc_out = pseudo_sdc(df = df,
#'                      train.batch = "B3",
#'                      test.breaks = seq(2,3,1),
#'                      test.window = seq(1,3,2),
#'                      test.index = seq(2,3,1),
#'                      qc.label = "QC",
#'                      min.qc = 2)
#' list2env(sdc_out ,.GlobalEnv)

pseudo_sdc <- function(
  df = NULL,
  n.cores = 1,
  train.batch = NULL,
  test.breaks = NULL,
  test.window = NULL,
  test.index = NULL,
  criteria = "RSD",
  qc.label = NULL,
  qc.multibatch = FALSE,
  min.qc = 5,
  quantile.increment = 1,
  log_transform = TRUE,
  mad_outlier = TRUE,
  mad_threshold = 3){

  nc = 1:length(unique(df$compound))
  if (length(nc)>1) {
    print(paste0("df contains more than one compound. Running for ", max(nc), " compounds"))
  }
  comps_within = unique(df$compound)

  fun_within = function(x){
    df = df %>%
      filter(compound==all_of(x)) %>%
      group_by(batch, compound) %>%
      mutate(class = ifelse(grepl(paste0(qc.label),sample), "QC", "Sample"),
             index = batch_index,
             my_mad = mad(area, na.rm = TRUE),
             high_thresh = median(area, na.rm = TRUE) + all_of(mad_threshold)*my_mad,
             low_thresh = median(area, na.rm = TRUE) - all_of(mad_threshold)*my_mad,
             area_tmp = ifelse((area>=high_thresh | area<=low_thresh),NA,area))

    m = df %>%
      filter(batch%in%train.batch) %>%
      group_by(batch, compound) %>%
      mutate(my_mad = mad(area, na.rm = TRUE),
             high_thresh = median(area, na.rm = TRUE) + all_of(mad_threshold)*my_mad,
             low_thresh = median(area, na.rm = TRUE) - all_of(mad_threshold)*my_mad,
             area_tmp = ifelse((area>=high_thresh | area<=low_thresh),NA,area))

    if (mad_outlier) {
      m = m %>%
        mutate(area_og = area,
               area = area_tmp)
      df = df %>%
        mutate(area_og = area,
               area = area_tmp)
    }else{
      m = m %>%
        mutate(area_og = area,
               area = area)
      df = df %>%
        mutate(area_og = area,
               area = area)
    }

    vals_keep = m %>%
      group_by(batch) %>%
      mutate(pool_rank = percent_rank(area)) %>%
      filter(class == "QC")
    vals_keep = ceiling(range(vals_keep$pool_rank, na.rm = TRUE)/quantile.increment)*quantile.increment

    # Using what portion of the data and what partitioning of the batch gives best fit for the training data
    test_low = seq(0,vals_keep[1],quantile.increment)
    test_high = seq(vals_keep[2],1,quantile.increment)
    n_breaks = test.breaks
    k = test.window
    ind = test.index
    dat_portion = crossing(var1 = test_low,
                           var2 = test_high,
                           var3 = n_breaks,
                           var4 = k,
                           var5 = ind) %>%
      rowwise() %>%
      mutate(var1 = pmin(var1, var2),
             var2 = pmax(var1, var2)) %>%
      distinct() %>%
      filter(!var1==var2) %>%
      filter(!abs(var1-var2)<=0.50) %>%
      mutate(in_out = "inside")
    dat_portion1 = dat_portion %>%
      filter(!var1 == 0) %>%
      filter(!var2 == 1) %>%
      mutate(in_out = "outside") %>%
      bind_rows(dat_portion,.)
    print(paste0("Running...testing ", nrow(dat_portion1), " combinations of input parameters"))

    # QC-RSC with true QCs ----------------------------------------------------
    m_qcrsc = function(x, y, r){
      t_meta = colnames(x)
      t_meta = t_meta[!t_meta%in%c("name", "compound")]
      tqc = x %>%
        pivot_wider(id_cols = !all_of(t_meta), names_from = name, values_from = area)
      t_rn = tqc$compound
      tqc = as.matrix(tqc[,-1])
      rownames(tqc) = t_rn
      mc = suppressMessages(pmp::QCRSC(df=tqc,
                                       spar_lim = c(-2,2),
                                       minQC = 4,
                                       order=seq_along(x$index),
                                       batch=x$batch,
                                       classes=x$class,
                                       qc_label = paste0(y),
                                       log = log_transform))
      z = x %>%
        mutate(area_corrected = mc[1,], .before = area)

      m_qc = z %>%
        filter(class == "QC") %>%
        mutate(rsd_robust = mad(area, na.rm = TRUE)/abs(median(area, na.rm = TRUE)),
               rsd_tqc_robust = mad(area_corrected, na.rm = TRUE)/abs(median(area_corrected, na.rm = TRUE))) %>%
        mutate(rsd = sd(area, na.rm = TRUE)/abs(mean(area, na.rm = TRUE)),
               rsd_tqc = sd(area_corrected, na.rm = TRUE)/abs(mean(area_corrected, na.rm = TRUE)))

      if (r=="yes") {
        return(z)
      }else{
        return(m_qc)
      }
    }
    trueQC = m_qcrsc(m, "QC","nah")
    if (qc.multibatch) {
      trueQC_full = m_qcrsc(df, "QC","nah")
    }

    # pseudoQC-RSC ------------------------------------------------------------
    ssr = list()
    ## Training function for batch with QC samples available. Choose parameters which reduce total sum of sq error based on training batch
    m_fn = function(j, my_return){
      m_min = dat_portion1[j,1] %>% pull(var1)
      m_max = dat_portion1[j,2] %>% pull(var2)
      m_bre = dat_portion1[j,3] %>% pull(var3)
      m_kkk = dat_portion1[j,4] %>% pull(var4)
      m_ind = dat_portion1[j,5] %>% pull(var5)
      m_in_out = dat_portion1[j,6] %>% pull(in_out)

      if(m_in_out=="inside"){
        pseudo_qc = m %>%
          filter(!class=="QC") %>%
          mutate(vals_rank = percent_rank(area)) %>%
          filter(between(vals_rank, all_of(m_min), all_of(m_max))) %>%
          mutate(batch_position = ntile(index, all_of(m_bre))) %>%
          group_by(batch_position,batch) %>%
          mutate(n_per = n(),
                 index = min(index),
                 m_kkk = m_kkk)
        np = pseudo_qc %>%
          filter(n_per<m_kkk)
        np1 = list()
        for (i in unique(np$n_per)) {
          m_kkk1 = as.numeric(i)
          m_kkk1_test = 2*floor(m_kkk1/2)+1
          if(m_kkk1_test>=m_kkk1){
            m_kkk1_test = m_kkk1_test-2
          }
          if (m_kkk1_test<0) {
            m_kkk1_test = 1
          }
          np1[[i]] = np %>%
            filter(n_per==m_kkk1) %>%
            mutate(area = zoo::rollmedian(area, k = m_kkk1_test, fill = NA, align = "left"),
                   area = mean(area, na.rm = TRUE),
                   class = "Pseudo_QC"
            ) %>%
            ungroup() %>%
            select(area, index, class, batch, compound) %>%
            distinct() %>%
            drop_na(area)
        }
        np1 = bind_rows(np1)
        m_kkk2 = m_kkk
        pseudo_qc = pseudo_qc %>%
          filter(!name%in%np$name) %>%
          mutate(area = zoo::rollmedian(area, k = m_kkk2, fill = NA, align = "left"),
                 area = mean(area, na.rm = TRUE),
                 class = "Pseudo_QC"
          ) %>%
          ungroup() %>%
          select(area, index, class, batch, compound) %>%
          distinct() %>%
          drop_na(area) %>%
          bind_rows(., np1)
      }else{
        pseudo_qc = m %>%
          filter(!class=="QC") %>%
          mutate(vals_rank = percent_rank(area)) %>%
          filter(!between(vals_rank, all_of(m_min), all_of(m_max))) %>%
          mutate(batch_position = ntile(index, all_of(m_bre))) %>%
          group_by(batch_position,batch) %>%
          mutate(n_per = n(),
                 index = min(index),
                 m_kkk = m_kkk)
        np = pseudo_qc %>%
          filter(n_per<m_kkk)
        np1 = list()
        for (i in unique(np$n_per)) {
          m_kkk1 = as.numeric(i)
          m_kkk1_test = 2*floor(m_kkk1/2)+1
          if(m_kkk1_test>=m_kkk1){
            m_kkk1_test = m_kkk1_test-2
          }
          if (m_kkk1_test<0) {
            m_kkk1_test = 1
          }
          np1[[i]] = np %>%
            filter(n_per==m_kkk1) %>%
            mutate(area = zoo::rollmedian(area, k = m_kkk1_test, fill = NA, align = "left"),
                   area = mean(area, na.rm = TRUE),
                   class = "Pseudo_QC"
            ) %>%
            ungroup() %>%
            select(area, index, class, batch, compound) %>%
            distinct() %>%
            drop_na(area)
        }
        np1 = bind_rows(np1)
        m_kkk2 = m_kkk
        pseudo_qc = pseudo_qc %>%
          filter(!name%in%np$name) %>%
          mutate(area = zoo::rollmedian(area, k = m_kkk2, fill = NA, align = "left"),
                 area = mean(area, na.rm = TRUE),
                 class = "Pseudo_QC"
          ) %>%
          ungroup() %>%
          select(area, index, class, batch, compound) %>%
          distinct() %>%
          drop_na(area) %>%
          bind_rows(., np1)
      }
      if (nrow(pseudo_qc)>=min.qc) {
        tmp_pseudo_qc = pseudo_qc %>%
          mutate(index = index+m_ind+0.05)
        pseudo_qc = pseudo_qc %>%
          mutate(index = index+m_ind+0.10) %>%
          bind_rows(., tmp_pseudo_qc) %>%
          arrange(index) %>%
          mutate(name = paste0("Pseudo_QC", 1:n()))
        if(nrow(pseudo_qc)<4){
          tmp_pseudo_qc0 = tmp_pseudo_qc %>% mutate(index = index+0.05)
          tmp_pseudo_qc00 = tmp_pseudo_qc %>% mutate(index = index+0.10)
          pseudo_qc = pseudo_qc %>%
            bind_rows(., tmp_pseudo_qc0) %>%
            bind_rows(., tmp_pseudo_qc00) %>%
            arrange(index) %>%
            mutate(name = paste0("Pseudo_QC", 1:n()))
        }
        mm = m %>%
          bind_rows(pseudo_qc) %>%
          arrange(index)
        pseudoQC = m_qcrsc(mm, "Pseudo_QC", "nah")
        ssr = tibble(pred = pseudoQC$area_corrected,
                     obs = trueQC$area_corrected) %>%
          mutate(sse = (obs-pred)^2) %>%
          summarise(MSE = mean(sse, na.rm = TRUE),
                    TSS = sum(sse, na.rm = TRUE)) %>%
          mutate(RSD = unique(pseudoQC$rsd_tqc),
                 RSD_robust = unique(pseudoQC$rsd_tqc_robust))
      }else{
        ssr = tibble(MSE = NA,
                     TSS = NA,
                     RSD = NA,
                     RSD_robust = NA)
      }
      if (my_return=="yes") {
        return(list(m_min, m_max, m_bre, m_kkk, m_ind, m_in_out))
      }else{
        return(ssr)
      }
    }
    if (n.cores==1) {
      m_tuning = lapply(seq_along(dat_portion1$var1), m_fn, "no")
    }else{
      m_tuning = parallel::mclapply(seq_along(dat_portion1$var1), m_fn, "no", mc.cores = n.cores)
    }
    m_eval = m_tuning %>% map_dfr(~ .x %>% as_tibble(), .id = "name")
    m_min_ssr_keep = m_eval %>%
      select(name, all_of(criteria))
    colnames(m_min_ssr_keep)[2] = "value"
    m_min_ssr_keep = m_min_ssr_keep %>%
      slice_min(order_by = value, n = 1)

    if (criteria == "RSD" | criteria == "RSD_robust") {
      rsd_raw = m %>%
        group_by(compound) %>%
        filter(sample == all_of(qc.label)) %>%
        summarise(rsd_robust = mad(area, na.rm = TRUE)/abs(median(area, na.rm = TRUE)),
                  rsd = sd(area, na.rm = TRUE)/abs(mean(area, na.rm = TRUE)))

      if (criteria == "RSD") {
        rsd_red = rsd_raw$rsd - m_min_ssr_keep$value
      }
      if (criteria == "RSD_robust") {
        rsd_red = rsd_raw$rsd_robust - m_min_ssr_keep$value
      }

      if (rsd_red<=0) {
        print(paste0(criteria, " NOT reduced, but instead increased by ", round(-1*rsd_red, 2), " Consider changing parameters or not applying correction to this compound"))
      }
      print(paste0(criteria, " reduced by ", round(rsd_red, 2), " relative to no correction"))
    }

    if (nrow(m_min_ssr_keep)==0) {
      stop(" No model could be fit using those parameters. \n
           Try adjusting min.qc or test.breaks parameter. \n
           min.qc should be >= max(test.breaks)")
    }
    m_min_ssr = m_min_ssr_keep[1,1]
    ## Minimums with all criteria
    mse_min = m_eval %>%
      select(name, MSE) %>%
      slice_min(order_by = MSE, n=1) %>%
      left_join(.,  dat_portion1 %>% rownames_to_column(), by = c("name"="rowname")) %>%
      select(-c(name)) %>%
      mutate(criteria = "MSE")
    colnames(mse_min)[1] = "value"
    tss_min = m_eval %>%
      select(name, TSS) %>%
      slice_min(order_by = TSS, n=1) %>%
      left_join(.,  dat_portion1 %>% rownames_to_column(), by = c("name"="rowname")) %>%
      select(-c(name)) %>%
      mutate(criteria = "TSS")
    colnames(tss_min)[1] = "value"
    rsd_min = m_eval %>%
      select(name, RSD) %>%
      slice_min(order_by = RSD, n=1) %>%
      left_join(.,  dat_portion1 %>% rownames_to_column(), by = c("name"="rowname")) %>%
      select(-c(name)) %>%
      mutate(criteria = "RSD")
    colnames(rsd_min)[1] = "value"
    Rrsd_min = m_eval %>%
      select(name, RSD_robust) %>%
      slice_min(order_by = RSD_robust, n=1) %>%
      left_join(.,  dat_portion1 %>% rownames_to_column(), by = c("name"="rowname")) %>%
      select(-c(name)) %>%
      mutate(criteria = "RSD_robust")
    colnames(Rrsd_min)[1] = "value"
    all_min = bind_rows(mse_min, tss_min, rsd_min, Rrsd_min)
    colnames(all_min) = c("value","quantile_min","quantile_max", "test.breaks","test.window","test.index","in_out","criteria")
    m_fit = lapply(as.numeric(m_min_ssr$name), m_fn, "yes")
    m_min = m_fit[[1]][[1]]
    m_max = m_fit[[1]][[2]]
    m_bre = m_fit[[1]][[3]]
    m_kkk = m_fit[[1]][[4]]
    m_ind = m_fit[[1]][[5]]
    m_in_out = m_fit[[1]][[6]]
    metab = unique(m$compound)
    all_min = all_min %>%
      mutate(data.train = train.batch,
             compound = all_of(metab))
    test_dat = m %>%
      filter(!class=="QC")
    n_samps = nrow(test_dat)
    print(paste0("The best fit for ", metab, " in training batch ", train.batch, " is:"))
    if (m_kkk == 1) {
      print(paste0("With MEAN (opposed to median) smoothing every ", ceiling(n_samps/m_bre), " samples to generate pseudo-pools."))
    }else{
      print(paste0("With median smoothing every ", m_kkk, " samples to generate pseudo-pools."))
    }
    print(paste0("Index (injection order) offset by ",m_ind, ", and splitting the batch into ", m_bre, " sections."))
    if(as.numeric(m_min_ssr$name)<=nrow(dat_portion)){
      print(paste0("Using peak area values between quantiles ", m_min, " and ", round(m_max,2)))
    }else{
      print(paste0("Using peak area values outside the quantiles ", m_min, " and ", round(m_max,2)))
    }
    ## Get pseudo-QCs from other batches based on training batch
    batch_split = df %>%
      filter(!class=="QC") %>%
      group_by(batch) %>%
      summarise(samps = n_distinct(name)) %>%
      mutate(bre = round(samps*all_of(m_bre)/all_of(n_samps)))

    if(m_in_out=="inside"){
      pseudo_qc1 = df %>%
        filter(!class=="QC") %>%
        left_join(., batch_split, by = "batch") %>%
        group_by(batch) %>%
        mutate(vals_rank = percent_rank(area)) %>%
        filter(between(vals_rank, all_of(m_min), all_of(m_max))) %>%
        mutate(batch_position = ntile(index, bre)) %>%
        group_by(batch_position,batch) %>%
        mutate(n_per = n(),
               index = min(index),
               m_kkk = m_kkk)
      np = pseudo_qc1 %>%
        filter(n_per<m_kkk)
      np1 = list()
      for (i in unique(np$n_per)) {
        m_kkk1 = as.numeric(i)
        m_kkk1_test = 2*floor(m_kkk1/2)+1
        if(m_kkk1_test>=m_kkk1){
          m_kkk1_test = m_kkk1_test-2
        }
        if (m_kkk1_test<0) {
          m_kkk1_test = 1
        }
        np1[[i]] = np %>%
          filter(n_per==m_kkk1) %>%
          mutate(area = zoo::rollmedian(area, k = m_kkk1_test, fill = NA, align = "left"),
                 area = mean(area, na.rm = TRUE),
                 class = "Pseudo_QC"
          ) %>%
          ungroup() %>%
          select(area, index, class, batch, compound) %>%
          distinct() %>%
          drop_na(area)
      }
      np1 = bind_rows(np1)
      m_kkk2 = m_kkk
      pseudo_qc1 = pseudo_qc1 %>%
        filter(!name%in%np$name) %>%
        mutate(area = zoo::rollmedian(area, k = m_kkk2, fill = NA, align = "left"),
               area = mean(area, na.rm = TRUE),
               class = "Pseudo_QC"
        ) %>%
        ungroup() %>%
        select(area, index, class, batch, compound) %>%
        distinct() %>%
        drop_na(area) %>%
        bind_rows(., np1)
    }else{
      pseudo_qc1 = df %>%
        filter(!class=="QC") %>%
        left_join(., batch_split, by = "batch") %>%
        group_by(batch) %>%
        mutate(vals_rank = percent_rank(area)) %>%
        filter(!between(vals_rank, all_of(m_min), all_of(m_max))) %>%
        mutate(batch_position = ntile(index, bre)) %>%
        group_by(batch_position,batch) %>%
        mutate(n_per = n(),
               index = min(index),
               m_kkk = m_kkk)
      np = pseudo_qc1 %>%
        filter(n_per<m_kkk)
      np1 = list()
      for (i in unique(np$n_per)) {
        m_kkk1 = as.numeric(i)
        m_kkk1_test = 2*floor(m_kkk1/2)+1
        if(m_kkk1_test>=m_kkk1){
          m_kkk1_test = m_kkk1_test-2
        }
        if (m_kkk1_test<0) {
          m_kkk1_test = 1
        }
        np1[[i]] = np %>%
          filter(n_per==m_kkk1) %>%
          mutate(area = zoo::rollmedian(area, k = m_kkk1_test, fill = NA, align = "left"),
                 area = mean(area, na.rm = TRUE),
                 class = "Pseudo_QC"
          ) %>%
          ungroup() %>%
          select(area, index, class, batch, compound) %>%
          distinct() %>%
          drop_na(area)
      }
      np1 = bind_rows(np1)
      m_kkk2 = m_kkk
      pseudo_qc1 = pseudo_qc1 %>%
        filter(!name%in%np$name) %>%
        mutate(area = zoo::rollmedian(area, k = m_kkk2, fill = NA, align = "left"),
               area = mean(area, na.rm = TRUE),
               class = "Pseudo_QC"
        ) %>%
        ungroup() %>%
        select(area, index, class, batch, compound) %>%
        distinct() %>%
        drop_na(area) %>%
        bind_rows(., np1)
    }
    tmp_pseudo_qc1 = pseudo_qc1 %>%
      mutate(index = index+m_ind+0.05)
    pseudo_qc1 = pseudo_qc1 %>%
      mutate(index = index+m_ind+0.10) %>%
      bind_rows(., tmp_pseudo_qc1) %>%
      arrange(batch, index) %>%
      mutate(name = paste0("Pseudo_QC", 1:n()))
    n_qc = pseudo_qc1 %>%
      filter(batch%in%train.batch)
    if(nrow(n_qc)<4){
      tmp_pseudo_qc01 = tmp_pseudo_qc1 %>% mutate(index = index+0.05)
      tmp_pseudo_qc001 = tmp_pseudo_qc1 %>% mutate(index = index+0.10)
      pseudo_qc1 = pseudo_qc1 %>%
        bind_rows(., tmp_pseudo_qc01) %>%
        bind_rows(., tmp_pseudo_qc001) %>%
        arrange(index) %>%
        mutate(name = paste0("Pseudo_QC", 1:n()))
    }
    df_psuedoQC = df %>%
      bind_rows(pseudo_qc1) %>%
      arrange(batch, index) %>%
      mutate(area = ifelse(is.na(area), area_og, area)) %>%
      select(-c(area_og, area_tmp, my_mad, high_thresh, low_thresh)) %>%
      ungroup()
    df_psuedoQC_corrected = m_qcrsc(df_psuedoQC, "Pseudo_QC", "yes")
    df_psuedoQC = df_psuedoQC
    df = df %>%
      mutate(area = area_og) %>%
      select(-c(area_og, area_tmp, my_mad, high_thresh, low_thresh)) %>%
      ungroup()

    dd <- c(" ", " -------------- ", " Done! ", " --------------")
    cat(dd, sep = "\n")
    return(list(df = df, df_pseudoQC = df_psuedoQC, df_pseudoQC_corrected = df_psuedoQC_corrected, criteria_table = all_min))
  }

  tt = lapply(comps_within, fun_within)

  df_pseudoQC = list()
  df_pseudoQC_corrected = list()
  criteria_table = list()
  for (i in seq_along(tt)) {
    df_pseudoQC[[i]] = tt[[i]]$df_pseudoQC
    df_pseudoQC_corrected[[i]] = tt[[i]]$df_pseudoQC_corrected
    criteria_table[[i]] = tt[[i]]$criteria_table
  }
  df_pseudoQC = bind_rows(df_pseudoQC) %>%
    arrange(batch, index) %>%
    select(-c(index))
  df_pseudoQC_corrected = bind_rows(df_pseudoQC_corrected) %>%
  arrange(batch, index) %>%
    select(-c(index))
  criteria_table = bind_rows(criteria_table)

  return(list(df = df, df_pseudoQC = df_pseudoQC, df_pseudoQC_corrected = df_pseudoQC_corrected, criteria_table = criteria_table))
}
