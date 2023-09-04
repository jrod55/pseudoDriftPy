#' @title simulate_data
#' @name simulate_data
#' @description Simulate LC-MS/MS data or any other spectra from sdf db files.
#' @param db_ids \code{character()} ID(s) of compounds to be read from a validated sdf file.
#' @param compound_names \code{character()} Name of the compound(s) to be queried. Not case sensitive, but will only search compounds which begin with the string entered.
#' @param xls_file_name \code{character(1)} Index file in .xls format. \cr
#' This file can be produced using \link[pseudoDrift]{sdf2Index}.
#' @param valid_sdf_file \code{character(1)} sdf file to be read. \cr
#' This file can also be produced using \link[pseudoDrift]{sdf2Index}.
#' @param nbatch \code{numeric()} Number of batches to simulate.
#' @param nsamps_per_batch \code{numeric()} A numeric vector the same length as 1:nbatch. Number of samples included in each batch to be simulated.
#' @param QC_freq \code{numeric()} A numeric vector the same length as 1:nbatch. Frequency of QC samples.
#' @param multiplyer \code{numeric(1)} value to multiply exact m/z value from MoNA database to have an area measure (default is 1e2)
#' @param sim_sd \code{numeric()} Standard deviation of simulated data. Default value is 0.25*m/z value after multiplier.
#' @param m_eff \code{numeric()} monotonic effect (slope). Should be >=1.
#' @param b_eff_pos \code{numeric()} batch effect in the positive direction (between batch mean differences). Should be <=1.
#' @param b_eff_neg \code{numeric()} batch effect in the negative direction (between batch mean differences). Should be >=1
#' @param seed \code{numeric()} the seed to be used for reproducibility..
#' @param save_rds \code{logical(1)} To write to an rds file.
#' @param rds_name \code{character(1)} Name of the .rds file if save_rds is set to TRUE (default is FALSE).
#' @return A nested tibble. The first column is the db id of the compounds simulated. The second column is the simulated matrix. The following four columns represent four batch simulation scenarios, which are: \cr
#' \itemize{
#' \item 1) t1_sim_mat - a monotonic (up/down effect) by batch.
#' \item 2) t2_sim_mat - a batch to batch block effect.
#' \item 3) t3_sim_mat - random, no systematic change.
#' \item 4) t4_sim_mat - monotonic and batch to batch block effect.
#' }
#' @import dplyr vroom
#' @export
#' @examples
#' sim1 = simulate_data(compound_names = "tricin",
#'                      xls_file_name = system.file("extdata", "Index.xls", package = "pseudoDrift"),
#'                      valid_sdf_file = system.file("extdata", "valid-test.sdf", package = "pseudoDrift"))
#'
#' sim2 = simulate_data(db_ids = c("FIO00738", "FIO00739","FIO00740"),
#'                      xls_file_name = system.file("extdata", "Index.xls", package = "pseudoDrift"),
#'                      valid_sdf_file = system.file("extdata", "valid-test.sdf", package = "pseudoDrift"))
#'
#' ## Providing name or compound db IDs gives the same result
#' identical(sim1,sim2)
simulate_data <- function(db_ids = NULL,
                          compound_names = NULL,
                          xls_file_name = NULL,
                          valid_sdf_file = NULL,
                          nbatch = 3,
                          nsamps_per_batch = c(100,200, 300),
                          QC_freq = c(25,25,25),
                          multiplyer = 1e2,
                          seed = 123,
                          sim_sd = NULL,
                          m_eff = 1.25,
                          b_eff_pos = 1.25,
                          b_eff_neg = 0.75,
                          save_rds = FALSE,
                          rds_name = NULL){
  xls_file <- xls_file_name
  sdf_file <- valid_sdf_file
  compound_name <- compound_names
  set.seed(seed)
  if (is.null(db_ids)&is.null(compound_name)) {
    stop(paste0("Need to provide either a db_id or compound_name"))
  }
  if (!is.null(db_ids)&!is.null(compound_name)) {
    stop(paste0("Need to provide either a db_id or compound_name, NOT both"))
  }
  if (!grepl("\\.xls$", xls_file)) {
    stop(paste0("Check your output file name: ", out_name), " needs to end in .xls")
  }
  ## Read in the sdf index .xls file
  idx <- read.delim(paste0(xls_file), nrows = 1)
  idx <- names(idx)
  db_dat <- suppressMessages(vroom::vroom(paste0(xls_file), col_names = idx,delim = "\t",skip = 1))
  if(!is.null(compound_name)){
    db_dat <- filter(db_dat, grepl(paste0("^",compound_name, collapse = "|"), NAME, ignore.case = TRUE))
  }
  if(!is.null(db_ids)){
    db_dat <- suppressWarnings(filter(db_dat, ID %in% all_of(db_ids)))
  }
  if (nrow(db_dat)==0) {
    stop(paste0("Check your db_id or compound name. There are no matches in the xls_file index file"))
  }
  if (length(QC_freq)!=length(nsamps_per_batch)) {
    stop(paste0("Check your QC_freq. It should be the same length as "))
  }

  ## create sdf object from indices then pull data matrix
  sdfset = ChemmineR::read.SDFindex(file = paste0(sdf_file), index = data.frame(db_dat[,2:3]))
  db_dat = as_tibble(ChemmineR::datablock2ma(datablocklist = ChemmineR::datablock(x = sdfset))) %>%
    janitor::clean_names()
  mz_sim = db_dat %>%
    mutate(mz = as.numeric(exact_mass)) %>%
    summarise(id = id,
              mz_sim = mz*multiplyer) %>%
    distinct()
  b = sort(rep(1:nbatch, nsamps_per_batch))
  nsamps = length(b)
  m_names = paste0("S",1:nsamps, "_B",b)
  pp = list()
  for (i in seq_along(table(b))) {
    pp[[i]] = tibble::tibble(batch_index = seq(QC_freq[i],table(b)[i],QC_freq[i]),
                             batch = paste0("B", all_of(i)),
                             is_QC = TRUE)
  }
  p_idx = bind_rows(pp)
  c = toupper(janitor::make_clean_names(db_dat$name))
  mz_sim = mz_sim %>%
    mutate(compound = paste0(all_of(c),"_SIMULATED"))

  ## Generate the matrices as nested tibbles by mapping functions
  sim_fun = function(x){
    set.seed(seed)
    if (is.null(sim_sd)) {
      m_sd <- x$mz_sim*0.25
    }else{
      m_sd <- sim_sd
    }
    tibble::tibble(tmp = 1,
                   name = m_names,
                   sample = gsub("_.*","",m_names),
                   batch = gsub(".*_","",m_names),
                   compound = x$compound,
                   area = abs(stats::rnorm(nsamps,
                                           mean = x$mz_sim,
                                           sd = m_sd))
    ) %>%
      group_by(tmp) %>%
      mutate(experiment_index = 1:n()
      ) %>%
      group_by(batch) %>%
      mutate(batch_index = 1:n()) %>%
      ungroup() %>%
      left_join(., p_idx, by = c("batch", "batch_index")) %>%
      mutate(area = ifelse(is.na(is_QC), area, x$mz_sim),
             sample = ifelse(is.na(is_QC),sample, "QC")
      ) %>%
      select(-c(tmp, is_QC))
  }
  ## Functions for each batch mode
  ## Monotonic
  t1 = function(x){
    set.seed(seed)
    x %>%
      group_by(batch) %>%
      mutate(up_down = sample(c(T, F),1),
             up = sort(abs(runif(n(),min = 1, max = m_eff))),
             down = sort(abs(runif(n(),min = 1, max = m_eff)), decreasing = TRUE),
             area = ifelse(up_down, area*up, area*down)
      ) %>%
      select(-c(up_down, up, down)) %>%
      ungroup()
  }
  ## Batch-to-batch
  t2 = function(x){
    set.seed(seed)
    x %>%
      group_by(batch) %>%
      mutate(blk = runif(1,min = b_eff_neg, max = b_eff_pos),
             area = area*blk
      ) %>%
      select(-c(blk)) %>%
      ungroup()
  }
  ## Random
  t3 = function(x){
    set.seed(seed)
    x %>%
      group_by(batch) %>%
      mutate(area = area*abs(stats::runif(n()))
      ) %>%
      ungroup()
  }
  ## Monotonic with Batch-to-batch
  t4 = function(x){
    set.seed(seed)
    mm_eff <- sort(c(m_eff,1))
    bb_eff <- sort(c(b_eff_neg, b_eff_pos))

    x %>%
      group_by(batch) %>%
      mutate(n_tile = ntile(batch_index, sum(sample=="QC")*2)) %>%
      mutate(up_down = sample(c(T, F),1)) %>%
      mutate(up = sort(abs(runif(n(),min = mm_eff[1], max = mm_eff[2]))),
             down = sort(abs(runif(n(),min = mm_eff[1], max = mm_eff[2])), decreasing = TRUE),
             area = ifelse(up_down, area*up, area*down)
      ) %>%
      group_by(batch, n_tile) %>%
      mutate(up_down1 = sample(c(T, F),1)) %>%
      group_by(batch) %>%
      mutate(n_tile1 = rleid(up_down1)) %>%
      group_by(batch, n_tile1) %>%
      mutate(blk = runif(1,min = bb_eff[1], max = bb_eff[2]),
             area = area*blk
      ) %>%
      ungroup() %>%
      select(-c(up_down,up_down1, up, down, n_tile, n_tile1, blk))
  }
  ## Nested tibble to be returned
  db_dat = db_dat %>%
    mutate(name_in_tibble = mz_sim$compound, .before = 1)
  sim_dat = mz_sim %>%
    group_by(id) %>%
    tidyr::nest() %>%
    mutate(sim_mat = purrr::map(data, sim_fun),
           t1_sim_mat = purrr::map(sim_mat, t1),
           t2_sim_mat = purrr::map(sim_mat, t2),
           t3_sim_mat = purrr::map(sim_mat, t3),
           t4_sim_mat = purrr::map(sim_mat, t4)
    ) %>%
    select(-c(data)) %>%
    left_join(.,db_dat, by = "id")
  if (save_rds) {
    saveRDS(sim_dat, "simulated_data.rds")
  }
  return(sim_dat)
}
