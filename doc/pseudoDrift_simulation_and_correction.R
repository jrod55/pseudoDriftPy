## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=TRUE, message=FALSE------------------------------------------------
# install.packages("remotes")
# remotes::install_github("jrod55/pseudoDrift")

library(pseudoDrift)
library(tidyverse)
library(data.table)

## ---- eval=FALSE--------------------------------------------------------------
#  # sdf2index(sdf_input = "/path_to_your/desired.sdf", out_name = "Index.xls") # Example
#  # sdf2index(sdf_input = system.file("extdata", "test.sdf", package = "pseudoDrift"), out_name = "Index.xls") # Can run, but not necessary for this vignette

## ---- eval=TRUE---------------------------------------------------------------
sim1 = simulate_data(compound_names = "tricin",
                     xls_file_name = system.file("extdata", "Index.xls", package = "pseudoDrift"),
                     valid_sdf_file = system.file("extdata", "valid-test.sdf", package = "pseudoDrift"),
                     # nbatch = 3,
                     # nsamps_per_batch = c(100, 200, 300),
                     # QC_freq = c(25, 25, 25),
                     # seed = 123
                     )

## ---- eval=TRUE---------------------------------------------------------------
names(sim1)

## ---- eval=TRUE---------------------------------------------------------------
sim1_sub = sim1 %>% 
  filter(id=="FIO00738")

sim_dat = sim1_sub$sim_mat[[1]]         # Simulated no effects
mono = sim1_sub$t1_sim_mat[[1]]         # Monotonic
b2b = sim1_sub$t2_sim_mat[[1]]          # Batch-to-batch
rando = sim1_sub$t3_sim_mat[[1]]        # Random
mono_b2b = sim1_sub$t4_sim_mat[[1]]     # Monotonic and batch-to-batch

## ---- eval=TRUE---------------------------------------------------------------
# Small plotting function
plt_fun = function(x,dat_lab){
  qc_samps = x %>% filter(sample=="QC")
  plt = ggplot(x, aes(batch_index, area))+
    geom_point(size = 0.75)+
    geom_point(size = 2, data = qc_samps, aes(color=sample))+
    geom_line(data = qc_samps, aes(color=sample))+
    facet_grid(cols = vars(batch), scales = "free")+
    labs(title = paste0(dat_lab), y = "area_simulated")
  return(plt)
}

## ----fig.align = "center", eval=TRUE------------------------------------------
plt_fun(sim_dat, "Simulated")
plt_fun(mono, "Monotonic")
plt_fun(b2b, "Batch-to-batch")
plt_fun(rando, "Random")
plt_fun(mono_b2b, "Monotonic and batch-to-batch")

## ---- eval=TRUE, fig.align = "center", fig.height=6, fig.width=10.5-----------
# Assign replicates to simulated data. Note no technical replicates are included here, assuming only biological replicates.
n_reps = 3
qcs = mono_b2b %>%
  filter(sample=="QC")
tmp = mono_b2b %>%
  filter(!sample=="QC") %>%
  mutate(n_samps = n()/all_of(n_reps))
dat_rep = rep(1:n_reps, unique(tmp$n_samps))
dat_sam = paste0("S",rep(1:unique(tmp$n_samps), each = n_reps))

# The Monotonic + batch-to-batch simulated data with replicateds assigned 
mono_b2b_WR = tmp %>%
  mutate(rep = all_of(dat_rep),
         sample = all_of(dat_sam),
         rep_tech = 1) %>%
  bind_rows(.,qcs) %>%
  arrange(experiment_index)

pw_out = pw_outlier(df = mono_b2b_WR, 
                    return_plot = TRUE, 
                    samps_exclude = "QC")

pw_out$batch_plots

## ---- eval=TRUE, fig.align = "center"-----------------------------------------
# Make a plot of the data removed by pairwise outlier removal
df_rm = pw_out$df_rm %>%
  mutate(sample = "Outlier")
plt_fun(pw_out$df_cleaned, "Monotonic and batch-to-batch pairwise cleaned")+
  geom_point(data = df_rm, size = 1.25, aes(fill=sample), color = "blue")

## ---- eval=TRUE---------------------------------------------------------------
# Use the pw cleaned data
mono_b2b_cleaned = pw_out$df_cleaned

# A good starting point for parameters to test
train.batch = "B3"
df_param = mono_b2b_cleaned %>%
  filter(!sample=="QC") %>%
  group_by(batch) %>%
  summarise(samps = n_distinct(name))

# Sample size of smallest batch
s_perBatch = min(df_param$samps)

# Test breaks in smallest batch, using the number of QC samples simulated there
mnqc = 25
t.b = round(s_perBatch/mnqc)+1

# Proportional in larger training set
df_param = df_param %>%
  mutate(bre = round(samps*all_of(t.b)/all_of(s_perBatch)))
t.b.train = df_param %>%
  filter(batch == all_of(train.batch)) %>%
  pull(bre)
t.b.train = seq(t.b.train-3, t.b.train+3, 1)

# Test window for median smoothing. Should be a vector of odd numbers
w.n.max = round(min(s_perBatch/t.b))/2
w.n = 2*floor(w.n.max/2)+1
if (w.n>w.n.max) {
  w.n = w.n-2
}
w.n = seq(1,w.n,2)

# Test index offset. Larger values can be tested with larger datasets
ti.max = round(s_perBatch*0.15)
t.i = seq(0,ti.max,1)

# test.breaks = t.b.train     # Test breaking up the training batch into 10 to 12 equally sized sections
# test.window = w.n           # Test taking the median every 1 (reduces to the mean per test.break), to 9 samples. This sequence of numbers should be all odd.
# test.index = t.i            # Test offsetting the injection order of pseudoQC samples estimated
# criteria = "RSD"            # This can be one of: "RSD" (relative standard deviation), "MSE" (mean squared error), or "TSS" (total sum of squares). This is the criteria to be minimized.
# n.cores = 15                # Number of cores to use if your machine has the ability to multi-thread processes.
# quantile_increment = 0.05   # Quantile values (above/below increment) of the data surrounding QCs to use to estimate pseudoQCs.

# To reduce computing time, we are going to use values previously obtained by running the exhaustive set of tests with the inputs commented out above.

# For MSE and TSS criteria:
# test.breaks = 12 ; test.window = 7; test.index = 5; criteria = "MSE"; quantile_increment = 0.05; n.cores = 1

# For RSD:
test.breaks = 13 ; test.window = 9; test.index = 14; criteria = "RSD"; quantile_increment = 0.05; n.cores = 1

sdc_out = pseudo_sdc(df = mono_b2b_cleaned,
                     n.cores = n.cores,
                     train.batch = train.batch,
                     test.breaks = test.breaks,
                     test.window = test.window,
                     test.index = test.index,
                     criteria = criteria,
                     qc.label = "QC",
                     min.qc = min(test.breaks),
                     quantile_increment = quantile_increment)


## ---- eval=TRUE, fig.height=7, fig.width=8.5, fig.align = "center"------------
# slightly larger plotting function
plt_fun1 = function(x,train.batch){
  metab = unique(x$df$compound)
  x1 = x$df_pseudoQC %>% 
    group_by(batch) %>% 
    mutate(index = 1:n(),
           class = factor(class, levels = c("QC", "Pseudo_QC", "Sample")))
  
  qcs1 = x1 %>% 
    filter(class%in%c("QC","Pseudo_QC")) %>% 
    mutate(sample = class)
  
  plt1 = ggplot(x1, aes(index, area))+
    geom_point()+
    geom_point(data = qcs1, aes(color=sample))+
    geom_line(data = qcs1, aes(color=sample))+
    facet_grid(cols = vars(batch), scales = "free")+
    labs(title = paste0(metab," raw data using ",train.batch, " data to train pseudo-QC"))
  
  legend = cowplot::get_legend(plt1)
  
  x2 = x$df_pseudoQC_corrected %>% 
    group_by(batch) %>% 
    mutate(index = 1:n(),
           class = factor(class, levels = c("QC", "Pseudo_QC", "Sample")))
  
  qcs2 = x2 %>% 
    filter(class%in%c("QC","Pseudo_QC")) %>% 
    mutate(sample = class)
  
  plt2 = ggplot(x2, aes(index, area_corrected))+
    geom_point()+
    geom_point(data = qcs2, aes(color=sample))+
    geom_line(data = qcs2, aes(color=sample))+
    facet_grid(cols = vars(batch), scales = "free")+
    labs(title = paste0(metab," pseudo-QC corrected data"))
  
  cplt = cowplot::plot_grid(plt1, plt2, ncol = 1, labels = c("A","B"))
  return(cplt)
}

plt_fun1(sdc_out,train.batch)

## ----eval=TRUE,  fig.align = "center"-----------------------------------------
cor_dat = sdc_out$df_pseudoQC_corrected %>%
  select(name,area_corrected, class) %>%
  left_join(., sim_dat) %>%
  drop_na(area)

cor_dat %>% 
  ggplot(., aes(area_corrected, area, color = batch))+
  geom_point(alpha=0.5)+
  labs(x = "area_corrected_with_pseudoQC",
       y = "original_area_simulated_without_any_effects")


## ----eval=TRUE, fig.height=10, fig.width=10, fig.align = "center"-------------
library(caret)
fc = trainControl(method = "cv", number = 10)
fit = train(area ~ area_corrected, 
             data = cor_dat, 
             method = "lm", 
             trControl = fc)
fit
par(mfrow = c(2, 2))
plot(fit$finalModel)

## ----eval=TRUE, fig.height=10, fig.width=10, fig.align = "center"-------------
## A short function to format the data how pmp expects it
m_qcrsc = function(x, y){
  x = x %>% 
    mutate(class = ifelse(sample%in%all_of(y), "QC", "Sample"))
  t_meta = colnames(x)
  t_meta = t_meta[!t_meta%in%c("name", "compound")]
  tqc = x %>%
    pivot_wider(id_cols = !all_of(t_meta), names_from = name, values_from = area)
  t_rn = tqc$compound
  tqc = as.matrix(tqc[,-1])
  rownames(tqc) = t_rn
  mc = pmp::QCRSC(df=tqc,
                  minQC = 4,
                  order=seq_along(x$experiment_index),
                  batch=x$batch,
                  classes=x$class,
                  qc_label = paste0(y))
  z = x %>%
    mutate(area_corrected = mc[1,], .before = area)
  return(z)
  }

## Run the function then conduct the same diagnostics tests done with pseudoQC above. 
true_QC = m_qcrsc(mono_b2b_cleaned, "QC")

## ----eval=TRUE,  fig.align = "center"-----------------------------------------
cor_dat_true_QC = true_QC %>%
  select(name,area_corrected, class) %>%
  left_join(., sim_dat) %>%
  drop_na(area)

cor_dat_true_QC %>% 
  ggplot(., aes(area_corrected, area, color = batch))+
  geom_point(alpha=0.5)+
  labs(x = "area_corrected_with_true_QC",
       y = "original_area_simulated_without_any_effects")


## ----eval=TRUE, fig.height=10, fig.width=10, fig.align = "center"-------------
fit_true_QC = train(area ~ area_corrected, 
             data = cor_dat_true_QC, 
             method = "lm", 
             trControl = fc)
fit_true_QC
par(mfrow = c(2, 2))
plot(fit_true_QC$finalModel)

## ----eval=TRUE----------------------------------------------------------------
sessionInfo()

