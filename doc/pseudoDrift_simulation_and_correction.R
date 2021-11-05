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
#  sdf2index(sdf_input = system.file("extdata", "test.sdf", package = "pseudoDrift"), out_name = "Index.xls") # Can run, but not necessary for this vignette

## ---- eval=TRUE---------------------------------------------------------------
sim1 = simulate_data(compound_names = "tricin",
                     xls_file_name = system.file("extdata", "Index.xls", package = "pseudoDrift"),
                     valid_sdf_file = system.file("extdata", "valid-test.sdf", package = "pseudoDrift"),
                     nbatch = 3,
                     nsamps_per_batch = c(25, 50, 100),
                     QC_freq = c(5,15,20),
                     seed = 123)

## ---- eval=TRUE---------------------------------------------------------------
names(sim1)

## ---- eval=TRUE---------------------------------------------------------------
sim1_sub = sim1 %>% 
  filter(id=="FIO00738")

sim_dat = sim1_sub$sim_mat[[1]]         # simulated no effects
mono = sim1_sub$t1_sim_mat[[1]]         # monotonic
b2b = sim1_sub$t2_sim_mat[[1]]          # batch-to-batch
rando = sim1_sub$t3_sim_mat[[1]]        # random
mono_b2b = sim1_sub$t4_sim_mat[[1]]     # monotonic and batch-to-batch

## ---- fig.height=5, fig.width=6.5, eval=TRUE----------------------------------
# small plotting function
plt_fun = function(x,dat_lab){
  qc_samps = x %>% filter(sample=="QC")
  plt = ggplot(x, aes(batch_index, area))+
    geom_point()+
  geom_point(data = qc_samps, aes(color=sample))+
  geom_line(data = qc_samps, aes(color=sample))+
  facet_grid(cols = vars(batch), scales = "free")+
  labs(title = paste0(dat_lab))
  return(plt)
}

plt_fun(sim_dat, "Simulated")
plt_fun(mono, "Monotonic")
plt_fun(b2b, "Batch-to-batch")
plt_fun(rando, "Random")
plt_fun(mono_b2b, "Monotonic and batch-to-batch")


## ---- eval=TRUE---------------------------------------------------------------
train.batch = "B3"         # Batch to be used as training set
test.breaks = seq(10,12,1) # Test breaking up the training batch into 10 to 12 equally sized sections
test.window = seq(1,9,2)  # Test taking the median every 1 (reduces to the mean per test.break), to 9 samples. This sequence of numbers should be all odd.
test.index = seq(0,5,5)   # Test offsetting the injection order of pseudoQC samples estimated
criteria = "MSE"           # This can be one of: "RSD" (relative standard deviation), "MSE" (mean squared error), or "TSS" (total sum of squares). This is the criteria to be minimized.
n.cores = 10               # Number of cores to use if your machine has the ability to multi-thread processes.


sdc_out = pseudo_sdc(df = mono_b2b,
                     train.batch = train.batch,
                     test.breaks = test.breaks,
                     test.window = test.window,
                     test.index = test.index,
                     criteria = criteria,
                     qc.label = "QC",
                     n.cores = n.cores,
                     min.qc = min(test.breaks))


## ---- eval=TRUE, echo = FALSE, fig.height=8.5, fig.width=6.5------------------
# slightly larger plotting function
plt_fun1 = function(x,train.batch){
  metab = paste0(unique(x$df$compound),"_SIMULATED")
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

## ----eval=TRUE, fig.height=8.5, fig.width=6.5---------------------------------
cor_dat = sdc_out$df_pseudoQC_corrected %>% 
  select(name:area_corrected, class) %>% 
  left_join(., sim_dat) %>% 
  drop_na(area)

cor_dat %>% 
  ggplot(., aes(area_corrected, area, color = batch))+
  geom_point()+
  labs(x = "area_corrected_with_pseudoQC",
       y = "original_area_simulated_without_any_effects")

  
library(caret)
fc = trainControl(method = "cv", number = 10)
fit = train(area ~ area_corrected, 
             data = cor_dat, 
             method = "lm", 
             trControl = fc)

plot(fit$finalModel)

## ----eval=TRUE----------------------------------------------------------------
sessionInfo()

