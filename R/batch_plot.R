#' @title batch_plot
#' @name batch_plot
#' @description plot sample values along batch
#' Plotting function
#' @param df dataframe containing the data

batch_plot = function(df, qc.label = "QC", xl, yl, tit){
  df <- df %>%
    mutate(class = ifelse(grepl(paste0(qc.label),sample), "QC", "Sample"))
  my_pools <- df %>%
    filter(class=="QC")
  b <- df %>%
    ggplot(., aes(x = batch_index, y = area, col = batch, shape = class)) +
    geom_point(size = 0.5)+
    scale_shape_manual(values=c(19, 19))+
    geom_point(data = my_pools, size = 3.0) +
    geom_line(data = my_pools, aes(x = batch_index, y = area))+
    scale_x_continuous(expand=expansion(mult=c(0.01, 0.03)))+
    facet_grid(cols = vars(batch), scales = "free", space = "free")+
    theme_classic()+
    theme(
      axis.text.x = element_text(vjust = 0.5, hjust=1, angle = 90, size = 12,face = "bold", color = "black"),
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),
      axis.title.y = element_text(size = 16, color = "black", face = "bold"),
      axis.title.x = element_text(size = 16, color = "black", face = "bold"),
      strip.text.x = element_text(size = 14, color = "black", face = "bold"),
      plot.title = element_text(hjust = 0.5,color = "black", face = "bold"),
      panel.spacing.x = unit(0.9, "lines"),
      legend.position="top",
      legend.box="vertical",
      legend.margin=margin(),
      legend.justification='center',
      legend.title = element_blank(),
      legend.direction='horizontal'
    )+
    guides(colour = guide_legend(nrow = 1,override.aes = list(size=5)),
           shape = guide_legend(nrow = 1,override.aes = list(size=c(5,2))))+
    labs(x = paste0(xl), y = paste0(yl),title = paste0(tit))
  return(b)
}
