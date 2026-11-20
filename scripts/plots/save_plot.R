# script for saving plots
theme_fct <- function(font_size1=11,  font_size2=7.5){
  theme_set(theme_bw() +
              theme(
                plot.background = element_blank(),
                panel.background = element_blank(),
                legend.background = element_blank(),
                strip.background = element_rect(fill = "white"),
                plot.caption=element_text(size=font_size2, hjust=0,
                                          margin=margin(t=15)),
                text = element_text(size = font_size1),
                axis.ticks = element_blank(),
                axis.text = element_text(size = font_size1),
                panel.grid.major = element_line(linewidth = 0.25)
              )
            #+
            # theme_cowplot(font_size = 11)
  )
}



save_myplot <- function(plt, plt_nm,
                        width, height,
                        width_pdf = 50, height_pdf = 50) {

  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }

  ggsave(plt_nm,
         egg::set_panel_size(p = plt,
                             width = unit(width, "in"),
                             height = unit(height, "in")),
         width = width_pdf, height = height_pdf,
         limitsize = FALSE, units = c("in"))

  knitr::plot_crop(plt_nm)
}
