library(grf)
library(engression)
library(evd)
library(rgl)
library(htmlwidgets)
library(mvtnorm)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
#library(progression)


unif <- function (x, na.rm = FALSE){
    if (na.rm) {
        na.last <- NA
    }
    else {
        na.last <- "keep"
    }
    rank(x, ties.method = "first", na.last = na.last)/(sum(!is.na(x)) +
        1)
}




lst_methods <- list("RF" = "RF",
                    "LLF" = "LLF",
                    "engression" = "engression",
                    "progression" = "progression"
)

my_palette <- list(
  "red" = "#D55E00",
  "blue" = "#0072B2",
  "green" = "#009E73",
  "yellow" = "#E69F00",
  "pink" = "#CC79A7",
  "light_blue" = "#56B4E9",
  "grey" = "#999999",
  "background" = "#332288"
)



my_palette_methods <- list(
  c("reg_method" = lst_methods$RF, "color" = my_palette$red, "fill" = "white"),
  c("reg_method" = lst_methods$LLF, "color" = my_palette$blue, "fill" = "white"),
  c("reg_method" = lst_methods$engression, "color" = my_palette$green, "fill" = "white"),
  c("reg_method" = lst_methods$progression, "color" = my_palette$yellow, "fill" = "white")
) %>%
  purrr::transpose() %>%
  as_tibble() %>%
  unnest(cols = c(reg_method, color, fill))

my_col <-  my_palette_methods %>%
  dplyr::select(reg_method, color) %>%
  deframe()

my_fill <-  my_palette_methods %>%
  dplyr::select(reg_method, fill) %>%
  deframe()





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

theme_fct()



save_myplot <- function(plt, plt_nm,
                        width, height,
                        width_pdf = 50, height_pdf = 50,
                        crop = TRUE, cairo = TRUE) {

  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }

  if (cairo) {
    ggsave(plt_nm,
           egg::set_panel_size(p = plt,
                               width = unit(width, "in"),
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"),
           device = cairo_pdf, family = "Arial")
  } else {
    ggsave(plt_nm,
           egg::set_panel_size(p = plt,
                               width = unit(width, "in"),
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"))
  }

  if (crop){
    knitr::plot_crop(plt_nm)
  }
}


texify_column <- function(column, letter) {
  ## vector character -> factor
  ## paste latex formula of the form "$letter = column$"
  factor(column,
         levels = unique(column),
         labels = TeX(paste0(
           "$", letter, " = ",
           unique(column), "$"
         ))
  )
}

