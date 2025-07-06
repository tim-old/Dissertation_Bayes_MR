# plot template
plot_template <- function(data)  {
  data %>%
    ggplot() +
    ## set generic theme
    scale_fill_edin_d("digital") +
    scale_colour_edin_d("digital") +
    scale_fill_edin_c("blue_ramp") +
    scale_colour_edin_c("blue_ramp") +
    theme_bw() +
    theme(
      plot.title = element_text(colour = edin_spruce_grey_hex,
                                size = 14),
      plot.subtitle = element_text(colour = edin_spruce_grey_hex,
                                   size = 10),
      plot.caption = element_text(face = "italic",
                                  colour = edin_spruce_grey_hex,
                                  size = 10),
      axis.title = element_text(colour = edin_spruce_grey_hex,
                                size = 12),
      axis.text = element_text(colour = edin_spruce_grey_hex,
                               size = 10),
      panel.grid = element_line(colour = edin_light_blue_hex),
      strip.text = element_text(colour = edin_spruce_grey_hex,
                                size = 10),
      strip.background = element_rect(fill = edin_light_blue_hex,
                                      colour = edin_spruce_grey_hex)
    )
}
