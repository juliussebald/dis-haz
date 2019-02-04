# Load packages and set options -------------------------------------------

{
  library(tidyverse) # version 1.2.1
  library(patchwork)  # version 0.0.1
}

# Load data ---------------------------------------------------------------

estimates <- list.files("../results/", glob2rx("estimates*.csv"), recursive = TRUE, full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows() %>%
  mutate(name = factor(name, levels = c("Artificial", "Elevation","Area", 
                                        "Elevation ratio", "Circularity", "Melton ratio", "Elongation",    
                                        "Extent x Type", "Extent", "Type", "Patch density", "Forest")))

p_estimate <- list(a = estimates %>% 
                     mutate(type = case_when(name %in% c("Forest", "Patch density") ~ "Forest",
                                             name %in% c("Extent x Type", "Extent", "Type") ~ "Disturbance",
                                             TRUE ~ type)) %>%      
                     mutate(type = factor(type, levels = c("General", "Geomorphology", "Forest", "Disturbance"))) %>% 
                     mutate(model = factor(model, labels = c("Occurrence", "Frequency"))) %>% 
       split(.$type),
     b = list(NULL, NULL, NULL, "Effect size"),
     c = list(element_blank(), element_blank(), element_blank(), element_text()),
     d = list(element_text(size = 10), element_blank(), element_blank(), element_blank()),
     e = list(element_blank(), element_blank(), element_blank(), element_line()),
     f = list("right", "none", "none", "none")) %>%
  pmap(.l = ., .f = function (a, b, c, d, e, f) {
    ggplot(a, aes(x = name, y = value)) +
      geom_violin(aes(fill = model)) +
      theme_bw() +
      theme(legend.position = f,
            legend.title = element_text(size = 9),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = d,
            strip.text.y = element_text(size = 10),
            axis.ticks.x = e,
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.text.x = c,
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed", col = scales::muted("red")) +
      labs(y = b, x = NULL, fill = "Model") +
      facet_grid(type~process) +
      ylim(-0.9, 0.9) +
      scale_fill_brewer(palette = "Greys") +
      guides(fill = guide_legend(ncol = 1, 
                                 keywidth = 0.15,
                                 keyheight = 0.1,
                                 default.unit = "inch"))
  }) %>%
  wrap_plots(ncol = 1, heights = c(1, 1.2, 0.6, 1))

ggsave("../results/estimates_combined.pdf", p_estimate, width = 5.5, height = 4.5)
