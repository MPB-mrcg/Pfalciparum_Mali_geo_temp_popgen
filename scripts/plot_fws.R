
rm(list = ls())
library(tidyverse)

fws <- read_delim("Fws.txt", delim = "\t", col_names = c("Samples", "Fws", "Locations", "COI"), skip = 1)

Fws <- fws %>% 
    group_by(Locations) %>% 
    arrange(Fws) %>% 
    mutate(Cum = row_number()/n()) %>% 
    arrange(Locations)

Fws %>% 
    ggplot(aes(x = Cum, y = Fws, color = Locations)) +
    geom_point() + 
    geom_hline(yintercept = 0.95, linetype="dashed", color = "red") +
    labs(x="Proportion", y= "Fws")
