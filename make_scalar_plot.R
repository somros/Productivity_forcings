library(tidyverse)
library(PNWColors)


ratios <- readRDS("output/roms_change_lm.RDS")
grps <- read.csv("data/GOA_Groups.csv")

# Create plotting data: two points per row (2020 = 1, 2100 = relative_change)
plot_data <- ratios %>%
  ungroup() %>%
  select(Code, simulation, relative_change) %>%
  mutate(
    `2020` = 1,
    `2100` = relative_change
  ) %>%
  pivot_longer(cols = c(`2020`, `2100`), names_to = "year", values_to = "value") %>%
  mutate(year = as.integer(year)) %>%
  left_join(grps %>% select(Code, LongName))

# Plot
all_cap_col <- pnw_palette(name = "Sunset2", n = 4, type = "discrete")
cap_col     <- all_cap_col[2:4]

p <- ggplot(plot_data, aes(x = year, y = value, color = simulation, group = simulation)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  facet_wrap(~ LongName, nrow = 1) +
  scale_x_continuous(breaks = seq(2020, 2100, 10)) +
  scale_y_continuous(limits = c(0,1))+
  scale_color_manual(
    values = cap_col,
    labels = c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5", ssp585 = "SSP5-8.5")
  ) +
  labs(
    x = NULL,
    y = "Growth rate scalar (2020 = 1)",
    color = "Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p

ggsave("scalar.png", p, dpi = 600, width = 8, height = 4)
