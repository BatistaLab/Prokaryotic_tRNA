library(tidyverse)
library(viridis)
library(patchwork)
theme_set(theme_classic())

tRNA_order <- c("transcript_reference")
#fill out a list of transcripts in order that they should appear on a heatmap
x_pos <- c(
  "-1",
  "0",
  "1",
  "2",
  "3",
  "4",
  "5",
  "5a",
  "6",
  "7",
  "8",
  "9",
  "10",
  "11",
  "12",
  "13",
  "14",
  "15",
  "16",
  "17",
  "17a",
  "18",
  "19",
  "20",
  "20a",
  "20b",
  "21",
  "22",
  "23",
  "24",
  "25",
  "26",
  "27",
  "28",
  "29",
  "30",
  "31",
  "32",
  "33",
  "34",
  "35",
  "36",
  "37",
  "38",
  "39",
  "40",
  "41",
  "42",
  "43",
  "44",
  "45",
  "V1",
  "V2",
  "V3",
  "V4",
  "V5",
  "V6",
  "V7",
  "V8",
  "V9",
  "V10",
  "V11",
  "V12",
  "V13",
  "V14",
  "V15",
  "V16",
  "V17",
  "V18",
  "V19",
  "V20",
  "46",
  "47",
  "48",
  "49",
  "50",
  "51",
  "52",
  "53",
  "54",
  "55",
  "56",
  "57",
  "58",
  "59",
  "60",
  "61",
  "62",
  "63",
  "64",
  "65",
  "66",
  "67",
  "67a",
  "68",
  "69",
  "70",
  "71",
  "72",
  "73",
  "74",
  "75",
  "76"
)

input_dataframe <- read_csv("collated_freqs.csv")
#takes *_collated_freqs_CCA_NEW.csv files
input_dataframe <-
  input_dataframe %>%
  mutate(gene_code = str_extract(tRNA_index, pattern = "^[^_]+(?=_)")) %>%
  mutate(read_position = as.numeric(str_extract(
    tRNA_index, pattern = "_(.*)", group = 1
  )))

input_dataframe <- input_dataframe %>%
  arrange(gene_code, read_position) %>%
  group_by(gene_code) %>%
  mutate(
    rep1_misincorp = case_when(
      ref_base == "A" ~ 1 - A_freq_rep1,
      ref_base == "C" ~ 1 - C_freq_rep1,
      ref_base == "G" ~ 1 - G_freq_rep1,
      ref_base == "T" ~ 1 - T_freq_rep1,
      TRUE ~ NA
    )
  ) %>%
  mutate(
    rep2_misincorp = case_when(
      ref_base == "A" ~ 1 - A_freq_rep2,
      ref_base == "C" ~ 1 - C_freq_rep2,
      ref_base == "G" ~ 1 - G_freq_rep2,
      ref_base == "T" ~ 1 - T_freq_rep2,
      TRUE ~ NA
    )
  ) %>%
  mutate(
    rep3_misincorp = case_when(
      ref_base == "A" ~ 1 - A_freq_rep3,
      ref_base == "C" ~ 1 - C_freq_rep3,
      ref_base == "G" ~ 1 - G_freq_rep3,
      ref_base == "T" ~ 1 - T_freq_rep3,
      TRUE ~ NA
    )
  ) %>%
  mutate(
    avg_misincorp = case_when(
      ref_base == "A" ~ (1 - A_avg_freq),
      ref_base == "C" ~ (1 - C_avg_freq),
      ref_base == "G" ~ (1 - G_avg_freq),
      ref_base == "T" ~ (1 - T_avg_freq),
      TRUE ~ NA
    )
#calculating misincorporation frequencies
  ) %>%
  mutate(rep1_termination = ((lead(
    depth_rep1, default = 0
  ) - depth_rep1) / (lead(
    depth_rep1, default = 0
  )))) %>%
  mutate(rep2_termination = ((lead(
    depth_rep2, default = 0
  ) - depth_rep2) / (lead(
    depth_rep2, default = 0
  )))) %>%
  mutate(rep3_termination = ((lead(
    depth_rep3, default = 0
  ) - depth_rep3) / (lead(
    depth_rep3, default = 0
  )))) %>%
#calculating termination scores
  mutate(avg_termination = ((lead(avg_depth, default = 0) - avg_depth) / (lead(avg_depth, default = 0))))  %>%
  mutate(rep1_misincorp = case_when(depth_rep1 < 100 ~ NA,
                                    TRUE ~ rep1_misincorp)) %>%
  mutate(rep2_misincorp = case_when(depth_rep2 < 100 ~ NA,
                                    TRUE ~ rep2_misincorp)) %>%
  mutate(rep3_misincorp = case_when(depth_rep3 < 100 ~ NA,
                                    TRUE ~ rep3_misincorp)) %>%
  mutate(avg_misincorp = case_when(avg_depth < 100 ~ NA,
                                   TRUE ~ avg_misincorp)) %>%
  #filters based on a set minimum read depth
  mutate(
    rep1_termination = case_when(
      depth_rep1 < 100 ~ NA,
      rep1_termination < 0 ~ 0,
      TRUE ~ rep1_termination
    )
  ) %>%
  mutate(
    rep2_termination = case_when(
      depth_rep2 < 100 ~ NA,
      rep2_termination < 0 ~ 0,
      TRUE ~ rep2_termination
    )
  ) %>%
  mutate(
    rep3_termination = case_when(
      depth_rep3 < 100 ~ NA,
      rep3_termination < 0 ~ 0,
      TRUE ~ rep3_termination
    )
  ) %>%
  mutate(avg_termination = case_when(avg_depth < 100 ~ NA,
                                     avg_termination < 0 ~ 0,
                                     TRUE ~ avg_termination))
#filters based on a set minimum read depth

write_csv(input_dataframe, "processed.csv")
#writing a file to quickly generate figures
input <- read_csv("processed.csv")

ggplot(data = input, aes(
  x = factor(input$`Annot. Pos`, levels = x_pos),
  y = fct_rev(factor(tRNA_name, levels = tRNA_order)),
  fill = avg_misincorp
)) +
  geom_tile() +
  theme(aspect.ratio = 100 / 150,
        panel.background = element_rect(fill = "#e4f5e5")) +
  ggtitle("Misincorporation") +
  scale_fill_viridis(option = "mako",
                     direction = -1,
                     limits = c(0, 1.25))

ggsave("misincorporation.pdf",
       dpi = 800,
       scale = 5)
#misincorporation heatmap

ggplot(data = input, aes(
  x = factor(input$`Annot. Pos`, levels = x_pos),
  y = fct_rev(factor(tRNA_name, levels = tRNA_order)),
  fill = avg_termination
)) +
  geom_tile() +
  theme(aspect.ratio = 100 / 150,
        panel.background = element_rect(fill = "#e4f5e5")) +
  ggtitle("Termination") +
  scale_fill_viridis(option = "mako",
                     direction = -1,
                     limits = c(0, 1.25))

ggsave("termination.pdf",
       dpi = 800,
       scale = 5)
#termination heatmap
