library(tidyverse)
library(reshape2)
library(umap)
library(patchwork)
library(scales)

file_fly_biq_dump <- "./kiq_dump_fly.tsv.gz"
file_fly_exps <- "./metadata_fly.tsv"
file_encode_biq_dump <- "kiq_dump_ENCODE.tsv.gz"
file_encode_exps <- "./metadata_ENCODE.tsv"

# ------------------------------------------------------------------------
# FLY
# ------------------------------------------------------------------------

fly.exps <- read_tsv(file_fly_exps)
names(fly.exps) <- make.names(names(fly.exps))
fly.exps <- rename(fly.exps, ID = SRA.id)

fly.dump <- read_tsv(file_fly_biq_dump, col_names = c("kmer", "ID", "count", "rpm"))
# filter out kmers that occur in cDNAs:
linear_kmers <- fly.dump %>% filter(ID == "_LINEAR") %>% pull(kmer)
fly.dump <- fly.dump %>% filter(!kmer %in% linear_kmers)

# Number of k-mers also found in linear transcripts:
length(linear_kmers)
# Number of counted k-mers:
fly.dump %>% summarize(totalcount = sum(count)) %>% pull(totalcount)
# Number of unique k-mers found across all samples:
fly.dump %>% pull(kmer) %>% sort %>% unique %>% length

# ------------------------------------------------------------------------
# Figure 2a: UMAP dimension reduction

m <- acast(fly.dump, ID ~ kmer, value.var = "rpm", fill = 0.0)
set.seed(12222213)
umap <- umap(m)
df.umap <- as.data.frame(umap$layout)
colnames(df.umap) <- c("Dim1", "Dim2")
df.umap$ID <- rownames(df.umap)
df.join <- df.umap %>%
  inner_join(fly.exps, by = "ID") %>%
  mutate(cell_line = ifelse(str_detect(group, "BG3|CME|S1|S2|Kc167"), TRUE, FALSE)) %>%
  mutate(
    coarse_group = case_when(
      str_detect(group, "whole embryo") ~ "whole embryo",
      str_detect(group, "CNS") ~ "larvae / pupae CNS",
      str_detect(group, "head") ~ "head",
      str_detect(group, "ovaries|testes") ~ "ovaries / testes",
      str_detect(group, "digestive") ~ "digestive",
      str_detect(group, "gut") ~ "gut",
      str_detect(group, "carcass") ~ "carcass",
      cell_line == TRUE ~ "cell line",
      TRUE ~ "other"
    )
  )
blues <- brewer_pal(palette = "RdYlBu")(11)[9:11]
acc <- brewer_pal(palette = "Accent")(8)
plot.cols <- c(
  "whole embryo" = blues[3], "larvae / pupae CNS" = blues[2], "head" = blues[1], "other" = acc[8],
  "ovaries / testes" = acc[1], "digestive" = acc[3], "gut" = acc[7], "cell line" = acc[6], "carcass" = acc[2]
)
lvls <- c("whole embryo", "larvae / pupae CNS", "head", "carcass", "digestive", "gut", "ovaries / testes", "other", "cell line")
p.umap <- ggplot(df.join) +
  geom_point(aes(x = Dim1, y = Dim2, color = fct_relevel(coarse_group, lvls))) +
  xlab("") +
  ylab("") +
  theme_bw() +
  scale_x_continuous(position = "top") +
  scale_color_manual(values = plot.cols) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  theme(
    legend.justification = c(1, 0),
    legend.position = c(1, 0),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.key.height = unit(11, "pt")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = rel(0.8)),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1)) +
  theme(aspect.ratio = 1)

# ------------------------------------------------------------------------
# Figure 2b: RPM plot over time, as in Westholm2014 Figure 7a

timepoints <- fly.exps %>%
  filter(str_detect(description, "after egg laying")) %>%
  mutate(time = ordered(str_replace(description, "embryos, (.+) hr after egg laying", "\\1"))) %>%
  separate(time, c("starttime", "endtime"), sep = "-", remove = FALSE) %>%
  mutate(endtime = as.integer(endtime), time = fct_reorder(time, endtime), time = factor(paste(time, "hr embryo"))) %>%
  pull(time)

df.plot <- fly.dump %>%
  inner_join(fly.exps, by = "ID") %>%
  filter(group == "whole embryo" | str_detect(group, "heads|CNS")) %>%
  mutate(label = case_when(
    str_detect(description, "CNS") ~ paste(tolower(str_extract(description, "(?i)larvae|pupae")), "CNS"),
    str_detect(description, "after egg laying") ~ paste(str_extract(description, "\\d+\\-\\d+"), "hr embryo"),
    str_detect(description, "female.*day") ~ paste(str_extract(description, "\\d+ day"), "F head"),
    str_detect(description, "male.*day") ~ paste(str_extract(description, "\\d+ day"), "M head"),
    TRUE ~ "NA"
  )) %>%
  mutate(label = factor(label, levels = rev(c(as.character(unique(timepoints)), "larvae CNS", "pupae CNS", "1 day F head", "1 day M head", "4 day F head", "4 day M head", "20 day F head", "20 day M head")))) %>%
  group_by(label, ID, raw_reads) %>%
  summarise(n = n()) %>%
  mutate(coarse_group = case_when(
    str_detect(label, "embryo") ~ "whole embryo",
    str_detect(label, "CNS") ~ "larvae / pupae CNS",
    str_detect(label, "day") ~ "head"
  )) %>%
  mutate(jpm = n / raw_reads * 1e6)

p.rpm <- df.plot %>%
  ggplot() +
  geom_point(aes(y = label, x = jpm, color = coarse_group)) +
  scale_colour_manual(values = plot.cols, guide = FALSE) +
  scale_x_continuous(position = "top") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title.x = element_text(size = rel(0.8)), axis.ticks = element_blank()) +
  ylab("") +
  xlab("BSJ k-mers / million reads") +
  theme(aspect.ratio = 1)

# combine Figure 2 a and b and save to PDF
p.fig2 <- p.umap + labs(tag = "a") + p.rpm + labs(tag = "b") + plot_annotation(tag_levels = "a")
ggsave(p.fig2, file = "figure2.pdf", width = 11, height = 5)

# ------------------------------------------------------------------------
# Suppl Figure 1a: cumulative number of unique k-mers across all samples

fly.dump$ID <- reorder(fly.dump$ID, fly.dump$ID, FUN = length)

df.fly.cumulative_kmers <- fly.dump %>%
  arrange(desc(ID)) %>%
  mutate(cum_unique_kmers = cumsum(!duplicated(kmer))) %>%
  group_by(ID) %>%
  summarise(cumcount = last(cum_unique_kmers)) %>%
  select(ID, cumcount) %>%
  mutate(ID = as.integer(ID))

p.fly.cumulative <- df.fly.cumulative_kmers %>%
  ggplot() +
  geom_line(aes(x = rev(ID), y = cumcount)) +
  xlab("SRA files") +
  ylab("Cumulative count of unique BSJ k-mers") +
  theme_bw()

# ------------------------------------------------------------------------
# Suppl Figure 2a: Abundance of BSJ k-mers across all samples

p.suppfig2a <- fly.dump %>%
  group_by(kmer) %>%
  summarise(totalcount = sum(count)) %>%
  group_by(totalcount) %>%
  summarise(kmers_per_totalcount = n_distinct(kmer)) %>%
  ggplot() +
  geom_point(aes(x = totalcount, y = kmers_per_totalcount)) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks() +
  theme_bw() +
  ylab("Number of BSJ k-mers") +
  xlab("Abundance")

# ------------------------------------------------------------------------
# Suppl Figure 3a: Number of samples in which BSJ k-mers are observed

p.suppfig3a <- fly.dump %>%
  group_by(kmer) %>%
  summarise(n_samples = n_distinct(ID)) %>%
  group_by(n_samples) %>%
  summarise(n_kmer = n_distinct(kmer)) %>%
  ggplot() +
  geom_point(aes(y=n_kmer,x=n_samples)) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks() +
  theme_bw() +
  ylab("Number of BSJ k-mers") +
  xlab("Number of samples")

# ------------------------------------------------------------------------
# Suppl Figure 4: plot number of unique k-mers per sample and cumulative number
#                 of unique k-mers for samples from Figure 2b.

df.plot <- fly.dump %>%
  mutate(ID = as.character(ID)) %>%
  inner_join(fly.exps, by = "ID") %>%
  filter(group == "whole embryo" | str_detect(group, "heads|CNS")) %>%
  mutate(label = case_when(
    str_detect(description, "CNS") ~ paste(tolower(str_extract(description, "(?i)larvae|pupae")), "CNS"),
    str_detect(description, "after egg laying") ~ paste(str_extract(description, "\\d+\\-\\d+"), "hr embryo"),
    str_detect(description, "female.*day") ~ paste(str_extract(description, "\\d+ day"), "F head"),
    str_detect(description, "male.*day") ~ paste(str_extract(description, "\\d+ day"), "M head"),
    TRUE ~ "NA"
  )) %>%
  mutate(label = factor(label, levels = rev(c(as.character(unique(timepoints)), "larvae CNS", "pupae CNS", "1 day F head", "1 day M head", "4 day F head", "4 day M head", "20 day F head", "20 day M head")))) %>%
  arrange(desc(label)) %>%
  mutate(cc = cumsum(!duplicated(kmer))) %>%
  group_by(label, ID, raw_reads) %>%
  summarise(ccmax = max(cc), n = n(), u = n_distinct(kmer)) %>%
  mutate(coarse_group = case_when(
    str_detect(label, "embryo") ~ "whole embryo",
    str_detect(label, "CNS") ~ "larvae / pupae CNS",
    str_detect(label, "day") ~ "head"
  ))

p.sum <- df.plot %>%
  mutate(g = "1") %>%
  ggplot() +
  geom_point(aes(y = label, x = u, color = coarse_group)) +
  geom_line(aes(y = as.numeric(label), x = ccmax, color = coarse_group, group = g), lwd = 1, lineend = "round") +
  scale_colour_manual(values = plot.cols, guide = FALSE) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title.x = element_text(size = rel(0.8)), axis.ticks = element_blank()) +
  ylab("") +
  xlab("number of distinct BSJ k-mers and cumulative count") +
  theme(aspect.ratio = 1)
ggsave(p.sum, file = "suppl_figure4.pdf", width = 6, height = 6)

# ------------------------------------------------------------------------
# ENCODE
# ------------------------------------------------------------------------

encode.exps <- read_tsv(file_encode_exps)
names(encode.exps) <- make.names(names(encode.exps))
encode.exps <- encode.exps %>%
  filter(Paired.end == 1) %>% # remove all second rows from paired-end samples
  mutate(ID = str_remove(dbxrefs, "SRA:"))

organs <- read_tsv("./cl_uberon2organs.tsv", quote = "", col_names = c("ID", "organs"))
systems <- read_tsv("./cl_uberon2systems.tsv", quote = "", col_names = c("ID", "systems"))

encode.dump <- read_tsv(file_encode_biq_dump, col_names = c("kmer", "ID", "count", "rpm"))
linear_kmers <- encode.dump %>% filter(ID == "_LINEAR") %>% pull(kmer)
encode.dump <- encode.dump %>% filter(!kmer %in% linear_kmers)

# Number of k-mers also found in linear transcripts:
length(linear_kmers)
# Number of counted k-mers:
encode.dump %>% summarize(totalcount = sum(count)) %>% pull(totalcount)
# Number of unique k-mers found across all samples:
encode.dump %>% pull(kmer) %>% sort %>% unique %>% length

# ------------------------------------------------------------------------
# Suppl Figure 1b: cumulative number of unique k-mers across all samples

encode.dump$ID <- reorder(encode.dump$ID, encode.dump$ID, FUN = length)
df.encode.cumulative_kmers <- encode.dump %>%
  arrange(desc(ID)) %>%
  mutate(cum_unique_kmers = cumsum(!duplicated(kmer))) %>%
  group_by(ID) %>%
  summarise(cumcount = last(cum_unique_kmers)) %>%
  select(ID, cumcount) %>%
  mutate(ID = as.integer(ID))

p.encode.cumulative <- df.encode.cumulative_kmers %>%
  ggplot() +
  geom_line(aes(x = rev(ID), y = cumcount)) +
  xlab("SRA files") +
  ylab("Cumulative count of unique BSJ k-mers") +
  theme_bw()

# ------------------------------------------------------------------------
# combine Suppl Figure 1 a and b and save to PDF

p.suppfig1 <- p.fly.cumulative + labs(tag = "a") + p.encode.cumulative + labs(tag = "b") + plot_annotation(tag_levels = "a")
ggsave(p.suppfig1, file = "suppl_figure1.pdf", width = 11, height = 5)

# ------------------------------------------------------------------------
# Suppl Figure 2b: Abundance of BSJ k-mers across all samples

p.suppfig2b <- encode.dump %>%
  group_by(kmer) %>%
  summarise(totalcount = sum(count)) %>%
  group_by(totalcount) %>%
  summarise(kmers_per_totalcount = n_distinct(kmer)) %>%
  ggplot() +
  geom_point(aes(x = totalcount, y = kmers_per_totalcount)) +
  scale_y_log10() + scale_x_log10() +
  annotation_logticks() +
  theme_bw() +
  ylab("Number of BSJ k-mers") +
  xlab("Abundance")

# ------------------------------------------------------------------------
# combine Suppl Figure 2 a and b and save to PDF

p.suppfig2 <- p.suppfig2a + labs(tag = "a") + p.suppfig2b + labs(tag = "b") + plot_annotation(tag_levels = "a")
ggsave(p.suppfig2, file = "suppl_figure2.pdf", width = 11, height = 5)

# ------------------------------------------------------------------------
# Suppl Figure 3b: Number of samples in which BSJ k-mers are observed

p.suppfig3b <- encode.dump %>%
  group_by(kmer) %>%
  summarise(n_samples = n_distinct(ID)) %>%
  group_by(n_samples) %>%
  summarise(n_kmer = n_distinct(kmer)) %>%
  ggplot() +
  geom_point(aes(y=n_kmer,x=n_samples)) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks() +
  theme_bw() +
  ylab("Number of BSJ k-mers") +
  xlab("Number of samples")

# ------------------------------------------------------------------------
# combine Suppl Figure 3 a and b and save to PDF

p.suppfig3 <- p.suppfig3a + labs(tag = "a") + p.suppfig3b + labs(tag = "b") + plot_annotation(tag_levels = "a")
ggsave(p.suppfig3, file = "suppl_figure3.pdf", width = 11, height = 5)

# ------------------------------------------------------------------------
# Suppl Figure 5: UMAP dimension reduction

encode.dump.matrix <- acast(encode.dump, ID ~ kmer, value.var = "rpm", fill = 0.0)
set.seed(6423)
encode.umap <- umap(encode.dump.matrix)
df.encode.umap <- as.data.frame(encode.umap$layout)
colnames(df.encode.umap) <- c("Dim1", "Dim2")
df.encode.umap$ID <- rownames(df.encode.umap)

df.join <- df.encode.umap %>%
  inner_join(encode.exps, by = "ID") %>%
  inner_join(organs, by = c("Biosample.term.id" = "ID")) %>%
  inner_join(systems, by = c("Biosample.term.id" = "ID")) %>%
  mutate(
    system1 = case_when(
      str_detect(systems, ",") ~ "multiple systems",
      TRUE ~ systems
    ),
    organ1 = case_when(
      str_detect(organs, ",") ~ "multiple organs",
      TRUE ~ organs
    )
  )

names_systems <- df.join %>% filter(!system1 %in% as.list(c(NA, "multiple systems"))) %>% arrange(system1) %>% pull(system1) %>% unique()
names_systems <- c(unlist(names_systems), "multiple systems")
cols <- c(hue_pal()(length(names_systems) - 1), "grey50")
names(cols) <- names_systems

p.encode.umap <- df.join %>%
  ggplot() +
  geom_point(aes(x = Dim1, y = Dim2, shape = Biosample.type, color = system1)) +
  scale_color_manual(values = cols, na.value = "black", breaks = c(names_systems, NA), name = "System") +
  scale_shape_discrete(name = "Type") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.direction = "vertical"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = rel(0.8)),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  guides(color = guide_legend(ncol = 3)) +
  theme(aspect.ratio = 1)
ggsave(p.encode.umap, file = "suppl_figure5.pdf", width = 7, height = 9)
