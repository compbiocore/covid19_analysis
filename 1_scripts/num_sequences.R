library(tidyverse)

setwd("/gpfs/data/ris3/dev/20230312/3_results/20230313/")
ri <- read_csv("qc-passed.csv") %>%
  mutate(
    step=1,
    Source=as.factor(case_when(
      startsWith(strain, "hCoV-19/USA/RI-Broad") ~ "Broad",
      startsWith(strain, "hCoV-19/USA/RI-CDC")   ~ "CDC",
      startsWith(strain, "hCoV-19/USA/RI-RISHL") ~ "RISHL",
      startsWith(strain, "hCoV-19/USA/RI_RKL")   ~ "Kantor Lab",
      TRUE                                       ~ "Other"))
  ) %>%
  arrange(date)
print(select(filter(ri, Source=="Other"), strain), n=100)
nseq <- nrow(ri)
earliest <- min(ri$date)
latest <- max(ri$date)

ri <- group_by(ri, date, Source) %>%
  summarise(step=sum(step)) %>%
  ungroup() %>%
  group_by(Source) %>%
  mutate(Cumulative=cumsum(step)) %>%
  ungroup()

write_csv(ri, "num-sequences.csv")

# Summarize by week
ri <- mutate(ri, week=lubridate::floor_date(date, unit="week")) %>%
  group_by(week, Source) %>%
  summarise(Cumulative=max(Cumulative)) %>%
  ungroup() %>%
  pivot_wider(names_from=Source, values_from=Cumulative)

print(ri)

ri <- tibble(week=seq.Date(from=lubridate::floor_date(earliest, unit="week"), to=lubridate::floor_date(latest, unit="week"), by="week")) %>% 
  left_join(ri, on="week") %>%
  fill(everything()) %>%
  replace(is.na(.), 0)

print(ri)

ri <- ri %>%
  pivot_longer(!week, names_to="Source", values_to="Cumulative") %>%
  mutate(Source=factor(Source, levels=c("CDC", "Broad", "Kantor Lab", "RISHL", "Other")))

print(ri)

g <- ggplot(data=ri) +
  geom_area(aes(x=week, y=Cumulative, fill=Source)) +
  geom_hline(yintercept=nseq, linetype="dashed", color="black", size=0.25) +
  geom_text(
    data=data.frame(x=earliest, y=c(1.02*nseq), label=c(paste(scales::comma(nseq), "RI SARS-CoV-2 sequences as of", strftime(latest, format="%B %-m, %Y")))),
    aes(x=x, y=y, label=label, vjust=0, hjust=0),
    size=2.5
  ) +
  labs(
    x="Date of Sample",
    y="Cumulative Number of Sequences",
    fill="Source Lab"
  ) +
  scale_x_date(
    breaks=waiver(),
    date_breaks="month",
    labels=waiver(),
    date_labels="%b %Y"
  ) +
  scale_y_continuous(breaks=seq(0, length(nseq), 500), position="right") +
  scale_fill_manual(
    values=c(
      "CDC"="#e41a1c",
      "Broad"="#377eb8",
      "Kantor Lab"="#4daf4a",
      "RISHL"="#984ea3",
      "Other"="darkgray"
    )
  ) +
  theme_classic() +
  theme(
    title=element_text(size=9),
    legend.position=c(0, 0.5),
    legend.justification=c("left", "center"),
    legend.title=element_text(size=9),
    legend.text=element_text(size=8),
    legend.key.size=unit(0.1, "in"),
    axis.line=element_blank(),
    axis.ticks.x=element_line(size=0.25),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=8, color="black", angle=90, hjust=1, vjust=0.5),
    axis.text.y=element_text(size=8, color="black"),
    panel.grid.major.y=element_line(color="gray", size=0.1)
  )

f_out <- paste("Fig_Cumulative_num_seqs_by_source_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 6, height = 6, dpi = 300)

