ri <- read_csv("concern-long.csv")### This is the wrong input file I think - look at the initial mutation file

tot_muts <- unique(ri$mutation)


mut_concern_df <- as.data.frame(matrix(NA,0,6))

for(z in 1:length(tot_muts)){ # filter through list of mutations of concern
  mut <- ri %>%           # assign the file dat1 to mut 
    filter(mutation == tot_muts[z]) %>%  #filter lineage for particular mutation of concern
    group_by(date) %>%
    summarise(n=n()) %>%
    arrange(date) 
  
  if(nrow(mut) > 0){
    moc_df <- as.data.frame(matrix(NA,nrow(mut),5))
    moc_df$V1 <- tot_muts[z]

    moc_df$V2[1] <- sum(mut$n)
    moc_df$V3 <- mut$date
    if(nrow(mut) > 1){
      moc_df$V4 <- paste(stamp("Mar 2")(ymd(mut$date[1])),"to", stamp("Mar 2, 2020")(ymd(mut$date[nrow(mut)])))
    } else {
      moc_df$V4 <- paste0(" ", stamp("Mar 2, 2020")(ymd(mut$date[1])))
    }
    moc_df$V5 <- mut$n
  } else {
    moc_df <- as.data.frame(matrix(NA,1,5))
    moc_df$V1[1] <- tot_muts[z]
    moc_df$V2[1] <- 0
    moc_df$V3[1] <- NA
    moc_df$V4[1] <- "-"
    moc_df$V5[1] <- NA
  }
  mut_concern_df <- rbind(mut_concern_df,moc_df)
}

colnames(mut_concern_df) <- c("mutation","total","date","date_range","n_per_date")
mut_concern_df$mutation <- factor(mut_concern_df$mutation, level=unique(mut_concern_df$mutation))
mut_concern_df <- mut_concern_df[!is.na(mut_concern_df$date),]
mocdf <- mut_concern_df %>%  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  mutate(ym = as.Date(paste(year(date), month(date), "01", sep = "-"))) 

umd <- as.data.frame(unique(mocdf$ym))
colnames(umd) <- c("ym")
umd$ym <- umd[order(as.Date(umd$ym, format="%Y/%m/%d")),]
umd <- umd %>% mutate(m = c(1:length(unique(mocdf$ym))))

mut_concern_df <- left_join(mocdf,umd, by = "ym")
mut_concern_df$ym <- (format(as.Date(mut_concern_df$ym), "%Y-%m"))

# numbers per variant per week
mut_by_month <- mut_concern_df %>%
  mutate(month = month(date)) %>%
  group_by(mutation,m) %>%
  summarise(n=sum(n_per_date)) %>%
  filter(!is.na(n))

mut_by_month$n <- as.numeric(mut_by_month$n)

#### graphs with variants of concern and variants of interest

tpm <- ri %>%
  filter(year(date) >= 2021) %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  mutate(ym = as.Date(paste(year(date), month(date), "01", sep = "-"))) %>%
  group_by(ym,year, month) %>%
  summarise(nd=n())
tpm$ym <- (format(as.Date(tpm$ym), "%Y-%m"))
tpm$m <- 1:nrow(tpm)
tpm$seqs <- "Total"
tpm

# create a data frame with existing variants using total_per_week as a base
mut_per_month <- as.data.frame(matrix(NA,0,4))
for(i in 1:length(tot_muts)){
  local_df <- mut_by_month %>%
    filter(mutation == tot_muts[i])
  loc_df_per_var <- as.data.frame(matrix(NA,nrow(tpm),4))
  loc_df_per_var$V1 <- tot_muts[i]
  loc_df_per_var$V2 <- tpm$m
  loc_df_per_var$V4 <- tpm$ym
  if(nrow(local_df) > 0){
    for(q in 1:nrow(tpm)){
      for(l in 1:nrow(local_df)){
        if(local_df$m[l] == q){
          loc_df_per_var$V3[q] <- local_df$n[l]
          break
        } else {
          loc_df_per_var$V3[q] <- 0
        }
      }
    }
  } else {
    loc_df_per_var$V3 <- 0
  }
  mut_per_month <- rbind(mut_per_month,loc_df_per_var)
}
colnames(mut_per_month) <- c("mutation","m","n","ym")
#mut_per_month$ym <- (format(as.Date(mut_per_month$ym), "%Y-%m"))

# Figure 2: stacked bars with % of VOC, VBM and non-VOC/non-VBM 


mut_per_month$perc <- NA
percent_mut_per_month <- as.data.frame(matrix(NA,0,5))
for(i in 1:nrow(tpm)){
  moc_df_perc <- mut_per_month %>%
    filter(m == i)
  mutation <- "mutation"
  m <- i
  n <- tpm$nd[i] - sum(moc_df_perc$n)
  dates <- as.character(tpm$ym[i])
  perc <- round(n/tpm$nd[i]*100,1)
 # non_voc_voi <- cbind(mutation,m,n,dates,perc)
  #moc_df_perc$perc[26] <- round(moc_df_perc$n[26]/sum(moc_df_perc$n)*100,1)
  for(k in 1:nrow(moc_df_perc)){
    moc_df_perc$perc[k] <- round(moc_df_perc$n[k]/sum(moc_df_perc$n)*100,1)
 #   moc_df_perc$perc[k] <- round(moc_df_perc$n[k]/tpm$nd[13]*100,1)
  }
  #colnames(non_voc_voi) <- c("mutation","m","n","ym","perc")
  #moc_df_perc <- rbind(moc_df_perc,non_voc_voi)
  percent_mut_per_month <- rbind(percent_mut_per_month,moc_df_perc)
}

percent_mut_per_month$perc <- as.numeric(percent_mut_per_month$perc)
percent_mut_per_month$n <- as.numeric(percent_mut_per_month$n)
percent_mut_per_month$mutation <- factor(percent_mut_per_month$mutation, level=unique(percent_mut_per_month$mutation))

# stacked bars with %
col51 <- c("#899b34",
          "#7d71f0",
          "#b5bf21",
          "#a24ac5",
          "#57bd4a",
          "#d73bab",
          "#7eb72d",
          "#5453c0",
          "#ddb026",
          "#9778e2",
          "#b0b034",
          "#c370d9",
          "#4ab25d",
          "#b23f9a",
          "#54c38b",
          "#e13d7c",
          "#8ebc68",
          "#e173d0",
          "#597624",
          "#5484e3",
          "#eb9519",
          "#655aa9",
          "#cf9e32",
          "#5790d3",
          "#db7821",
          "#4cb6db",
          "#d14e2a",
          "#46ab94",
          "#d23f8c",
          "#43854b",
          "#d63f53",
          "#6e81b6",
          "#bb8b37",
          "#b297e1",
          "#935e24",
          "#da86ce",
          "#818b51",
          "#86579c",
          "#b59e57",
          "#a25190",
          "#e38853",
          "#bc85bb",
          "#c65543",
          "#b47699",
          "#d19f73",
          "#c16093",
          "#925937",
          "#cb4868",
          "#d07669",
          "#bd5c7d",
          "#bc6871")

ggplot(percent_mut_per_month) + 
  geom_bar(aes(x = ym, y = perc, fill = mutation),color="white",stat = "identity") +
  scale_fill_manual(values = col51) +
  scale_x_discrete(labels = mut_per_month$ym) +
  #geom_text(data=tpm,aes(x = ym, y = 100, label=nd ), angle=30, size = 3.5, hjust = -0.5 ) +
  xlab("Collection date, months starting January 2021") +
  ylab("Percent of spike protein mutations, %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))
      #  legend.position="top")


f_out <- paste("Fig_Percent_RI_Spike_protein_mutations_by_month_",format(Sys.Date(),"%Y%b%d"),".png",sep = "")
ggsave(f_out, device = "png",width = 18, height = 7, dpi = 300)
f_out <- paste("Fig_Percent_RI_Spike_protein_mutations_by_month_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 12, height = 10, dpi = 300)


