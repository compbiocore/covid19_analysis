library(tidyverse)

  ### still to do:
  ### Defining lineage mutations: https://github.com/cov-lineages/pango-designation/tree/master/lineage_constellations
  
  ### Define specific mutation groups for each lineage and use those to produce plots
  
  ### Plot of specific lineages - redo plot with mutations of concern including all omicron mutations as well.

day=$(date "+%Y%m%d")

ridata <- read_csv("../3_results/${day}/qc-passed.csv") %>%
  rename(seqName=strain, lineage=pangolin.lineage, coll_date=date, mutations=nextclade.aaSubstitutions) %>%
  select(seqName, lineage, coll_date, mutations) %>% 
  separate_rows(mutations, sep = ",")

concernMuts <- read.table("../2_metadata/mutations_of_concern.txt", sep="\t", header=TRUE) %>% # table is all lineages copied from https://www.ecdc.europa.eu/en/covid-19/variants-concern with headers manually assigned for ease
  select(WHO_lineage, lineage, mutations) %>%
  separate_rows(mutations, sep = ", ") %>%
  group_by(mutations) %>%
  summarise(lineage=paste0(lineage, collapse = ', '), .groups = 'drop')
concernMuts$mutations <- paste("S", concernMuts$mutations, sep = ":")
write.table(concernMuts, file = "../3_results/${day}/concernMuts.txt", sep = "\t",
            row.names = FALSE)

all_total_mutations <- unique(ridata$mutations)

keyMuts <- semi_join(ridata, concernMuts, by = "mutations")
total_mutations <- unique(keyMuts$mutations)

mut_concern_df <- as.data.frame(matrix(NA,0,6))

for(z in 1:length(total_mutations)){ 
  mut <- keyMuts %>%          
    filter(mutations == total_mutations[z]) %>% 
    group_by(coll_date) %>%
    summarise(n=n()) %>%
    arrange(coll_date) 
  
  if(nrow(mut) > 0){
    moc_df <- as.data.frame(matrix(NA,nrow(mut),5))
    moc_df$V1 <- total_mutations[z]
    
    moc_df$V2[1] <- sum(mut$n)
    moc_df$V3 <- mut$coll_date
    if(nrow(mut) > 1){
      moc_df$V4 <- paste(stamp("Mar 2")(ymd(mut$coll_date[1])),"to", stamp("Mar 2, 2020")(ymd(mut$coll_date[nrow(mut)])))
    } else {
      moc_df$V4 <- paste0(" ", stamp("Mar 2, 2020")(ymd(mut$coll_date[1])))
    }
    moc_df$V5 <- mut$n
  } else {
    moc_df <- as.data.frame(matrix(NA,1,5))
    moc_df$V1[1] <- total_mutations[z]
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

tpm <- keyMuts %>%
  mutate(year = year(coll_date)) %>%
  mutate(month = month(coll_date)) %>%
  mutate(ym = as.Date(paste(year(coll_date), month(coll_date), "01", sep = "-"))) %>%
  group_by(ym,year, month) %>%
  summarise(nd=n())
tpm$ym <- (format(as.Date(tpm$ym), "%Y-%m"))
tpm$m <- 1:nrow(tpm)
tpm$seqs <- "Total"
tpm

# create a data frame with existing mutations using total_per_week as a base
mut_per_month <- as.data.frame(matrix(NA,0,4))
for(i in 1:length(total_mutations)){
  local_df <- mut_by_month %>%
    filter(mutation == total_mutations[i])
  loc_df_per_var <- as.data.frame(matrix(NA,nrow(tpm),4))
  loc_df_per_var$V1 <- total_mutations[i]
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
  for(k in 1:nrow(moc_df_perc)){
    moc_df_perc$perc[k] <- moc_df_perc$n[k]/sum(moc_df_perc$n)*100
  }
  percent_mut_per_month <- rbind(percent_mut_per_month,moc_df_perc)
}

percent_mut_per_month$perc <- as.numeric(percent_mut_per_month$perc)
percent_mut_per_month$n <- as.numeric(percent_mut_per_month$n)
percent_mut_per_month$mutation <- factor(percent_mut_per_month$mutation, level=unique(percent_mut_per_month$mutation))

# stacked bars with %
col67 <- c("#3f912e",
           "#8740bd",
           "#6eb729",
           "#b465e5",
           "#5bc453",
           "#516fef",
           "#a4be2b",
           "#3e56c6",
           "#bdad30",
           "#717aee",
           "#d9992a",
           "#6c4fb5",
           "#4bbb6d",
           "#d94db1",
           "#4ec88c",
           "#a43aa0",
           "#8f9c36",
           "#c973dc",
           "#45813b",
           "#9480e8",
           "#e2751e",
           "#4691eb",
           "#d5482e",
           "#34b9e1",
           "#d44248",
           "#45c4a8",
           "#dd3e77",
           "#87b25f",
           "#7b51a7",
           "#bd903b",
           "#476ac2",
           "#c16726",
           "#66a1e5",
           "#d26e45",
           "#3bbcc3",
           "#ce3b54",
           "#61a874",
           "#d34d92",
           "#307646",
           "#d680d4",
           "#6c671e",
           "#aa8edd",
           "#bba25b",
           "#505da4",
           "#ad6d37",
           "#4177b9",
           "#ce7058",
           "#5a9ac9",
           "#905734",
           "#a8a1de",
           "#808d52",
           "#b3579b",
           "#459576",
           "#88539a",
           "#dea47a",
           "#6a679f",
           "#a08253",
           "#d28ec9",
           "#d0706b",
           "#935f93",
           "#dc7c88",
           "#944471",
           "#eb85a6",
           "#8b4a5f",
           "#ce83ad",
           "#a24654",
           "#ad486c")

percent_mut_per_month_start_2021 <- percent_mut_per_month %>% filter(ym >= "2021-01")
plen <- unique(percent_mut_per_month_start_2021$ym)
percent_mut_per_month_start_2021$m <- 1:length(plen)

ggplot(percent_mut_per_month_start_2021) + 
  geom_bar(aes(x = ym, y = perc, fill = mutation),color="white",stat = "identity") +
  scale_fill_manual(values = col67) +
  scale_x_discrete(labels = plen) +
  #geom_text(data=tpm,aes(x = ym, y = 100, label=nd ), angle=30, size = 3.5, hjust = -0.5 ) + this would put total numbers at the top of each month. Dropped because #'s are <50,000 some months
  xlab("Collection date, months starting January 2021") +
  ylab("Percent of spike protein mutations, %") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size=12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))
#  legend.position="top")


f_out <- paste("../3_results/${day}/Fig_Percent_RI_Spike_protein_mutations_by_month_",format(Sys.Date(),"%Y%b%d"),".png",sep = "")
ggsave(f_out, device = "png",width = 18, height = 7, dpi = 300)
f_out <- paste("../3_results/${day}/Fig_Percent_RI_Spike_protein_mutations_by_month_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 12, height = 10, dpi = 300)
