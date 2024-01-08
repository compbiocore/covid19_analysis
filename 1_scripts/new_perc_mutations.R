library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
day = args[1]
pth = args[2]

ridata <- read_csv(paste0(pth , "/3_results/", day , "/qc-passed.csv")) %>%
  rename(seqName=strain, lineage=pangolin.lineage, coll_date=date, mutations=nextclade.aaSubstitutions) %>%
  select(seqName, lineage, coll_date, mutations) %>% 
  separate_rows(mutations, sep = ",")


setwd(paste0(pth, "/2_metadata"))
concernMuts <- read.table("mutations_of_concern.txt", sep="\t", header = TRUE) %>% # table is all lineages copied from https://www.ecdc.europa.eu/en/covid-19/variants-concern with headers manually assigned for ease
  select(WHO_lineage, lineage, mutations) %>%
  separate_rows(mutations, sep = ", ") %>%
  group_by(mutations) %>%
  summarise(lineage=paste(lineage, collapse = ', '), .groups = 'drop')
concernMuts$mutations <- paste("S", concernMuts$mutations, sep = ":")
setwd(paste0(pth, "/3_results/", day))
write.table(concernMuts, file = paste0(pth , "/3_results/", day , "/concernMuts.txt"), sep = "\t",
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

# Figure 1: stacked bars with % of spike mutations of interest over time


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
col68 <- c("#8d476f", "#52c249", "#a053d2", "#8fbd39", "#5156cd", "#c4b83a", "#6e77f3", "#3e8f2b", "#dc63d3", "#61c36c", "#e0378f", "#4ed09b", "#a93b9c", "#8ea845", "#6d53b8", "#de8f2b", "#3c66c4", "#9e8823", "#6089ec", "#ea6c39", "#659ce8", "#cb3c28", "#3bbcc3", "#de3657", "#4aa364", "#d484e0", "#628124", "#9d82e1", "#606818", "#6d4599", "#96be77", "#a62c69", "#61c3a6", "#bb3860", "#3d7f43", "#df67af", "#306a3c", "#eb678d", "#379577", "#b14140", "#5bc0eb", "#ad5121", "#4797c8", "#dda955", "#505099", "#aa742c", "#5061a8", "#bbaf6c", "#9461b1", "#818b4b", "#784389", "#56642b", "#cd90cf", "#206e54", "#e46e67", "#4c78b7", "#db966a", "#68649c", "#866a32", "#a7a0e1", "#8a5229", "#9d5e92", "#ae6053", "#e292b6", "#9a4657", "#e18889", "#8b4a5f", "#c3618b")

percent_mut_per_month_start_2021 <- percent_mut_per_month %>% filter(ym >= "2021-01")
plen <- unique(percent_mut_per_month_start_2021$ym)
percent_mut_per_month_start_2021$m <- 1:length(plen)

ggplot(percent_mut_per_month_start_2021) + 
  geom_bar(aes(x = ym, y = perc, fill = mutation),color="white",stat = "identity") +
  scale_fill_manual(values = col68) +
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

state_name <- Sys.getenv("GISAIDR_STATE")
state_abbr <- state.abb[grep(state_name, state.name)]
f_name <- paste("/Fig_Percent_", state_abbr, "_Spike_protein_mutations_by_month_", sep="")

f_out <- paste0(pth , "/3_results/", day , f_name ,format(Sys.Date(),"%Y%b%d"),".png",sep = "")
ggsave(f_out, device = "png",width = 18, height = 7, dpi = 300)
f_out <- paste0(pth , "/3_results/", day , f_name,format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 18, height = 10, dpi = 300)
