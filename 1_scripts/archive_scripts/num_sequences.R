library(tidyverse)

day=$(date "+%Y%m%d")

setwd("../3_results/${day}/")
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






#########################################################################################################################


library(tidyverse)
library(scales)
library(lubridate)
library(readxl)
library(RColorBrewer)

##### Extract dates and enumerate variants of concern from pangolin results

### date format: x1 <- stamp("Mar 2, 2020")(ymd("2021-04-25")) # Apr 25, 2021

## Without BA.3 since it was never found in the data
var_of_concern <- c("B.1.1.529", "BA.1 Lineage", "BA.2 Lineage",  
                    "BA.4 Lineage", "BA.5 Lineage", "Other Omicron Lineage", 
                    "Recombinant Omicron Lineage", "Alpha Lineage", "Beta Lineage",
                    "Gamma Lineage", "Delta Lineage", "Eta Lineage", "Theta Lineage", "Kappa Lineage",
                    "Iota Lineage", "Zeta Lineage", "Mu Lineage")

set1mod20 <- c("#FF0000", "#FD4E4E", "#D55B5B", "#CD1F03", "#C70039",
               "#AA3927","#CC0066","#6600CC","#D8AAFB","#0000FF","#F88520","#00CCFF",
               "#F8F120","#93F820","#009900","#00CC33","#33FF99","#989898")

var_regions <- c("South Africa","South Africa","South Africa","South Africa",
                 "South Africa","South Africa","South Africa","UK","South Africa","Brazil",
                 "India", "Nigeria", "Philippines", "India", "USA", "Brazil", "Colombia")

# Pangolin QC only
# pangolin results from Oscar
pangolin_res <- read_csv("/gpfs/data/ris3/dev/20230312/3_results/20230313/qc-passed.csv") %>%
  rename(seqName=strain, lineage=pangolin.lineage, coll_date=date) %>%
  select(seqName, lineage, coll_date)

# combine Variants into larger categories
sum_df <- pangolin_res  %>%
  mutate(lineage = str_replace(string = lineage, pattern = "BA\\.1\\..+", replacement = "BA.1 Lineage"))  %>%
  mutate(lineage = str_replace(string = lineage, pattern = "BA\\.2\\..+", replacement = "BA.2 Lineage"))   %>%
  mutate(lineage = str_replace(string = lineage, pattern = "BA\\.3\\..+", replacement = "BA.3 Lineage"))  %>% 
  mutate(lineage = str_replace(string = lineage, pattern = "BA\\.4\\..+", replacement = "BA.4 Lineage"))  %>%
  mutate(lineage = str_replace(string = lineage, pattern = "BA\\.5\\..+", replacement = "BA.5 Lineage"))  %>% 
  mutate(lineage = str_replace(string = lineage, pattern = "AY\\..+", replacement = "Delta Lineage"))  %>%
  mutate(lineage = str_replace(string = lineage, pattern = "^X.+", replacement = "Recombinant Omicron Lineage"))  %>%
  mutate(lineage = str_replace(string = lineage, pattern = "^[B-D][C-Z]\\..+", replacement = "Other Omicron Lineage")) %>%
  mutate(lineage = str_replace(string = lineage, pattern = "^Q\\..+", replacement = "Alpha Lineage")) %>%
  mutate(lineage = if_else((lineage == "B.1.1.7"),"Alpha Lineage",lineage)) %>%
  mutate(lineage = str_replace(string = lineage, pattern = "^B.1.351.+", replacement = "Beta Lineage")) %>%
  mutate(lineage = str_replace(string = lineage, pattern = "^B.1.621.+", replacement = "Mu Lineage")) %>%
  mutate(lineage = if_else((lineage == "B.1.427" | lineage == "B.1.429"),"Epsilon Lineage",lineage)) %>%
  mutate(lineage = if_else((lineage == "B.1.525"),"Eta Lineage",lineage)) %>%  
  mutate(lineage = if_else((lineage == "B.1.526"),"Iota Lineage",lineage)) %>%  
  mutate(lineage = if_else((lineage == "B.1.617.1"),"Kappa Lineage",lineage)) %>% 
  mutate(lineage = if_else((lineage == "P.1"),"Gamma Lineage",lineage)) %>%
  mutate(lineage = if_else((lineage == "P.2"),"Zeta Lineage",lineage)) %>%
  mutate(lineage = if_else((lineage == "P.3"),"Theta Lineage",lineage))



roundUp <- function(x) 10^ceiling(log10(x))

dat1 <- sum_df

var_concern_df <- as.data.frame(matrix(NA,0,6))


for(z in 1:length(var_of_concern)){ # filter through list of variants of concern
  loc_variant <- dat1 %>%           # assign the file dat1 to loc_variant 
    filter(lineage == var_of_concern[z]) %>%  #filter lineage for particular variant of concern
    group_by(coll_date) %>%
    summarise(n=n()) %>%
    arrange(coll_date) 
  
  if(nrow(loc_variant) > 0){
    loc_df <- as.data.frame(matrix(NA,nrow(loc_variant),6))
    loc_df$V1 <- var_of_concern[z]
    loc_df$V2 <- var_regions[z]
    loc_df$V3[1] <- sum(loc_variant$n)
    loc_df$V4 <- loc_variant$coll_date
    if(nrow(loc_variant) > 1){
      loc_df$V5 <- paste(stamp("Mar 2")(ymd(loc_variant$coll_date[1])),"to", stamp("Mar 2, 2020")(ymd(loc_variant$coll_date[nrow(loc_variant)])))
    } else {
      loc_df$V5 <- paste0(" ", stamp("Mar 2, 2020")(ymd(loc_variant$coll_date[1])))
    }
    loc_df$V6 <- loc_variant$n
  } else {
    loc_df <- as.data.frame(matrix(NA,1,5))
    loc_df$V1[1] <- var_of_concern[z]
    loc_df$V2 <- var_regions[z]
    loc_df$V3[1] <- 0
    loc_df$V4[1] <- NA
    loc_df$V5[1] <- "-"
    loc_df$V6[1] <- NA
  }
  var_concern_df <- rbind(var_concern_df,loc_df)
}


colnames(var_concern_df) <- c("variant","region","total","coll_date","date_range","n_per_date")
var_concern_df$variant <- factor(var_concern_df$variant, level=unique(var_concern_df$variant))
var_concern_df <- var_concern_df[!is.na(var_concern_df$coll_date),]
vocdf <- var_concern_df %>%  mutate(year = year(coll_date)) %>%
  mutate(month = month(coll_date)) %>%
  mutate(ym = as.Date(paste(year(coll_date), month(coll_date), "01", sep = "-"))) 

umd <- as.data.frame(unique(vocdf$ym))
colnames(umd) <- c("ym")
umd$ym <- umd[order(as.Date(umd$ym, format="%Y/%m/%d")),]
umd <- umd %>% mutate(m = c(1:length(unique(vocdf$ym))))

var_concern_df <- left_join(vocdf,umd, by = "ym")
var_concern_df$ym <- (format(as.Date(var_concern_df$ym), "%Y-%m"))



# number total with date range
report_total <- var_concern_df %>%
  filter(!is.na(total)) %>%
  select(variant,region,total,date_range)
colnames(report_total) <- c("Variant of concern/being monitored","Region Variant was Originally Identified","Identified total cases, n","Range of sampling dates")

f_var_total_out <- paste0("variants_of_concern_summary_",format(Sys.Date(),"%Y%b%d"),".xlsx")
openxlsx::write.xlsx(report_total,f_var_total_out)

####
# numbers per variant per week
rep_by_month <- var_concern_df %>%
  mutate(month = month(coll_date)) %>%
  group_by(variant,m) %>%
  summarise(n=sum(n_per_date)) %>%
  filter(!is.na(n))

rep_by_month$n <- as.numeric(rep_by_month$n)



#### sequences by date
#
#dat2 <- dat1 %>%
#  mutate(Source = if_else(str_detect(seqName,"_RKL_"),"Kantor Lab",
#                      if_else(str_detect(seqName,"RISHL"),"RISHL",
#                              if_else(str_detect(seqName,"Yale"),"Yale",
#                              if_else((str_detect(seqName,"Broad") | str_detect(seqName,"CDCBI")),"Broad",
#                                              if_else(str_detect(seqName,"ASC-"),"Aegis",
#                                                      if_else(str_detect(seqName,"CDC-"),"CDC","Others")))))),
#         value = sample(1:nrow(dat1)))
#
#total_num <- round(nrow(dat2)+50,-2)
#minor_breaks <- roundUp(total_num/10)/5
#
#### Total
#dat3 <- dat2 %>%  # total cumulative
#  group_by(coll_date) %>%
#  summarize(n=n()) %>%
#  mutate(N=cumsum(n),Source="Total") %>%
#  select(Source,coll_date,N)
#
#dat4 <- dat2 %>%   # by source cumulative
#  group_by(Source,coll_date) %>%
#  summarize(n=n()) %>%
#  mutate(N=cumsum(n)) %>%
#  select(-n)
#
#dat5 <- full_join(dat3,dat4)
#
#dat5$Source <- factor(dat5$Source,levels = c("Broad","CDC","Kantor Lab","RISHL","Others","Total"))
#
##
#starting_date <- "2020-12-01"
#
#for(i in 1:length(unique(dat5$Source))){
#  loc_df <- dat5 %>%
#    filter(Source == unique(dat5$Source)[i])
#  source_with_starting_date <- loc_df %>%
#    filter(coll_date == starting_date)
#  if(nrow(source_with_starting_date) >= 1){
#    next
#  } else {
#    loc_df_before <- loc_df %>%
#      filter(coll_date < starting_date)
#    if(nrow(loc_df_before) > 0){
#      starting_number <- loc_df_before$N[nrow(loc_df_before)]
#      new_line <- cbind(as.character(unique(dat5$Source)[i]),starting_date,starting_number)
#      colnames(new_line) <- c("Source","coll_date","N")
#      dat5 <- rbind(dat5,new_line)
#      
#    }
#  }
#}
#
#
#dat6 <- dat5 %>%
#  arrange(Source,coll_date) %>%
#  filter(coll_date >= starting_date)
#dat6$N <- as.numeric(dat6$N)
#
## Figure 1
#ggplot(dat6,aes(x=coll_date,y=N,group=Source)) +
#  geom_line(aes(col=Source),size=1.2) +
#  #scale_color_manual(values=c("#0095d6","#19b24b","#FF9900","#05618a","#d6bc00","#d60f63")) + # AI panel High Contrast 2
#  scale_x_date(breaks=date_breaks("1 months"),labels=date_format("%b %y")) +
#  scale_y_continuous(name="Cumulative Sequences, n", breaks=seq(0,total_num,minor_breaks)) +
#  xlab("Date of Sample") +
#  theme_classic() +
#  theme(axis.text= element_text(size=18),
#        axis.text.x = element_text(size = 14, angle = 90,hjust = 1, vjust = 0.5),
#        axis.title = element_text(size = 18),
#        legend.title = element_text(size=16),
#        legend.text = element_text(size = 15))
#
#f_out <- paste("Fig_Cumulative_seqs_MoreLabs_(n=",nrow(dat1),")_by_source_",format(Sys.Date(),"%Y%b%d"),".png",sep = "")
#ggsave(f_out, device = "png",width = 10, height = 10, dpi = 300)
#f_out <- paste("Fig_Cumulative_seqs_MoreLabs_(n=",nrow(dat1),")_by_source_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
#ggsave(f_out, device = "pdf",width = 10, height = 10, dpi = 300)
#
#
