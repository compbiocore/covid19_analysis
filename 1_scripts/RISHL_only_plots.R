library(tidyverse)
library(scales)
library(lubridate)
library(readxl)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
day = args[1]
pth = args[2]

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

pangolin_res <- read_csv(paste0(pth, "/3_results/", day, "/qc-passed.csv")) %>%
  rename(seqName=strain, lineage=pangolin.lineage, coll_date=date) %>%
  select(seqName, lineage, coll_date) %>%
  filter(grepl("hCoV-19/USA/RI-RISHL", seqName))

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

f_var_total_out <- paste0(pth, "/3_results/", day, "/RISHL_only_variants_of_concern_summary_",format(Sys.Date(),"%Y%b%d"),".xlsx")
openxlsx::write.xlsx(report_total,f_var_total_out)

####
# numbers per variant per week
rep_by_month <- var_concern_df %>%
  mutate(month = month(coll_date)) %>%
  group_by(variant,m) %>%
  summarise(n=sum(n_per_date)) %>%
  filter(!is.na(n))

rep_by_month$n <- as.numeric(rep_by_month$n)



### sequences by date

dat2 <- dat1 %>%
  mutate(Source = if_else(str_detect(seqName,"_RKL_"),"Kantor Lab",
                          if_else(str_detect(seqName,"RISHL"),"RISHL",
                                  if_else((str_detect(seqName,"Broad") | str_detect(seqName,"CDCBI")),"Broad",
                                          if_else(str_detect(seqName,"CDC-"),"CDC","Others")))),
         value = sample(1:nrow(dat1)))

total_num <- round(nrow(dat2)+50,-2)
minor_breaks <- roundUp(total_num/10)/5

### Total
dat3 <- dat2 %>%  # total cumulative
  group_by(coll_date) %>%
  summarize(n=n()) %>%
  mutate(N=cumsum(n),Source="Total") %>%
  select(Source,coll_date,N)

dat4 <- dat2 %>%   # by source cumulative
  group_by(Source,coll_date) %>%
  summarize(n=n()) %>%
  mutate(N=cumsum(n)) %>%
  select(-n)

dat5 <- full_join(dat3,dat4)

dat5$Source <- factor(dat5$Source,levels = c("Broad","CDC","Kantor Lab","RISHL","Others","Total"))

#
starting_date <- "2020-12-01"

for(i in 1:length(unique(dat5$Source))){
  loc_df <- dat5 %>%
    filter(Source == unique(dat5$Source)[i])
  source_with_starting_date <- loc_df %>%
    filter(coll_date == starting_date)
  if(nrow(source_with_starting_date) >= 1){
    next
  } else {
    loc_df_before <- loc_df %>%
      filter(coll_date < starting_date)
    if(nrow(loc_df_before) > 0){
      starting_number <- loc_df_before$N[nrow(loc_df_before)]
      new_line <- cbind(as.character(unique(dat5$Source)[i]),starting_date,starting_number)
      colnames(new_line) <- c("Source","coll_date","N")
      dat5 <- rbind(dat5,new_line)
      
    }
  }
}


dat6 <- dat5 %>%
  arrange(Source,coll_date) %>%
  filter(coll_date >= starting_date)
dat6$N <- as.numeric(dat6$N)

#### graphs with variants of concern and variants of interest

total_per_month_tmp <- dat1 %>%
  filter(year(coll_date) >= 2021) %>%
  mutate(year = year(coll_date)) %>%
  mutate(month = month(coll_date)) %>%
  mutate(ym = as.Date(paste(year(coll_date), month(coll_date), "01", sep = "-"))) %>%
  group_by(ym,year, month) %>%
  summarise(nd=n())
total_per_month_tmp

total_per_month_tmp2 <- dat1 %>%
  filter(year(coll_date) >= 2021) %>%
  mutate(year = year(coll_date)) %>%
  mutate(month = month(coll_date)) %>%
  mutate(ym = as.Date(paste(year(coll_date), month(coll_date), "01", sep = "-"))) %>%
  group_by(ym) %>%
  summarise(nm=n())
total_per_month_tmp2

total_per_month_tmp2 <- total_per_month_tmp %>% select(ym)
total_per_month_tmp2$m <- 1:nrow(total_per_month_tmp2)
total_per_month_tmp2

total_per_month <- left_join(total_per_month_tmp, total_per_month_tmp2, by = "ym")
total_per_month$variant <- "Total"
total_per_month$ym <- (format(as.Date(total_per_month$ym), "%Y-%m"))


# create a data frame with existing variants using total_per_week as a base
var_per_month <- as.data.frame(matrix(NA,0,4))
for(i in 1:length(var_of_concern)){
  local_df <- rep_by_month %>%
    filter(variant == var_of_concern[i])
  loc_df_per_var <- as.data.frame(matrix(NA,nrow(total_per_month),4))
  loc_df_per_var$V1 <- var_of_concern[i]
  loc_df_per_var$V2 <- total_per_month$m
  loc_df_per_var$V4 <- total_per_month$ym
  if(nrow(local_df) > 0){
    for(q in 1:nrow(total_per_month)){
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
  var_per_month <- rbind(var_per_month,loc_df_per_var)
}
colnames(var_per_month) <- c("variant","m","n","ym")
#var_per_month$ym <- (format(as.Date(var_per_month$ym), "%Y-%m"))

# Figure 2: stacked bars with % of VOC, VBM and non-VOC/non-VBM 


var_per_month$perc <- NA
percent_var_per_month <- as.data.frame(matrix(NA,0,5))
for(i in 1:nrow(total_per_month)){
  loc_df_perc <- var_per_month %>%
    filter(m == i)
  variant <- "non-VOC/non-VBM"
  m <- i
  n <- total_per_month$nd[i] - sum(loc_df_perc$n)
  dates <- as.character(total_per_month$ym[i])
  perc <- round(n/total_per_month$nd[i]*100,1)
  non_voc_voi <- cbind(variant,m,n,dates,perc)
  for(k in 1:nrow(loc_df_perc)){
    loc_df_perc$perc[k] <- round(loc_df_perc$n[k]/total_per_month$nd[i]*100,1)
  }
  colnames(non_voc_voi) <- c("variant","m","n","ym","perc")
  loc_df_perc <- rbind(loc_df_perc,non_voc_voi)
  percent_var_per_month <- rbind(percent_var_per_month,loc_df_perc)
}

percent_var_per_month$perc <- as.numeric(percent_var_per_month$perc)
percent_var_per_month$n <- as.numeric(percent_var_per_month$n)
percent_var_per_month$variant <- factor(percent_var_per_month$variant, level=unique(percent_var_per_month$variant))

# stacked bars with %
ggplot(percent_var_per_month) + 
  geom_bar(aes(x = ym, y = perc, fill = variant),color="white",stat = "identity") +
  scale_fill_manual(values = set1mod20) +
  scale_x_discrete(labels = var_per_month$ym) +
  geom_text(data=total_per_month,aes(x = ym, y = 100, label=nd ), size = 3, vjust = -0.5 ) +
  xlab("Collection date, months starting January 2021") +
  ylab("RIDOH Percent of identified cases, %") +
  theme_classic() +
  theme(axis.text= element_text(size=6),
        axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=8),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position="top")

f_out <- paste0(pth, "/3_results/", day, "/Fig_RISHL_Percent_variants_by_month_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 12, height = 10, dpi = 300)
f_out <- paste0(pth, "/3_results/", day, "/Fig_RISHL_Percent_variants_by_month_",format(Sys.Date(),"%Y%b%d"),".png",sep = "")
ggsave(f_out, device = "png",width = 12, height = 10, dpi = 300)


# Figure 3: stacked bars with actual numbers
ggplot(percent_var_per_month) + 
  geom_bar(aes(x = ym, y = n, fill = variant),color="white",stat = "identity") +
  scale_fill_manual(values = set1mod20) +
  scale_x_discrete(labels = var_per_month$ym) +
  xlab("Collection date, months starting February 2021") +
  ylab("RIDOH Number of identified cases") +
  theme_classic() +
  theme(#axis.text= element_text(size=8),
    axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size=12, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.position="top")

f_out <- paste0(pth, "/3_results/", day, "/Fig_RISHL_Total_variants_by_month_",format(Sys.Date(),"%Y%b%d"),".png",sep = "")
ggsave(f_out, device = "png",width = 18, height = 8, dpi = 300)
f_out <- paste0(pth, "/3_results/", day, "/Fig_RISHL_Total_variants_by_month_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 12, height = 10, dpi = 300)

# Figure 4: stacked bars with actual number but without non-VOC/non-VBM
df_var_per_month <- percent_var_per_month %>%
  filter(variant != "non-VOC/non-VBM") %>%
  select(-perc)

tmp_df <- df_var_per_month %>%
  group_by(m) %>%
  summarise(m_n=sum(n))

max_y <- round(max(tmp_df$m_n)+10,-1)
minor_br_y <- roundUp(max_y/10)/5

ggplot(df_var_per_month) + 
  geom_bar(aes(x = ym, y = n, fill = variant),color="white",stat = "identity") +
  scale_fill_manual(values = set1mod20) +
  scale_x_discrete(labels = df_var_per_month$ym) +
  scale_y_continuous(name="RIDOH Number of identified cases", breaks=seq(0,max_y,minor_br_y)) +
  xlab("Collection date, months starting February 2021") +
  theme_classic() +
  theme(axis.text= element_text(size=18),
        axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=8),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position="top")

f_out <- paste0(pth, "/3_results/", day, "/Fig_RISHL_VOC-VBM_by_month_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 12, height = 10, dpi = 300)
f_out <- paste0(pth, "/3_results/", day, "/Fig_RISHL_VOC-VBM_by_month_",format(Sys.Date(),"%Y%b%d"),".png",sep = "")
ggsave(f_out, device = "png",width = 12, height = 10, dpi = 300)
