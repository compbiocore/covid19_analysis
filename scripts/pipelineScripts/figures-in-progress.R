library(tidyverse)
library(scales)
library(lubridate)
library(readxl)
library(RColorBrewer)

##### Extract dates and enumerate variants of concern from pangolin results

### date format: x1 <- stamp("Mar 2, 2020")(ymd("2021-04-25")) # Apr 25, 2021

### pick up on this - first start with Figuring out which are the most abundant variants and use those for plotting
### make all other omicrons one color, all other VBMs another and then everything else is grey
### Or perhaps make it starts.with("BA.1") etc to cover all the major subsets of lineages - pull the cumulative from the top-lineages script


var_of_concern <- c("B.1.1.529", "BA.1", "BA.1.1", "BA.1.1.1", "BA.1.1.10", "BA.1.1.11", "BA.1.1.12", 
                    "BA.1.1.13", "BA.1.1.14", "BA.1.1.15", "BA.1.1.16", "BA.1.1.17", "BA.1.1.18", "BA.1.1.2", 
                    "BA.1.1.3", "BA.1.1.4", "BA.1.1.5", "BA.1.1.6", "BA.1.1.7", "BA.1.1.8", "BA.1.1.9", 
                    "BA.1.10", "BA.1.12", "BA.1.13", "BA.1.13.1", "BA.1.14", "BA.1.14.1", "BA.1.14.2", 
                    "BA.1.15", "BA.1.15.1", "BA.1.15.2", "BA.1.15.3", "BA.1.16", "BA.1.16.1", "BA.1.16.2", 
                    "BA.1.17", "BA.1.17.1", "BA.1.17.2", "BA.1.18", "BA.1.19", "BA.1.20", "BA.1.21.1", "BA.1.23", 
                    "BA.1.3", "BA.1.4", "BA.1.5", "BA.1.8", "BA.1.9", "BA.2", "BA.2.1", "BA.2.10", "BA.2.10.1", 
                    "BA.2.10.2", "BA.2.10.3", "BA.2.10.4", "BA.2.11", "BA.2.12", "BA.2.12.1", "BA.2.12.2", 
                    "BA.2.13", "BA.2.13.1", "BA.2.14", "BA.2.15", "BA.2.17", "BA.2.18", "BA.2.19", "BA.2.2", 
                    "BA.2.2.1", "BA.2.20", "BA.2.21", "BA.2.22", "BA.2.23", "BA.2.23.1", "BA.2.24", "BA.2.25", 
                    "BA.2.25.1", "BA.2.26", "BA.2.27", "BA.2.28", "BA.2.29", "BA.2.3", "BA.2.3.1", "BA.2.3.10", 
                    "BA.2.3.11", "BA.2.3.12", "BA.2.3.13", "BA.2.3.14", "BA.2.3.15", "BA.2.3.16", "BA.2.3.17", 
                    "BA.2.3.18", "BA.2.3.2", "BA.2.3.20", "BA.2.3.21", "BA.2.3.4", "BA.2.3.5", "BA.2.3.6", 
                    "BA.2.3.7", "BA.2.3.8", "BA.2.3.9", "BA.2.30", "BA.2.31", "BA.2.31.1", "BA.2.32", 
                    "BA.2.33", "BA.2.34", "BA.2.35", "BA.2.36", "BA.2.37", "BA.2.38", "BA.2.38.1", "BA.2.38.2", 
                    "BA.2.38.3", "BA.2.39", "BA.2.4", "BA.2.40", "BA.2.40.1", "BA.2.41", "BA.2.42", "BA.2.43", 
                    "BA.2.44", "BA.2.45", "BA.2.46", "BA.2.47", "BA.2.48", "BA.2.49", "BA.2.5", "BA.2.50", "BA.2.51", 
                    "BA.2.52", "BA.2.53", "BA.2.54", "BA.2.55", "BA.2.56", "BA.2.57", "BA.2.58", "BA.2.59", "BA.2.6", 
                    "BA.2.60", "BA.2.61", "BA.2.62", "BA.2.63", "BA.2.64", "BA.2.66", "BA.2.67", "BA.2.68", "BA.2.7", 
                    "BA.2.70", "BA.2.71", "BA.2.72", "BA.2.73", "BA.2.74", "BA.2.75", "BA.2.75.2", "BA.2.75.3", 
                    "BA.2.75.4", "BA.2.75.5", "BA.2.75.9", "BA.2.76", "BA.2.76.1", "BA.2.77", "BA.2.78", "BA.2.79", 
                    "BA.2.79.1", "BA.2.8", "BA.2.80", "BA.2.81", "BA.2.82", "BA.2.9", "BA.2.9.1", "BA.2.9.2", 
                    "BA.2.9.3", "BA.2.9.4", "BA.2.9.5", "BA.2.9.6", "BA.2.9.7", "BA.3", "BA.3.1", "BA.4", "BA.4.1", 
                    "BA.4.1.1", "BA.4.1.10", "BA.4.1.2", "BA.4.1.3", "BA.4.1.4", "BA.4.1.5", "BA.4.1.6", "BA.4.1.8", 
                    "BA.4.1.9", "BA.4.2", "BA.4.3", "BA.4.4", "BA.4.5", "BA.4.6", "BA.4.6.1", "BA.4.6.2", "BA.4.6.3", 
                    "BA.4.6.4", "BA.4.6.5", "BA.4.7", "BA.5", "BA.5.1", "BA.5.1.1", "BA.5.1.10", "BA.5.1.11", 
                    "BA.5.1.12", "BA.5.1.15", "BA.5.1.18", "BA.5.1.2", "BA.5.1.20", "BA.5.1.21", "BA.5.1.22", 
                    "BA.5.1.23", "BA.5.1.24", "BA.5.1.25", "BA.5.1.27", "BA.5.1.28", "BA.5.1.29", "BA.5.1.3", 
                    "BA.5.1.30", "BA.5.1.31", "BA.5.1.32", "BA.5.1.4", "BA.5.1.5", "BA.5.1.6", "BA.5.1.7", 
                    "BA.5.1.8", "BA.5.1.9", "BA.5.10", "BA.5.10.1", "BA.5.11", "BA.5.2", "BA.5.2.1", "BA.5.2.11", 
                    "BA.5.2.12", "BA.5.2.13", "BA.5.2.14", "BA.5.2.16", "BA.5.2.18", "BA.5.2.19", "BA.5.2.2", 
                    "BA.5.2.20", "BA.5.2.21", "BA.5.2.22", "BA.5.2.23", "BA.5.2.24", "BA.5.2.25", "BA.5.2.26", 
                    "BA.5.2.27", "BA.5.2.28", "BA.5.2.29", "BA.5.2.3", "BA.5.2.31", "BA.5.2.32", "BA.5.2.33", 
                    "BA.5.2.34", "BA.5.2.35", "BA.5.2.37", "BA.5.2.4", "BA.5.2.42", "BA.5.2.43", "BA.5.2.44", 
                    "BA.5.2.45", "BA.5.2.47", "BA.5.2.48", "BA.5.2.49", "BA.5.2.6", "BA.5.2.7", "BA.5.2.8", 
                    "BA.5.2.9", "BA.5.3", "BA.5.3.1", "BA.5.3.2", "BA.5.3.3", "BA.5.3.4", "BA.5.3.5", "BA.5.5", 
                    "BA.5.5.1", "BA.5.5.2", "BA.5.5.3", "BA.5.6", "BA.5.6.1", "BA.5.6.2", "BA.5.6.3", "BA.5.6.4", 
                    "BA.5.8", "BA.5.9", "BC.1", "BC.2", "BD.1", "BE.1", "BE.1.1", "BE.1.1.1", "BE.1.1.2", "BE.1.2", 
                    "BE.1.2.1", "BE.1.4", "BE.1.4.1", "BE.1.4.2", "BE.1.4.3", "BE.1.4.4", "BE.10", "BE.2", "BE.3", 
                    "BE.4", "BE.4.1", "BE.4.1.1", "BE.4.2", "BE.5", "BE.6", "BF.1", "BF.1.1", "BF.10", "BF.10.1", 
                    "BF.11", "BF.11.1", "BF.11.2", "BF.11.3", "BF.11.4", "BF.11.5", "BF.12", "BF.13", "BF.14", 
                    "BF.15", "BF.16", "BF.17", "BF.19", "BF.2", "BF.21", "BF.22", "BF.25", "BF.26", "BF.27", "BF.28", "BF.29", "BF.30", "BF.31", "BF.31.1", "BF.32", "BF.34", "BF.4", "BF.5", "BF.5.1", "BF.6", "BF.7", "BF.7.1", "BF.7.10", "BF.7.11", "BF.7.12", "BF.7.13", "BF.7.13.1", "BF.7.13.2", "BF.7.14", "BF.7.14.1", "BF.7.3", "BF.7.4", "BF.7.4.1", "BF.7.4.2", "BF.7.5", "BF.7.5.1", "BF.7.6", "BF.7.7", "BF.7.8", "BF.7.9", "BF.8", "BF.9", "BG.1", "BG.2", "BG.3", "BG.4", "BG.5", "BH.1", "BJ.1", "BK.1", "BL.1", "BL.1.1", "BL.1.3", "BL.1.4", "BL.3", "BL.4", "BL.5", "BL.6", "BM.3", "BM.4", "BM.4.1", "BM.4.1.1", "BM.5", "BN.1", "BN.1.1", "BN.1.1.1", "BN.1.10", "BN.1.11", "BN.1.2", "BN.1.2.1", "BN.1.3", "BN.1.3.1", "BN.1.4", "BN.1.4.1", "BN.1.5", "BN.1.6", "BN.1.7", "BN.1.8", "BN.1.9", "BN.2", "BN.2.1", "BN.3", "BN.3.1", "BN.5", "BN.6", "BP.1", "BQ.1", "BQ.1.1", "BQ.1.1.1", "BQ.1.1.10", "BQ.1.1.11", "BQ.1.1.13", "BQ.1.1.14", "BQ.1.1.15", "BQ.1.1.16", "BQ.1.1.17", "BQ.1.1.18", "BQ.1.1.19", "BQ.1.1.2", "BQ.1.1.20", "BQ.1.1.21", "BQ.1.1.22", "BQ.1.1.23", "BQ.1.1.24", "BQ.1.1.25", "BQ.1.1.26", "BQ.1.1.27", "BQ.1.1.28", "BQ.1.1.3", "BQ.1.1.30", "BQ.1.1.31", "BQ.1.1.32", "BQ.1.1.35", "BQ.1.1.37", "BQ.1.1.38", "BQ.1.1.39", "BQ.1.1.4", "BQ.1.1.40", "BQ.1.1.41", "BQ.1.1.42", "BQ.1.1.43", "BQ.1.1.45", "BQ.1.1.47", "BQ.1.1.48", "BQ.1.1.49", "BQ.1.1.5", "BQ.1.1.6", "BQ.1.1.7", "BQ.1.1.8", "BQ.1.1.9", "BQ.1.10", "BQ.1.10.1", "BQ.1.11", "BQ.1.11.1", "BQ.1.12", "BQ.1.13", "BQ.1.13.1", "BQ.1.14", "BQ.1.15", "BQ.1.16", "BQ.1.17", "BQ.1.18", "BQ.1.2", "BQ.1.2.1", "BQ.1.20", "BQ.1.21", "BQ.1.23", "BQ.1.24", "BQ.1.25", "BQ.1.25.1", "BQ.1.26", "BQ.1.26.1", "BQ.1.3", "BQ.1.4", "BQ.1.5", "BQ.1.6", "BQ.1.7", "BQ.1.8", "BQ.1.8.2", "BQ.1.9", "BQ.2", "BR.1", "BR.1.2", "BR.2", "BR.2.1", "BR.3", "BR.4", "BR.5", "BT.1", "BT.2", "BU.1", "BU.2", "BU.3", "BV.1", "BW.1", "BW.1.1", "BZ.1", "CA.1", "CA.2", "CA.3", "CA.3.1", "CA.5", "CA.6", "CA.7", "CC.1", "CD.1", "CD.2", "CF.1", "CH.1", "CH.1.1", "CH.1.1.1", "CH.1.1.2", "CH.1.1.3", "CH.1.1.5", "CH.1.1.6", "CH.1.1.7", "CH.1.1.8", "CH.1.1.9", "CH.2", "CK.1", "CK.1.2", "CK.2", "CK.2.1", "CK.2.1.1", "CK.3", "CM.1", "CM.10", "CM.11", "CM.12", "CM.2", "CM.3", "CM.4", "CM.4.1", "CM.5", "CM.5.2", "CM.6", "CM.6.1", "CM.8", "CM.8.1", "CM.9", "CN.1", "CN.2", "CP.1", "CP.1.1", "CP.2", "CP.3", "CP.4", "CP.6", "CQ.1", "CQ.1.1", "CQ.2", "CR.1", "CR.1.1", "CR.1.3", "CR.2", "CS.1", "CY.1", "DB.2", "DC.1", "DE.1", "DE.2", "DF.1", "DF.1.1", "DH.1", "DJ.1", "DJ.1.1", "DJ.1.1.1", "DJ.1.2", "DJ.1.3", "DK.1", "DL.1", "DN.1", "DN.1.1", "DN.1.1.2", "DQ.1", "DR.1", "DS.1", "DS.2", "DT.1", "DU.1", "DV.1", "DV.2", "DV.3", "DW.1", "DY.1", "XAA", "XAB", "XAE", "XAG", "XAJ", "XAK", "XAL", "XAM", "XAN", "XAQ", "XAS", "XAU", "XAV", "XAY", "XAY.1", "XAZ", "XBB", "XBB.1", "XBB.1.1", "XBB.1.10", "XBB.1.11", "XBB.1.3", "XBB.1.4", "XBB.1.4.1", "XBB.1.5", "XBB.1.5.2", "XBB.1.5.3", "XBB.1.6", "XBB.1.7", "XBB.1.9", "XBB.2", "XBB.2.1", "XBB.2.2", "XBB.2.3", "XBB.3", "XBB.3.1", "XBB.3.2", "XBB.3.3", "XBB.4", "XBB.6", "XBB.6.1", "XBC", "XBC.1", "XBC.1.1", "XBD", "XBE", "XBF", "XBJ", "XBM", "XBP", "XD", "XE", "XF", "XG", "XH", "XL", "XN", "XQ", "XS", "XV", "XW", "XZ")
# manual colors:
set1mod14 <- c("#FF3300","#6600CC","#CC0066","#CC9900","#FF9900","#0000FF","#3399FF","#00CCFF",
               "#FF99FF","#006633","#009900","#00CC33","#33FF99","#989898")

var_regions <- c("UK","South Africa","Brazil","California, USA","California, USA","New York, USA",
                 "New York, USA","New York, USA","Brazil","India","India","India","India")

# Pangolin QC only (without nextstrain QC)
# pangolin results from Oscar
pangolin_res <- read_csv("results/qc-passed.csv") %>%
  rename(seqName=strain, lineage=pangolin.lineage, coll_date=date) %>%
  select(seqName, lineage, coll_date)

# combine B.1.526 with B.1.526.2
sum_df <- pangolin_res %>%
  mutate(lineage = if_else((lineage == "B.1.526.2" | lineage == "B.1.526.3"),"B.1.526",lineage))

roundUp <- function(x) 10^ceiling(log10(x))

dat1 <- sum_df

var_concern_df <- as.data.frame(matrix(NA,0,6))

for(z in 1:length(var_of_concern)){
  loc_variant <- dat1 %>%
    filter(lineage == var_of_concern[z]) %>%
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

# number total with date range
report_total <- var_concern_df %>%
  filter(!is.na(total)) %>%
  select(variant,region,total,date_range)
colnames(report_total) <- c("Variant of concern/interest","Region Variant was Originally Identified","Identified total cases, n","Range of sampling dates")

f_var_total_out <- paste0("results/variants_of_concern_summary_oscar_",format(Sys.Date(),"%Y%b%d"),".xlsx")
openxlsx::write.xlsx(report_total,f_var_total_out)

####
# numbers per variant per week
rep_by_week <- var_concern_df %>%
  mutate(week = week(coll_date)) %>%
  group_by(variant,week) %>%
  summarise(n=sum(n_per_date)) %>%
  filter(!is.na(n))

rep_by_week$n <- as.numeric(rep_by_week$n)


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

# Figure 1
ggplot(dat6,aes(x=coll_date,y=N,group=Source)) +
  geom_line(aes(col=Source),size=1.2) +
  scale_color_manual(values=c("#0095d6","#19b24b","#FF9900","#05618a","#d6bc00","#d60f63")) + # AI panel High Contrast 2
  scale_x_date(breaks=date_breaks("1 months"),labels=date_format("%b %y")) +
  scale_y_continuous(name="Cumulative Sequences, n", breaks=seq(0,total_num,minor_breaks)) +
  xlab("Date of Sample") +
  theme_classic() +
  theme(axis.text= element_text(size=18),
        axis.text.x = element_text(size = 14, angle = 90,hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 18),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 15))

f_out <- paste("results/Fig-1B_n",nrow(dat1),"_cumulative_seqs_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 10, height = 10, dpi = 300)


#### graphs with variants of concern and variants of interest

f2_in <- read_csv("src/weeks_2021.csv")
incl_weeks <- f2_in$Date

total_per_week_tmp <- dat1 %>%
  filter(year(coll_date) >= 2021) %>%
  mutate(week = week(coll_date)) %>%
  group_by(week) %>%
  summarise(n=n())

last_week <- total_per_week_tmp$week[nrow(total_per_week_tmp)]
week <- c(1:last_week)
weeks_df <- data.frame(week)
total_per_week <- left_join(weeks_df,total_per_week_tmp,by="week")
total_per_week$n[is.na(total_per_week$n)] <- 0

total_per_week$variant <- "Total"
total_per_week$date_range <- NA
for(k in 1:nrow(total_per_week)){
  total_df <- f2_in %>%
    filter(MMWR_week == total_per_week$week[k])
  total_per_week$date_range[k] <- total_df$Date[1]
}
total_per_week$date_range <- factor(total_per_week$date_range, levels = f2_in$Date)


# create a data frame with existing variants using total_per_week as a base
var_per_week <- as.data.frame(matrix(NA,0,4))
for(i in 1:length(var_of_concern)){
  local_df <- rep_by_week %>%
    filter(variant == var_of_concern[i])
  loc_df_per_var <- as.data.frame(matrix(NA,nrow(total_per_week),4))
  loc_df_per_var$V1 <- var_of_concern[i]
  loc_df_per_var$V2 <- total_per_week$week
  loc_df_per_var$V4 <- total_per_week$date_range
  if(nrow(local_df) > 0){
    for(m in 1:nrow(total_per_week)){
      for(l in 1:nrow(local_df)){
        if(local_df$week[l] == m){
          loc_df_per_var$V3[m] <- local_df$n[l]
          break
        } else {
          loc_df_per_var$V3[m] <- 0
        }
      }
    }
  } else {
    loc_df_per_var$V3 <- 0
  }
  var_per_week <- rbind(var_per_week,loc_df_per_var)
}
colnames(var_per_week) <- c("variant","week","n","date_range")


# Figure 2: stacked bars with % of VOC, VOI and non-VOC/non-VOI 
var_per_week$perc <- NA
percent_var_per_week <- as.data.frame(matrix(NA,0,5))
for(i in 1:nrow(total_per_week)){
  loc_df_perc <- var_per_week %>%
    filter(week == i)
  variant <- "non-VOC/non-VOI"
  week <- i
  n <- total_per_week$n[i] - sum(loc_df_perc$n)
  date_range <- as.character(total_per_week$date_range[i])
  perc <- round(n/total_per_week$n[i]*100,1)
  non_voc_voi <- cbind(variant,week,n,date_range,perc)
  for(k in 1:nrow(loc_df_perc)){
    loc_df_perc$perc[k] <- round(loc_df_perc$n[k]/total_per_week$n[i]*100,1)
  }
  loc_df_perc <- rbind(loc_df_perc,non_voc_voi)
  percent_var_per_week <- rbind(percent_var_per_week,loc_df_perc)
}

percent_var_per_week$perc <- as.numeric(percent_var_per_week$perc)
percent_var_per_week$n <- as.numeric(percent_var_per_week$n)
percent_var_per_week$variant <- factor(percent_var_per_week$variant, level=unique(percent_var_per_week$variant))

# stacked bars with %
ggplot(percent_var_per_week) + 
  geom_bar(aes(x = date_range, y = perc, fill = variant),color="white",stat = "identity") +
  scale_fill_manual(values = set1mod14) +
  scale_x_discrete(labels = incl_weeks) +
  geom_text(data=total_per_week,aes(x = date_range, y = 100, label=n), size = 8, vjust = -0.5) +
  xlab("Collection date, weeks during 2021") +
  ylab("Percent of identified cases, %") +
  theme_classic() +
  theme(axis.text= element_text(size=18),
        axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=16),
        axis.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position="top")

f_out <- paste("results/Fig-2_set1mod14_Lineges_of_concern_in_RI_by_week_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 12, height = 10, dpi = 300)

# Figure 3: stacked bars with actual numbers
ggplot(percent_var_per_week) + 
  geom_bar(aes(x = date_range, y = n, fill = variant),color="white",stat = "identity") +
  scale_fill_manual(values = set1mod14) +
  scale_x_discrete(labels = incl_weeks) +
  xlab("Collection date, weeks during 2021") +
  ylab("Identified cases, n") +
  theme_classic() +
  theme(axis.text= element_text(size=18),
        axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=16),
        axis.title = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position="top")

f_out <- paste("results/Fig-3_set1mo14_Lineges_of_concern_in_RI_by_week_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 12, height = 10, dpi = 300)

# Figure 4: stacked bars with actual number but without non-VOC/non-VOI
df_var_per_week <- percent_var_per_week %>%
  filter(variant != "non-VOC/non-VOI") %>%
  select(-perc)

tmp_df <- df_var_per_week %>%
  group_by(week) %>%
  summarise(m_n=sum(n))

max_y <- round(max(tmp_df$m_n)+10,-1)
minor_br_y <- roundUp(max_y/10)/5

ggplot(df_var_per_week) + 
  geom_bar(aes(x = date_range, y = n, fill = variant),color="white",stat = "identity") +
  scale_fill_manual(values = set1mod14) +
  scale_x_discrete(labels = incl_weeks) +
  scale_y_continuous(name="Identified cases, n", breaks=seq(0,max_y,minor_br_y)) +
  xlab("Collection date, weeks during 2021") +
  theme_classic() +
  theme(axis.text= element_text(size=18),
        axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=16),
        axis.title = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position="top")

f_out <- paste("results/Fig-4_set1mod14_VOC-VOI_in_RI_by_week_",format(Sys.Date(),"%Y%b%d"),".pdf",sep = "")
ggsave(f_out, device = "pdf",width = 12, height = 10, dpi = 300)