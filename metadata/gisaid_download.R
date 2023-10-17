library('GISAIDR')
library('dplyr')
library('httr')
library('XML')
library("gsubfn")
library('tidyr')
library('seqinr')
library('stringr')
library('collections')

#login into GISAID
username <- Sys.getenv("GISAIDR_USERNAME")
password <- Sys.getenv("GISAIDR_PASSWORD")
credentials <- login(username = username, password = password)
df_ids <- GISAIDR::query(credentials = credentials, location = "North America / USA / Rhode Island", nrows = 50000, fast = TRUE)

#download the dataframe
gisaid_id_split <- split(df_ids$accession_id, ceiling(seq_along(df_ids$accession_id) / 3000))

for (i in seq_along(gisaid_id_split)) {
  current_iteration <- i
  current_ids <- gisaid_id_split[i]
  print(current_iteration)
  full_df <- download(credentials = credentials, list_of_accession_ids = current_ids, get_sequence=TRUE)
  export_fasta(full_df, out_file_name = paste("gisaid_", current_iteration, ".fasta", sep=""))
  full_df <- dplyr::select(full_df, -c(sequence))
  write.csv(full_df, paste("gisaid_", current_iteration, ".csv", sep=""))
}

#combine csv files
files <-  list.files(pattern = "gisaid_[^.]+.csv$")
df <-  read.csv(files[1])
for (f in files[-1]) df <- rbind(df, read.csv(f))
write.csv(df, "gisaid.csv", row.names=FALSE, quote=TRUE)

project_details <- GET(url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA744530[BioProject]&retmax=5000')
project_details <- xmlParse(content(project_details))

ids <- xmlToDataFrame(nodes=getNodeSet(project_details, "//Id"))
ids <- ids$text

# add additional PRJNA911596
project_details <- GET(url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA911596[BioProject]&retmax=5000')
project_details <- xmlParse(content(project_details))
ids2 <- xmlToDataFrame(nodes=getNodeSet(project_details, "//Id"))
ids2 <- ids2$text

combined_ids <- c(ids, ids2)

id_split <- split(combined_ids, ceiling(seq_along(combined_ids) / 200))

total_ri_ids <- list()
total_df <- NULL
for (i in id_split) {
  first_ids <- paste(i, collapse=',')
  url <- paste('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=', first_ids, sep='')
  print("Fetching url:")
  print(url)

  sra_content <- GET(url = url)
  sra_content <- xmlParse(content(sra_content))
  sra_ids <- xmlToDataFrame(nodes=getNodeSet(sra_content, "//DocSum/Item[@Name='ExpXml']"))
  run_ids <- xmlToDataFrame(node=getNodeSet(sra_content, "//DocSum/Item[@Name='Runs']"))

  ri_ids <- list()
  biosample_ids <- list()
  biosample_host_subject_ids <- list()
  ri_run_names <- list()

  for (raw_text in sra_ids$text) {
    library_name <- strapplyc(raw_text, "<LIBRARY_NAME>(.*?)</LIBRARY_NAME>", simplify = c)
    biosample_name <- strapplyc(raw_text, "<Biosample>(.*?)</Biosample>", simplify = c)
    ri_ids <- append(ri_ids, library_name)
    biosample_ids <- append(biosample_ids, biosample_name)
  }

  for (raw_text in run_ids$text) {
    run_name <- strapplyc(raw_text, "<Run acc=\"(.*?)\"", simplify = c)
    ri_run_names <- append(ri_run_names, run_name)
    total_ri_ids <- append(total_ri_ids, run_name)
  }

  stripped_biosample_ids <- paste(str_sub(biosample_ids, 5), collapse = ",")
  biosample_url <- paste('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=biosample&id=', stripped_biosample_ids, sep='')
  biosample_content <- GET(url = biosample_url)
  print("Fetching url:")
  print(biosample_url)
  biosample_content <- xmlParse(content(biosample_content))
  biosample_datas <- xmlToDataFrame(nodes=getNodeSet(biosample_content, "//DocumentSummary/SampleData"))

  for (raw_text in biosample_datas$text) {
    host_subject_id <- strapplyc(raw_text, '<Attribute attribute_name="host_subject_id" harmonized_name="host_subject_id" display_name="host subject id">(.*?)</Attribute>', simplify = c)
    #print(host_subject_id)
    if (is.null(host_subject_id)) {
        biosample_host_subject_ids <- append(biosample_host_subject_ids, "")
    } else {
        biosample_host_subject_ids <- append(biosample_host_subject_ids, host_subject_id)
    }
  }

  biosample_host_subject_ids <- biosample_host_subject_ids <- vector("character" , length(biosample_ids))
  df = data.frame(unlist(ri_ids),unlist(ri_run_names),unlist(biosample_ids), unlist(biosample_host_subject_ids))
  names(df) = c("ri_library","sra_run","sra_biosample", "sra_biosample_host_subject_id")

  if (is.null(total_df)) {
    total_df = df
  } else {
    total_df = rbind(df,total_df)
  }

  Sys.sleep(10) ## Skeep 2 seconds
}

df <-  read.csv("gisaid.csv")
df <- df %>% separate(strain, sep="/", into = c(NA, NA, "ri_library", NA), remove=FALSE)
merged_df = merge(x = df, y = total_df, by.x = "ri_library", by.y="ri_library", all.x = TRUE)
merged_df <- dplyr::select(merged_df, -c(ri_library))

write.csv(merged_df, "gisaid.csv", row.names=FALSE, quote=TRUE)
write.table(merged_df, "gisaid.tsv", row.names=FALSE, quote=TRUE, sep="\t")

#combine fasta file
system("cat *.fasta > gisaid.fasta")

#write run file
write(unlist(total_ri_ids),"sra_run.txt",sep="\n")

#swap out accession id for strains in fasta file
fasta_df <- read.fasta(file = "gisaid.fasta", as.string = TRUE, seqtype = "DNA", forceDNAtolower=FALSE)
metadata <- read.delim(file='gisaid.tsv', sep='\t', header=TRUE)
metadata_dict <- dict(items=metadata$strain, keys=metadata$accession_id)

fun <- function(name) {
  accession <- str_split(name, "@", simplify=TRUE)[3]

  #strain_name <- metadata_dict$get(accession)
  #strain_name

  if (accession %in% metadata$accession_id ) {
    strain_name <- metadata_dict$get(accession)
    strain_name
  } else {
    print(name)
    print(accession)
    accession
  }
}

new_names <- lapply(names(fasta_df), fun)

write.fasta(sequences = fasta_df, names = new_names, file.out = "gisaid.fasta")