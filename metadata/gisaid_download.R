library('GISAIDR')
library('dplyr')

#login into GISAID
username <- Sys.getenv("GISAIDR_USERNAME")
password <- Sys.getenv("GISAIDR_PASSWORD")
credentials <- login(username = username, password = password)
df_ids <- query(credentials = credentials, location = "North America / USA / Rhode Island", nrows = 50000, fast = TRUE)

#download the dataframe
gisaid_id_split <- split(df_ids$accession_id, ceiling(seq_along(df_ids$accession_id) / 3000))

for (i in seq_along(gisaid_id_split)) {
  current_iteration <- i
  current_ids <- gisaid_id_split[i]
  print(current_iteration)
  full_df <- download(credentials = credentials, list_of_accession_ids = current_ids, get_sequence=TRUE)
  export_fasta(full_df, out_file_name = paste("gisaid_", current_iteration, ".fasta", sep=""))
  full_df <- select(full_df, -c(sequence))
  write.csv(full_df, paste("gisaid_", current_iteration, ".csv", sep=""))
}

#combine csv files
files <- list.files(pattern = "\\.csv$")
df <-  read.csv(files[1])
for (f in files[-1]) df <- rbind(df, read.csv(f))
write.csv(df, "gisaid.csv", row.names=FALSE, quote=FALSE)

#combine fasta file
system("cat *.fasta > gisaid.fasta")
