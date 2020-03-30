library(httr)
library(jsonlite)
library(xml2)


geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genesymbol, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

server <- "https://rest.ensembl.org"


ext <- "/sequence/id/ENSG00000157764?"

list_fasta <- list()

for (i in geneIDs2$GENEID) {
  ext <- paste("/sequence/id/", i, "?", sep = "")
  r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  fasta_seq <- content(r)
  list_fasta[[i]] <- fasta_seq
}

r <- GET(paste(server, ext, sep = ""), content_type("text/x-fasta"))

stop_for_status(r)


print(content(r))