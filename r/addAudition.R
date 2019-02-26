library(dplyr)

addAudition <- function(seqs, threshold, include_ncbi=F){
        lista2 = list()
        pb <- txtProgressBar(min = 0, max = length(seqs), style = 3, char = "*")

        for(i in 1:length(seqs)){

                Sys.sleep(0.0001)

                tmp = ID_engine(seqs[i], db = "COX1_SPECIES")


                tmp = tmp[[1]] %>%
                        select(ID, taxonomicidentification, similarity ) %>%
                        filter(grepl(" ", taxonomicidentification)) %>%
                        mutate(similarity = as.numeric(as.character(similarity)))

                if(!grepl("GenBank", as.character(tmp[1,]$ID))){
                        tmp = tmp %>%
                                filter(!grepl(" sp.",taxonomicidentification))
                }

                if(tmp[1,]$similarity < threshold || grepl("GenBank", as.character(tmp[1,]$ID)) ){

                        lista2[[i]] = data.frame(Match = "Any",
                                                 Species = "",
                                                 Grades = "",
                                                 Observations =  if(grepl("RID", as.character(tmp[1,]$ID))){
                                                         paste("RID not available.")}
                                                 else{
                                                         paste(
                                                                 if(grepl("GenBank", as.character(tmp[1,]$ID))){
                                                                         "GenBank: "}
                                                                 else{"BOLD: "},
                                                                 tmp$taxonomicidentification[1]," (", tmp$similarity[1], "), ",
                                                                 tmp$taxonomicidentification[2]," (", tmp$similarity[2], "), ",
                                                                 tmp$taxonomicidentification[3]," (", tmp$similarity[3],").",
                                                                 sep = "")}
                        )
                }
                else{
                        tmp = tmp %>% filter(similarity > threshold)
                        barcodes = sort(
                                table(as.character(tmp$taxonomicidentification)),
                                decreasing = T)
                        if(length(barcodes) > 1){
                                vec <- vector('character')
                                for(k in 1:length(barcodes)){
                                        vec[k] = paste(names(barcodes[k])," (",barcodes[k],")", sep = "")
                                }
                                lista2[[i]] = data.frame(Match = "Ambiguous",
                                                         Species = paste(vec, collapse = ", "),
                                                         Grades = paste(paste(AuditionBarcodes(species = names(barcodes),
                                                                                               matches = sum(barcodes),
                                                                                               include_ncbi = include_ncbi)$Grades,
                                                                              collapse = ", ")," respectively.",sep = ""),
                                                         Observations = "")
                        }
                        else{
                                lista2[[i]] = data.frame(Match = "Unique",
                                                         Species = paste(names(barcodes)),
                                                         AuditionBarcodes(species = names(barcodes),
                                                                          matches = sum(barcodes),
                                                                          include_ncbi = include_ncbi))
                        }
                        setTxtProgressBar(pb, i)
                }
        }
        close(pb)
        return(data.frame(Samples = names(seqs), do.call('rbind', lista2)))
}
