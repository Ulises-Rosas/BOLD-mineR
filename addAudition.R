library(dplyr)

AuditionBarcodes<- function(species, matches){ ##function for only using with public data
        frames = lapply(species, function(x){

                meta.by.barcodes1 = SpecimenData(taxon = x) %>%
                        dplyr::select(processid, bin_uri, species_name, institution_storing) %>%
                        dplyr::mutate_if(is.factor, as.character) %>%
                        dplyr::filter(grepl("BOLD", bin_uri),
                                      !grepl("Mined from GenBank, NCBI", institution_storing),
                                      !grepl("*unvouchered", institution_storing))


                js0 = getURL(
                        paste("http://www.boldsystems.org/index.php/API_Tax/TaxonSearch?taxName=",
                              gsub(" ","%20", x), sep = "")
                        ) %>%
                        gsub('.*\"taxid\":', "", x = .) %>%
                        gsub(',\"taxon\".*', "", x = .) %>%
                        paste("http://www.boldsystems.org/index.php/API_Tax/TaxonData?taxId=", . ,
                              "&dataTypes=all", sep = "") %>% getURL(url = .) %>%
                        gsub('.*\"depositry\":\\{', "", x = .) %>%
                        gsub('\\}.*', "", x = .) %>% gsub('\"', "", x = .) %>%
                        strsplit(x = ., split = ",") %>% .[[1]] %>%
                        strsplit(x = ., split = "\\:") %>%
                        lapply(., function(x){
                                tmp = x[!grepl("Mined from GenBank", x[1]) &
                                                !grepl(" NCBI", x[1]) &
                                                !grepl("*unvouchered", x[1])]
                                data.frame(institutions = tmp[1], records = as.numeric(tmp[2]))
                                }) %>%
                        do.call("rbind", .) %>%
                        .[!is.na(.$records),]

                if(nrow(meta.by.barcodes1) == 0 && sum(js0$records, na.rm = T) == 0){
                        data.frame(Grades = "F",
                                   Observations = paste("There were ",
                                                        matches ,
                                                        " matches. Barcodes mined from GenBank, NCBI.",
                                                        sep = ""))
                        }
                else if(nrow(meta.by.barcodes1) <= 3 && sum(js0$records, na.rm = T) != 0){
                        data.frame(Grades = "D",
                                   Observations = paste("There were ",
                                                        matches ,
                                                        " matches. Insufficient data. Institution storing: ",
                                                        length(js0$institutions),
                                                        ". Specimen records: ",
                                                        sum(js0$records, na.rm = T),
                                                        sep = ""))
                        }
                else{
                        if(length(unique(meta.by.barcodes1$bin_uri)) > 1){

                                bin = lapply(unique(meta.by.barcodes1$bin_uri), function(x){
                                        SpecimenData(bin = x) %>%
                                                dplyr::select(species_name, institution_storing)  %>%
                                                dplyr::filter(
                                                        grepl(" ",species_name),
                                                        !grepl("Mined from GenBank, NCBI", institution_storing),
                                                        !grepl("*unvouchered", institution_storing),
                                                        !grepl(" sp.",species_name)
                                                ) %>%
                                                dplyr::group_by(species_name) %>% dplyr::summarise()
                                })

                                if(length(unique(do.call('rbind', bin)$species_name)) > 1){


                                        data.frame(Grades = "E**",
                                                   Observations = paste("There were ",
                                                                        matches ,
                                                                        " matches. Mixtured BIN and it's composed by species such as: ",
                                                                        paste(unique(do.call('rbind', bin)$species_name), collapse = ", "),
                                                                        sep = ""))
                                }
                                else{
                                        data.frame(Grades = "C*",
                                                   Observations = paste("There were ", matches,
                                                                        " matches. Assessment of intraspecific divergences is still needed.",
                                                                        sep = ""))
                                }

                        }
                        else{
                                unique.bin = SpecimenData(bin = unique(meta.by.barcodes1$bin_uri)) %>%
                                        dplyr::select(species_name, institution_storing) %>%
                                        dplyr::filter(grepl(" ",species_name),
                                                      !grepl("Mined from GenBank, NCBI", institution_storing),
                                                      !grepl("*unvouchered", institution_storing),
                                                      !grepl(" sp.",species_name)) %>%
                                        dplyr::group_by(species_name, institution_storing) %>%
                                        dplyr::summarise()

                                if(length(unique(unique.bin$species_name)) == 1 &&
                                   length(unique(unique.bin$institution_storing)) > 1 ){
                                        data.frame(Grades = "A",
                                                   Observations = paste("There were ", matches ,
                                                                        " matches. External congruence.", sep = ""))
                                }
                                else if(length(unique(unique.bin$species_name)) == 1 &&
                                         length(unique(unique.bin$institution_storing)) == 1){
                                        data.frame(Grades = "B",
                                                   Observations = paste("There were ",
                                                                        matches ,
                                                                        " matches. Internal congruence.", sep = ""))
                                }
                                else{
                                        data.frame(Grades = "E*",
                                                   Observations = paste("There were ", matches,
                                                                        " matches. ", paste(unique(unique.bin$species_name),
                                                                                            collapse = ","),
                                                                        " shared the same BIN.",
                                                                        sep = ""))
                                }
                        }
                }

        })
        return(do.call('rbind', frames))
}
