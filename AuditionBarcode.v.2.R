library(RCurl)
library(dplyr)

library(reticulate)
use_python("/usr/local/bin/python3.6") #define your python version

py_run_string("
import re
import urllib
")


AuditionBarcodes <- function(species){ ##function for only using with public data


        frames = lapply(species, function(x){

                meta.by.barcodes1 = SpecimenData(taxon  = x) %>%
                        dplyr::select(processid, bin_uri, species_name, institution_storing) %>%
                        dplyr::mutate_if(is.factor, as.character) %>%
                        dplyr::filter(grepl("BOLD", bin_uri),
                                      !grepl("Mined from GenBank, NCBI", institution_storing),
                                      !grepl("*unvouchered", institution_storing))

                ## Total number of records and institutions storing barcodes of the query
                ## This includes either public or private data.
                js0 = getURL(
                        paste("http://www.boldsystems.org/index.php/API_Tax/TaxonSearch?taxName=",
                              gsub(" ","%20", x), sep = "")) %>%
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
                                   Observations = "Barcodes mined from GenBank or unvouchered")
                }
                else if(nrow(meta.by.barcodes1) <= 3 && sum(js0$records, na.rm = T) != 0){
                        data.frame(Grades = "D",
                                   Observations = paste("Insufficient data. Institution storing: ",
                                                        length(js0$institutions),
                                                        ". Total specimen records: ",
                                                        sum(js0$records, na.rm = T),
                                                        sep = ""))
                }
                else{
                        if(length(unique(meta.by.barcodes1$bin_uri)) > 1){

                                ##species and their number of records by bin:
                                bin = lapply(unique(meta.by.barcodes1$bin_uri), function(x){
                                        SpecimenData(bin = x) %>%
                                                dplyr::select(species_name, institution_storing)  %>%
                                                dplyr::filter(
                                                        grepl("[A-Z][a-z]+ [a-z]+$",species_name), #just considering species level
                                                        !grepl("Mined from GenBank, NCBI", institution_storing),
                                                        !grepl("*unvouchered", institution_storing),
                                                        !grepl("[A-Z][a-z]+ sp[p|\\.]{0,2}$",species_name) #just considering species level
                                                ) %>%
                                                dplyr::group_by(species_name) %>%
                                                dplyr::count() %>%
                                                mutate(bin = x)
                                })


                                names(bin) = unique(meta.by.barcodes1$bin_uri)

                                #table with accepted names per each species
                                table = sapply(unique(do.call('rbind', bin)$species_name),
                                       function(x){
                                               #it gets currently accepted names
                                               worms(x)$get_accepted_name() %>%
                                                       as.character(.)
                                               })

                                bin = lapply(bin, function(x){
                                        #it assumes that bold has species names correctly written
                                        validated_names = as.character(table[ names(table) %in% x$species_name])
                                        data.frame(species_name = validated_names,
                                                   x[,2:3]) %>%
                                                dplyr::group_by(species_name, bin) %>%
                                                dplyr::summarise(n = sum(n))
                                        }) %>%
                                        sapply(., function(x){length(x$species_name)}) %>%
                                        sort(., decreasing = T) %>%
                                        names(.) %>%
                                        bin[.]

                                bin_information_json <- lapply(bin, function(x){
                                        paste("'",
                                              unique(x$bin),
                                              "':{",
                                              paste(
                                                      paste("'",x$species_name,"':", x$n, sep = ""),
                                                      collapse = ","),
                                              "}",
                                              sep = "")}) %>%
                                        unlist(.) %>%
                                        paste(., collapse = ", ")


                                if(length(unique(do.call('rbind', bin)$species_name)) > 1){

                                        data.frame(Grades = "E**",
                                                   Observations = paste("Mixtured BIN and it's composition in a json-like table is: ",
                                                                        bin_information_json,
                                                                        sep = ""))
                                }
                                else{
                                        data.frame(Grades = "C*",
                                                   Observations = paste("Assessment of intraspecific divergences is still needed.",
                                                                        bin_information_json,
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
                                                   Observations ="External congruence.")
                                }
                                else if(length(unique(unique.bin$species_name)) == 1 &&
                                        length(unique(unique.bin$institution_storing)) == 1){
                                        data.frame(Grades = "B",
                                                   Observations = "Internal congruence.")
                                }
                                else{
                                        data.frame(Grades = "E*",
                                                   Observations = paste(paste(unique(unique.bin$species_name), collapse = ","),
                                                                        " are shared the same BIN.",
                                                                        sep = ""))
                                }
                        }
                }

        })
        return(do.call('rbind', frames))
}
