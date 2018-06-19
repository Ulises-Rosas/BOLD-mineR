library(RCurl)
library(dplyr)

library(reticulate)
use_python("/usr/local/bin/python3") #define your python version

#py_run_string("
#import re
#import urllib
#")

# Despite you have defined modules in you python script, regrettably 
# you will still need to call them via import when using reticulate package:

urllib <- reticulate::import("urllib", convert = F)
#futhermore, if there is a module inside a directory, you also must to define it 
# owing to it seems reticulate packages can deal with it directly
urllib$request

re <- reticulate::import("re", convert = F)

source_python("worms.py", convert = F)



AuditionBarcodes <- function(species){ ##function for only using with public data

  frames = lapply(species, function(x){
    
    meta.by.barcodes1 = SpecimenData(taxon  = x) %>%
      dplyr::select(processid, bin_uri, species_name, institution_storing) %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::filter(
        grepl("BOLD", bin_uri),
        !grepl("Mined from GenBank, NCBI", institution_storing),
        !grepl("*unvouchered", institution_storing)
      )
    
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
      
      ##species and their number of records by bin:
      bin = lapply(unique(meta.by.barcodes1$bin_uri), function(x){
        #x = "BOLD:ACE4593"
        SpecimenData(bin = x) %>%
          dplyr::select(species_name, institution_storing)  %>%
          dplyr::filter(
            grepl("[A-Z][a-z]+ [a-z]+$",species_name), #just considering species level
            !grepl("Mined from GenBank, NCBI", institution_storing),
            !grepl("*unvouchered", institution_storing),
            !grepl("[A-Z][a-z]+ sp[p|\\.]{0,2}$",species_name) #just considering species level
          ) %>%
          dplyr::group_by(species_name) %>%
          dplyr::summarise(institutes = length(unique(institution_storing)),
                           n = length(species_name))%>%
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
      
      # upon having accepted names into table, assess possible synonyms within
      # elements of the list bin
      bin = lapply(bin, function(x){
        #it assumes that bold has species names correctly written
        
        #validated_names contains the match between species names of each element
        # of the list bin and 'table'. It is ordenated according to position of
        # species name on each element of the list.
        validated_names = as.character(table[match(x$species_name, names(table))])
        data.frame(species_name = validated_names,
                   x[,2:4]) %>%
          dplyr::group_by(species_name, bin) %>%
          dplyr::summarise(n = sum(n),
                           institutes = sum(institutes))
      })
      
      # this new assignment of bin is about species number contained on list's nodes.
      # since it is ordened by their lenghts, merging status of bin would appear first
      bin = sapply(bin, function(x){length(x$species_name)}) %>%
        sort(., decreasing = T) %>%
        names(.) %>%
        bin[.]
      
      
      if(length(unique(meta.by.barcodes1$bin_uri)) > 1){
        
        if(length(unique(do.call('rbind', bin)$species_name)) > 1){
          
          # a compressed way to show bin's composition is the json format
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
          
          data.frame(Grades = "E**",
                     Observations = paste("Mixtured BIN and its composition in a json-like table is: ",
                                          bin_information_json,
                                          sep = ""))
        }
        else{
          data.frame(Grades = "C",
                     Observations = paste("Splited BIN and the assessment of intraspecific divergences is still needed.",
                                          bin_information_json,
                                          sep = ""))
        }
        
      }
      else{
        if(length(unique(bin[[1]]$species_name)) == 1 &&
           sum(bin[[1]]$institutes) > 1 ){
          data.frame(Grades = "A",
                     Observations ="Matched BIN with external congruence.")
        }
        else if(length(unique(bin[[1]]$species_name)) == 1 &&
                sum(bin[[1]]$institutes) == 1 ){
          data.frame(Grades = "B",
                     Observations = "Matched BIN with internal congruence only.")
        }
        else{
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
          
          data.frame(Grades = "E*",
                     Observations = paste("Merged BIN and its composition in a json-like table is:",
                                          bin_information_json)
                     )
        }
      }
    }
    
  })
  return(do.call('rbind', frames))
}
