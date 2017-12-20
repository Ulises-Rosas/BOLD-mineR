library(dplyr)

AuditionBarcodes<- function(species, matches){
  frames = lapply(species, function(x){
    meta.by.barcodes1 = SpecimenData(taxon = x) %>%
      select(processid, bin_uri, species_name) %>%
      mutate_if(is.factor, as.character) %>%
      filter(grepl("BOLD", bin_uri))
    if(nrow(meta.by.barcodes1) <= 3){
      data.frame(Grades = "D",
                 Observations = paste("There were ",
                                      matches ,
                                      " matches. Insufficient data.",
                                      sep = ""))
    }else{
      if(length(unique(meta.by.barcodes1$bin_uri)) > 1){
        bin = lapply(unique(meta.by.barcodes1$bin_uri), function(x){
          SpecimenData(bin = x) %>%
            select(species_name) %>%
            filter(grepl(" ",species_name), !grepl(" sp.",species_name)) %>%
            group_by(species_name) %>% summarise()
        })
        
        if(length(unique(do.call('rbind', bin)$species_name)) > 1){
          data.frame(Grades = "E**",
                     Observations = paste("There were ",
                                          matches ,
                                          " matches. More than one BIN and these are shared with different species.",
                                          sep = ""))
        }else{
          data.frame(Grades = "C*",
                     Observations = paste("There were ", matches,
                                          " matches. Assessment of intraspecific divergences is still needed.",
                                          sep = ""))
        }
      }
      else{
        unique.bin = SpecimenData(bin = unique(meta.by.barcodes1$bin_uri)) %>%
          select(species_name, institution_storing) %>%
          filter(grepl(" ",species_name), !grepl(" sp.",species_name)) %>%
          group_by(species_name, institution_storing) %>% summarise()
        
        if(length(unique(unique.bin$species_name)) == 1 &&
           length(unique(unique.bin$institution_storing)) > 1){
          data.frame(Grades = "A",
                     Observations = paste("There were ", matches ,
                                          " matches. External congruence.", sep = ""))
        }else if(length(unique(unique.bin$species_name)) == 1){
          data.frame(Grades = "B",
                     Observations = paste("There were ",
                                          matches ,
                                          " matches. Internal congruence.", sep = ""))
        }else{
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