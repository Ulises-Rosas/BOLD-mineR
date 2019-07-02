library(dplyr)
library(RCurl)


addAudition <- function(seqs, threshold,
                        include_ncbi=F,
                        just_ID = F,
                        make_blast = F ){
  # 
  # seqs = read.FASTA("subset.fasta")
  # threshold = 0.99

  lista2 = list()
  pb <- txtProgressBar(min = 0, max = length(seqs), style = 3, char = "*")
  
  for( i in 1:length(seqs) ){
    # i = 1
    tmp = ID_engine(seqs[i], db = "COX1_SPECIES", make_blast)
    
    
    if( grepl("Unavailable", tmp[[1]]$taxonomicidentification[1]) ){
      
      if(!just_ID){
        data.frame(Match   = "Any",
                   Species = "",
                   Grades  = "",
                   Observations = tmp[[1]]$taxonomicidentification) -> lista2[[i]]
      }else{
        data.frame(Match   = "Any", Species = tmp[[1]]$taxonomicidentification) -> lista2[[i]]
      }
      setTxtProgressBar(pb, i)
      next
    }
    
    tmp = tmp[[1]] %>%
      dplyr::select(ID, taxonomicidentification, similarity ) %>%
      dplyr::filter(grepl("[A-Z][a-z]+ [a-z]+$", taxonomicidentification)) %>%
      dplyr::mutate(similarity = as.numeric(as.character(similarity)))
    
    if(!grepl("GenBank", as.character(tmp[1,]$ID) )){
      
      tmp = tmp %>%
        dplyr::filter(!grepl(" sp.",taxonomicidentification))
    }
    
    if( tmp[1,]$similarity < threshold ){
      
      paste0(
        if(grepl("GenBank", as.character(tmp[1,]$ID))) "GenBank: " else "BOLD: ",
        paste(
          head(
            apply(tmp,MARGIN = 1, function(x){paste0(x[2], " (sim. = ",x[3], ")") }),
            n = 3),
          collapse = ", ")
        ) -> best_matches
      
      if(!just_ID){
        
        if(grepl("RID", as.character(tmp[1,]$ID))){
          paste("RID not available.") -> obs
        }else{
          best_matches -> obs
        }
        
        data.frame(Match   = "Any",
                   Species = "",
                   Grades  = "",
                   Observations = obs ) -> lista2[[i]]
      }else{
        data.frame(Match = "Any",Species = best_matches) -> lista2[[i]] 
        }
      
      }else{
        
        tmp = tmp %>% 
          dplyr::filter(similarity > threshold)
        
        barcodes = sort(
          table(as.character(tmp$taxonomicidentification)),
          decreasing = T)
        
      if(length(barcodes) > 1){
        
        vec <- vector('character')
        for(k in 1:length(barcodes)){
          vec[k] = paste(names(barcodes[k])," (n = ",barcodes[k],")", sep = "")
          }
        
        if(!just_ID){
          data.frame(Match = "Ambiguous",
                     Species = paste(vec, collapse = ", "),
                     Grades = paste(
                       paste(
                         AuditionBarcodes(species = names(barcodes),
                                          matches = sum(barcodes),
                                          include_ncbi = include_ncbi)$Grades,
                         collapse = ", "),
                       " respectively.",sep = ""),
                     Observations = "") -> lista2[[i]]
          }else{
            data.frame(Match = "Ambiguous", Species = paste(vec, collapse = ", ") ) -> lista2[[i]] }
        
        }else{
          
          if(!just_ID){
            
            data.frame(Match = "Unique",
                       Species = paste(names(barcodes)),
                       AuditionBarcodes(species = names(barcodes),
                                        matches = sum(barcodes),
                                        include_ncbi = include_ncbi )) -> lista2[[i]]
          }else{
            
            data.frame(Match   = "Unique",
                       Species = paste(names(barcodes))) -> lista2[[i]]
            }
          }
        setTxtProgressBar(pb, i)
      }
    }
  close(pb)
  return(data.frame(Samples = names(seqs), do.call('rbind', lista2)))
}