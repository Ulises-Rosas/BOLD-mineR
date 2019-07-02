library(ape)
library(XML)

ID_engine<- function(query, db, make_blast = T, ...){
  
  # query = read.FASTA("subset.fasta")[1]
  # db    = "COX1_SPECIES"
  
  seqs <- lapply(query, function(x){
    paste(as.character.DNAbin(x), collapse = "")
    })

  lapply(names(seqs), function(y){
    
    x = seqs[[y]]
    
    data <- xmlParse( paste("http://www.boldsystems.org/index.php/Ids_xml?db=",
                            db, "&sequence=", x, sep = ""))
    
    bold.results = xmlToDataFrame(data)
    
    if( nrow(bold.results) == 0 && make_blast ){
      
      rid <- RCurl::getURL(
        paste("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&QUERY=",
              x,"&WORD_SIZE=28&HITLIST_SIZE=3", sep = ""))  %>%
        gsub(".*\"RID\" value=\"", "", x =.) %>%
        gsub("\".*","", x =.) 
      
      if(rid == "")
        return(
          data.frame(ID = "GenBank: RID not available",
                          taxonomicidentification = "RID not available.",
                          similarity = 0))
      
      hits = list()
      
      while( all(is.na(hits)) ){
        
        tmp.output = RCurl::getURL(
          paste("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID=",
                rid, sep = "")) %>%
          strsplit(x = ., split = "<Hit>\n")
        hits = tmp.output[[1]][2:4]
        Sys.sleep(5)
      }
      
      genbank.results = lapply(hits, function(e){
        
        data.frame( ID = gsub(".*<Hit_id>", "", e) %>%
                      gsub("</Hit_id>\n.*", "", x =.) %>%
                      gsub("\\|", "", x =.) %>%
                      strsplit(x =., split = "gb") %>%
                      tail(x =.[[1]], n = 1) %>%
                      paste("GenBank: ", . ,sep = ""),
                    
                    taxonomicidentification = gsub(".*<Hit_def>", "", e) %>%
                      gsub("</Hit_def>.*", "", x =.) %>%
                      strsplit(x =. , split = " ") %>%
                      head(x =.[[1]], n = 2) %>%
                      paste(., collapse = " "),
                    
                    similarity = round(
                      gsub(".*<Hsp_identity>", "", e) %>%
                        gsub("</Hsp_identity>.*", "", x =.) %>%
                        as.numeric(.) / gsub(".*<Hsp_align-len>", "", e) %>%
                        gsub("</Hsp_align-len>.*", "", x =.) %>%
                        as.numeric(.) ,
                      digits = 4)
                    )
        })
      
      return(do.call("rbind", genbank.results))
      
    }else if(  nrow(bold.results) == 0 && !make_blast ){
      return(data.frame(ID = y,
                        taxonomicidentification = "Unavailable with BOLD",
                        similarity = 0))
      
    }else{
      return(bold.results)
      }
  })
}
