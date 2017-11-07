library(ape)
library(XML)

ID_engine <- function(query, db, ...){

        seqs <- lapply(query, function(x){
                paste(as.character.DNAbin(x), collapse = "")})

        lapply(seqs, function(x){
                data <- xmlParse( paste("http://www.boldsystems.org/index.php/Ids_xml?db=",
                                        db, "&sequence=", x, sep = ""))
                xmlToDataFrame(data)
        })

}