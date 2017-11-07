library(ape)
library(data.table)

SpecimenData <- function(taxon, ids, bin, container,
                         institutions, researchers, geo, ...){
        input <- data.frame(
                names = names(as.list(environment())),
                args = sapply(as.list(environment()), paste),
                stringsAsFactors = F
        )

        #text <- RCurl::getURL(
        URLtxt <- paste(if(list(...)[1] == "only"){
                "http://www.boldsystems.org/index.php/API_Public/sequence?"}
                else{if(list(...)[1] == "combined"){
                        "http://www.boldsystems.org/index.php/API_Public/combined?"}
                        else{
                                "http://www.boldsystems.org/index.php/API_Public/specimen?"}
                },
                paste(
                        paste(input[!(input$args == ""),]$names,
                              "=",
                              sapply(input[!(input$args == ""),]$args, function(x){
                                      if(length(strsplit(x, " ")[[1]]) > 1){
                                              paste(gsub(" ", "%20", x), "&", sep = "")
                                      }else{paste(x, "&", sep = "")}}
                              ),
                              sep = ""),
                        collapse = ""),
                "format=tsv",
                sep = "")
        text <- RCurl::getURL(URLtxt)

        if(list(...)[1] == "only")
                return(ape::read.FASTA(textConnection(text)))

        data.table::fread(text)
}