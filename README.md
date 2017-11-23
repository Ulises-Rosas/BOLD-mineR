# BOLD-mineR



## Renamed chromatograms

We can take advantage of names stored in a metadata to rename by default names from a Sanger sequencing (i.e. \*.ab1 files). Therein we upload `metadata.txt` in the R console: 
```R
meta = read.table('metadata.txt', header = T, sep = "\t")
meta2 = data.frame(cod = meta$Code, pocillo = meta$well, stringsAsFactors = F)
head(meta2)
```
```
cod well
1 TP_001      A1
2 TP_002      A2
3 TP_003      A3
4 TP_004      A4
5 TP_005      A5
6 TP_006      A6
```
Then, we use this information to rename Forward chromatograms:
```R

stringF = list.files(pattern = "M13F-21.ab1") ## find a forward-based pattern in my working directory
newnamesF = vector("character") ##vector where new names will be stored
oldnamesF = vector("character") ##vector where original names will be stored while it runs the loop

for(i in 1:length(meta2$well)){ ##loop
        ##we modify the 'well' column to ensure matches
        pocillo = paste(as.character(meta2$well), "_", sep = "") 
        newnamesF[i] <- gsub(   
                paste("FISH_P1-", pocillo[i], sep = ""),
                paste(as.character(meta2$cod)[i], "_", sep = ""),
                grep(pocillo[i],stringF, value = T) ##find for matches between well and chromatograms' names
        )
        oldnamesF[i] <- grep(pocillo[i], stringF, value = T)
}

file.rename(from = oldnamesF, to = newnamesF)
```
and likewise for Reverse chromatograms:
```R
stringR = list.files(pattern = "M13R-27.ab1") ## find a forward-based pattern in my working directory

newnamesR = vector("character") ##vector where new names will be stored
oldnamesR = vector("character") ##vector where original names will be stored while it runs the loop

for(i in 1:length(meta2$pocillo)){##loop
        ##we modify the 'well' column to ensure matches
        pocillo = paste(as.character(meta2$pocillo), "_", sep = "") 
        newnamesR[i] <- gsub( 
                paste("FISH_P1-", pocillo[i], sep = ""),
                paste(as.character(meta2$cod)[i], "_", sep = ""),
                grep(pocillo[i],stringR, value = T)##find for matches between well and chromatograms' names
        )
        oldnamesR[i] <- grep(pocillo[i], stringR, value = T)
}

file.rename(from = oldnamesR, to = newnamesR)
```
if on both runs suddenly appear a message like this:
```
Error in newnamesF[i] <- gsub(paste("FISH_P1-", pocillo[i], sep = ""),  : 
  replacement has length zero
```
There was not a fully correspondence between all names from `meta2$well` vector and chromatograms available in the current working directory.

## SpecimenData 

This function **SpecimenData** let us mine associated metadata from any specimen according to following arguments:

* `taxon` (e.g. Aves|Elasmobranchii).
* `ids` (e.g. ANGBF12704-15).
* `bin` (e.g. BOLD:AAA4689).
* `container` (e.g. FIPP).
* `institution` (e.g. Smithsonian Institution).
* `researchers` (including identifiers and collectors).
* `geo` (e.g. Peru).

The function is the following:

```R
library(RCurl)
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
```
If we want to get information, for instance, from all specimens of elasmobranchs distributed in Peru and stored in BOLD, we can use the following line:
```R
specimendata <- SpecimenData(taxon = "Elasmobranchii", geo = "Peru")
```
Then, we use _tibble_ package to only assess its dimension:

```R
tibble::as.tibble(specimendata)
```
```
# A tibble: 10 x 68
       processid sampleid recordID catalognum fieldnum      institution_storing collection_code      bin_uri phylum_taxID phylum_name class_taxID
 *        <fctr>   <fctr>    <int>     <fctr>   <fctr>                   <fctr>           <lgl>       <fctr>        <int>      <fctr>       <int>
 1 ANGBF10913-15 KJ146022  5651960            KJ146022 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 2 ANGBF10914-15 KJ146023  5651961            KJ146023 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 3 ANGBF10915-15 KJ146024  5651962            KJ146024 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 4 ANGBF10916-15 KJ146025  5651963            KJ146025 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 5 ANGBF10917-15 KJ146026  5651964            KJ146026 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 6 ANGBF11015-15 KJ146027  5652062            KJ146027 Mined from GenBank, NCBI              NA BOLD:ACE6212           18    Chordata       34196
 7 ANGBF11043-15 KJ146021  5652090            KJ146021 Mined from GenBank, NCBI              NA BOLD:ABZ0850           18    Chordata       34196
 8 ANGBF11064-15 KJ146044  5652111            KJ146044 Mined from GenBank, NCBI              NA BOLD:AAA7096           18    Chordata       34196
 9 ANGBF11065-15 KJ146043  5652112            KJ146043 Mined from GenBank, NCBI              NA BOLD:AAA7096           18    Chordata       34196
10 ANGBF11066-15 KJ146042  5652113            KJ146042 Mined from GenBank, NCBI              NA BOLD:AAA7096           18    Chordata       34196
# ... with 57 more variables: class_name <fctr>, order_taxID <int>, order_name <fctr>, family_taxID <int>, family_name <fctr>, subfamily_taxID <lgl>,
#   subfamily_name <lgl>, genus_taxID <int>, genus_name <fctr>, species_taxID <int>, species_name <fctr>, subspecies_taxID <lgl>, subspecies_name <lgl>,
#   identification_provided_by <fctr>, identification_method <lgl>, identification_reference <fctr>, tax_note <lgl>, voucher_status <lgl>,
#   tissue_type <lgl>, collection_event_id <lgl>, collectors <fctr>, collectiondate_start <lgl>, collectiondate_end <lgl>, collectiontime <lgl>,
#   collection_note <lgl>, site_code <lgl>, sampling_protocol <fctr>, lifestage <lgl>, sex <lgl>, reproduction <lgl>, habitat <fctr>,
#   associated_specimens <lgl>, associated_taxa <lgl>, extrainfo <fctr>, notes <fctr>, lat <dbl>, lon <dbl>, coord_source <lgl>, coord_accuracy <lgl>,
#   elev <lgl>, depth <int>, elev_accuracy <lgl>, depth_accuracy <lgl>, country <fctr>, province_state <fctr>, region <fctr>, sector <fctr>,
#   exactsite <lgl>, image_ids <lgl>, image_urls <lgl>, media_descriptors <lgl>, captions <lgl>, copyright_holders <lgl>, copyright_years <lgl>,
#   copyright_licenses <lgl>, copyright_institutions <lgl>, photographers <lgl>
```
We can also add sequences, and its associated information, from each specimen by using the argument `seq = "combined"`:

```R
tibble::as.tibble(SpecimenData(taxon = "Elasmobranchii", geo = "Peru", seq = "combined"))
```

```
# A tibble: 47 x 80
       processid sampleid recordID catalognum fieldnum      institution_storing collection_code      bin_uri phylum_taxID phylum_name class_taxID
          <fctr>   <fctr>    <int>     <fctr>   <fctr>                   <fctr>           <lgl>       <fctr>        <int>      <fctr>       <int>
 1 ANGBF10913-15 KJ146022  5651960            KJ146022 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 2 ANGBF10914-15 KJ146023  5651961            KJ146023 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 3 ANGBF10915-15 KJ146024  5651962            KJ146024 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 4 ANGBF10916-15 KJ146025  5651963            KJ146025 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 5 ANGBF10917-15 KJ146026  5651964            KJ146026 Mined from GenBank, NCBI              NA BOLD:AAB0473           18    Chordata       34196
 6 ANGBF11015-15 KJ146027  5652062            KJ146027 Mined from GenBank, NCBI              NA BOLD:ACE6212           18    Chordata       34196
 7 ANGBF11043-15 KJ146021  5652090            KJ146021 Mined from GenBank, NCBI              NA BOLD:ABZ0850           18    Chordata       34196
 8 ANGBF11064-15 KJ146044  5652111            KJ146044 Mined from GenBank, NCBI              NA BOLD:AAA7096           18    Chordata       34196
 9 ANGBF11065-15 KJ146043  5652112            KJ146043 Mined from GenBank, NCBI              NA BOLD:AAA7096           18    Chordata       34196
10 ANGBF11066-15 KJ146042  5652113            KJ146042 Mined from GenBank, NCBI              NA BOLD:AAA7096           18    Chordata       34196
# ... with 37 more rows, and 69 more variables: class_name <fctr>, order_taxID <int>, order_name <fctr>, family_taxID <int>, family_name <fctr>,
#   subfamily_taxID <lgl>, subfamily_name <lgl>, genus_taxID <int>, genus_name <fctr>, species_taxID <int>, species_name <fctr>, subspecies_taxID <lgl>,
#   subspecies_name <lgl>, identification_provided_by <fctr>, identification_method <lgl>, identification_reference <fctr>, tax_note <lgl>,
#   voucher_status <lgl>, tissue_type <lgl>, collection_event_id <lgl>, collectors <fctr>, collectiondate_start <lgl>, collectiondate_end <lgl>,
#   collectiontime <lgl>, collection_note <lgl>, site_code <lgl>, sampling_protocol <fctr>, lifestage <lgl>, sex <lgl>, reproduction <lgl>,
#   habitat <fctr>, associated_specimens <lgl>, associated_taxa <lgl>, extrainfo <fctr>, notes <fctr>, lat <dbl>, lon <dbl>, coord_source <lgl>,
#   coord_accuracy <lgl>, elev <lgl>, depth <int>, elev_accuracy <lgl>, depth_accuracy <lgl>, country <fctr>, province_state <fctr>, region <fctr>,
#   sector <fctr>, exactsite <lgl>, image_ids <lgl>, image_urls <lgl>, media_descriptors <lgl>, captions <lgl>, copyright_holders <lgl>,
#   copyright_years <lgl>, copyright_licenses <lgl>, copyright_institutions <lgl>, photographers <lgl>, sequenceID <int>, markercode <fctr>,
#   genbank_accession <fctr>, nucleotides <fctr>, trace_ids <fctr>, trace_names <fctr>, trace_links <fctr>, run_dates <fctr>, sequencing_centers <fctr>,
#   directions <fctr>, seq_primers <fctr>, marker_codes <fctr>
```

However, if only sequences are desired, the argument `seq = "combined"` should change to `seq = "only"`:
```R
SpecimenData(taxon = "Elasmobranchii", geo = "Peru", seq = "only")
```
```
47 DNA sequences in binary format stored in a list.

Mean sequence length: 658.851 
   Shortest sequence: 586 
    Longest sequence: 712 

Labels:
ANGBF10913-15|Alopias pelagicus|COI-5P|KJ146022
ANGBF10914-15|Alopias pelagicus|COI-5P|KJ146023
ANGBF10915-15|Alopias pelagicus|COI-5P|KJ146024
ANGBF10916-15|Alopias pelagicus|COI-5P|KJ146025
ANGBF10917-15|Alopias pelagicus|COI-5P|KJ146026
ANGBF11015-15|Carcharhinus brachyurus|COI-5P|KJ146027
...

Base composition:
    a     c     g     t 
0.256 0.264 0.166 0.314 
```
Above sequence can be also stored by using `write.dna()` function from the _ape_ package:
```R
seqs <- SpecimenData(taxon = "Elasmobranchii", geo = "Peru", seqs = "only")[1:5] ##we only select five sequences
write.dna(seqs, 'secuencias.txt', format = 'fasta', nbcol = 1, colw = 90)
```
## ID_engine

**ID_engine** finds best matches between a query sequence and a database of BOLD by using BLASTn-based algorithms. Arguments of this function are `query` and `db`. First one are query sequences and the second one are one of avilable databases in BOLD (i.e. `COX1`, `COX1_SPECIES`, `COX1_SPECIES_PUBLIC` and `COX1_L640bp`).

The function is the following:
```R
library(ape)
library(XML)

ID_engine<- function(query, db, ...){
        seqs <- lapply(query, function(x){
                paste(as.character.DNAbin(x), collapse = "")})
        
        lapply(seqs, function(x){
                data <- xmlParse( paste("http://www.boldsystems.org/index.php/Ids_xml?db=",
                                        db, "&sequence=", x, sep = ""))
                xmlToDataFrame(data)
        })
}
```

If we want to identify at species level, for instance, all sequence samples stored in [secuencias.txt](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/secuencias.txt), we should run forthcoming lines:

```R
out <- ID_engine(query = read.FASTA('secuencias.txt'), db = "COX1_SPECIES")

```
First seven rows and first, fifth and sixth column of each element from the result list `out` are shown:

```R
lapply(out, function(x){
        x[1:7, c(1,5,6)]
})
```
```
$`ANGBF10913-15|Alopias pelagicus|COI-5P|KJ146022`
             ID taxonomicidentification similarity
1   PHANT458-08             Lamniformes          1
2   IRREK872-08             Lamniformes          1
3   IRREK873-08             Lamniformes          1
4   IRREK874-08             Lamniformes          1
5   IRREK876-08             Lamniformes          1
6 ANGBF12625-15       Alopias pelagicus          1
7 ANGBF12624-15       Alopias pelagicus          1

$`ANGBF10914-15|Alopias pelagicus|COI-5P|KJ146023`
             ID taxonomicidentification similarity
1 ANGBF10914-15       Alopias pelagicus          1
2   ESHKB029-07       Alopias pelagicus          1
3 ANGBF10913-15       Alopias pelagicus     0.9985
4   PHANT458-08             Lamniformes     0.9985
5   IRREK873-08             Lamniformes     0.9985
6   IRREK874-08             Lamniformes     0.9985
7   IRREK876-08             Lamniformes     0.9985

$`ANGBF10915-15|Alopias pelagicus|COI-5P|KJ146024`
             ID taxonomicidentification similarity
1 ANGBF12626-15       Alopias pelagicus          1
2 ANGBF12623-15       Alopias pelagicus          1
3 ANGBF11723-15       Alopias pelagicus          1
4 ANGBF10915-15       Alopias pelagicus          1
5   ESHKB036-07       Alopias pelagicus          1
6   ESHKB030-07       Alopias pelagicus          1
7   ESHKB031-07       Alopias pelagicus          1

$`ANGBF10916-15|Alopias pelagicus|COI-5P|KJ146025`
             ID taxonomicidentification similarity
1 ANGBF11731-15       Alopias pelagicus          1
2 ANGBF10916-15       Alopias pelagicus          1
3 ANGBF10915-15       Alopias pelagicus     0.9986
4   ESHKB031-07       Alopias pelagicus     0.9985
5   ESHKB034-07       Alopias pelagicus     0.9985
6 ANGBF11723-15       Alopias pelagicus     0.9984
7   ESHKB030-07       Alopias pelagicus     0.9984

$`ANGBF10917-15|Alopias pelagicus|COI-5P|KJ146026`
             ID taxonomicidentification similarity
1 ANGBF11729-15       Alopias pelagicus          1
2 ANGBF10917-15       Alopias pelagicus          1
3   ESHKB036-07       Alopias pelagicus          1
4 ANGBF10915-15       Alopias pelagicus     0.9986
5   ESHKB031-07       Alopias pelagicus     0.9985
6   ESHKB034-07       Alopias pelagicus     0.9985
7 ANGBF11723-15       Alopias pelagicus     0.9984
```
## addAudition

This function adds an **audition** step (Oliveira _et al._ 2016) to each specimen selected in the `ID_engine()` given a certain threshold. This function, in turn, uses another function called `AuditionBarcodes()`: 
```R
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
                                                                        " matches. ", paste(unique(bins$species_name),
                                                                                            collapse = ","),
                                                                        " shared the same BIN.",
                                                                        sep = ""))
                                }
                        }
                }
        })
        return(do.call('rbind', frames))
}
```
Finally,  `addAudition` is composed by these lines:

```R
library(dplyr)

addAudition <- function(seqs, threshold){
        lista2 = list()
        pb <- txtProgressBar(min = 0, max = length(seqs), style = 3, char = "*")
        for(i in 1:length(seqs)){
                Sys.sleep(0.0001)
                tmp = ID_engine(seqs[i], db = "COX1_SPECIES")
                tmp = tmp[[1]] %>%
                        select(ID, taxonomicidentification, similarity) %>%
                        filter(grepl(" ", taxonomicidentification), !grepl(" sp.",taxonomicidentification)) %>%
                        mutate(similarity = as.numeric(as.character(similarity)))
                if(tmp[1,]$similarity < threshold){
                        lista2[[i]] = data.frame(Match = "Any",
                                                 Species = "",
                                                 Grades = "",
                                                 Observations = paste("Three best matches: ",
                                                                      tmp$taxonomicidentification[1], " (", tmp$similarity[1],"), ",
                                                                      tmp$taxonomicidentification[2], " (", tmp$similarity[2],"), ",
                                                                      tmp$taxonomicidentification[3], " (", tmp$similarity[3],").", sep = ""))
                }else{
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
                                                                                               matches = sum(barcodes))$Grades,
                                                                              collapse = ", ")," respectively.",sep = ""),
                                                         Observations = "")
                        }else{
                                lista2[[i]] = data.frame(Match = "Unique",
                                                         Species = paste(names(barcodes)),
                                                         AuditionBarcodes(species = names(barcodes),
                                                                          matches = sum(barcodes)))
                        }
                        setTxtProgressBar(pb, i)
                }
        }
        close(pb)
        return(data.frame(Samples = names(seqs), do.call('rbind', lista2)))
}
```
In order to test its efficiency, species-level identification of samples stored in [secuencias.txt](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/secuencias.txt) using a `threshold = 0.99` is conducted:

```R
addAudition(seqs = read.FASTA('secuencias.txt'), threshold = 0.99)
#  |**********************************************************************************************************************************************| 100%
#                                          Samples  Match           Species Grades                                Observations
#1 ANGBF10913-15|Alopias pelagicus|COI-5P|KJ146022 Unique Alopias pelagicus      A There were 44 matches. External congruence.
#2 ANGBF10914-15|Alopias pelagicus|COI-5P|KJ146023 Unique Alopias pelagicus      A There were 34 matches. External congruence.
#3 ANGBF10915-15|Alopias pelagicus|COI-5P|KJ146024 Unique Alopias pelagicus      A There were 70 matches. External congruence.
#4 ANGBF10916-15|Alopias pelagicus|COI-5P|KJ146025 Unique Alopias pelagicus      A There were 34 matches. External congruence.
#5 ANGBF10917-15|Alopias pelagicus|COI-5P|KJ146026 Unique Alopias pelagicus      A There were 40 matches. External congruence.
```
