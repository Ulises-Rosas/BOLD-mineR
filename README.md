# BOLD-mineR

## Renombrado de cromatogramas

Se renombrarán los archivos en función a una metadata guardada en `metadata.txt`. Esta metadata se carga en la consola:

```R
meta = read.table('metadata.txt', header = T, sep = "\t")
meta2 = data.frame(cod = meta$Código, pocillo = meta$poccillo, stringsAsFactors = F)
head(meta2)
```
```
cod pocillo
1 TP_001      A1
2 TP_002      A2
3 TP_003      A3
4 TP_004      A4
5 TP_005      A5
6 TP_006      A6
```
Luego, usamos esta información para renombrar los cromatogramas Forward:

```R

stringF = list.files(pattern = "M13F-21.ab1") ##Obtenemos los nombres de los archivos que contengan ese patrón

newnamesF = vector("character") ##vector donde se guardará los nuevos nombres
oldnamesF = vector("character") ##vector donde se guardará los nombres orginales por match dentro del loop

for(i in 1:length(meta2$pocillo)){ ##loop
        pocillo = paste(as.character(meta2$pocillo), "_", sep = "") ##modficación de la columna 'pocillos'
        newnamesF[i] <- gsub(   ##función que subtituye globalmente
                paste("FISH_P1-", pocillo[i], sep = ""),
                paste(as.character(meta2$cod)[i], "_", sep = ""),
                grep(pocillo[i],stringF, value = T) ##hace el emparejamiento entre el pocillo y cromatogramas
        )
        oldnamesF[i] <- grep(pocillo[i], stringF, value = T)
}
```
y similarmente para los cromatogramas Reverse:

```R
stringR = list.files(pattern = "M13R-27.ab1") ##Obtenemos los nombres de los archivos que contengan ese patrón

newnamesR = vector("character") ##vector donde se guardará los nuevos nombres
oldnamesR = vector("character") ##vector donde se guardará los nombres orginales por match dentro del loop

for(i in 1:length(meta2$pocillo)){##loop
        pocillo = paste(as.character(meta2$pocillo), "_", sep = "") ##modficación de la columna 'pocillos'
        newnamesR[i] <- gsub( ##función que subtituye globalmente
                paste("FISH_P1-", pocillo[i], sep = ""),
                paste(as.character(meta2$cod)[i], "_", sep = ""),
                grep(pocillo[i],stringR, value = T) ##hace el emparejamiento entre el pocillo y cromatogramas
        )
        oldnamesR[i] <- grep(pocillo[i], stringR, value = T)
}
```
Si en ambas corridas aparece el siguiente mensaje:

```
Error in newnamesF[i] <- gsub(paste("FISH_P1-", pocillo[i], sep = ""),  : 
  replacement has length zero
```
Es porque no todos los nombres del vector `meta2$pocillo` han encontrado un emperajamiento dentro del loop. 

***

## SpecimenData 

La función **SpecimenData** nos permite minar la metadata asociada de cualquier espécimen según los argumentos `taxon` (e.g. Aves|Elasmobranchii), `ids` (e.g. ANGBF12704-15), `bin` (e.g. BOLD:AAA4689), `container` (e.g. FIPP), `institution` (e.g. Smithsonian Institution), `researchers` (incluye identificadores y colectores), `geo` (e.g. Peru).

La función es la siguiente:
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
Si queremos tener la información de especímenes de todos los elasmobranquios en el Peru depositados en BOLD, usamos el siguiente código:
```R
> specimendata <- SpecimenData(taxon = "Elasmobranchii", geo = "Peru")
```
Luego, solo para evaluar las dimensiones de la tabla obtenida usamos el paquete `tibble`
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

