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
        pocillo = paste(as.character(meta2$pocillo), "_", sep = "") ##modficación de la columna 'pocillos' para que corra en el loop
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
        pocillo = paste(as.character(meta2$pocillo), "_", sep = "") ##modficación de la columna 'pocillos' para que corra en el loop
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
````
>SpecimenData(taxon = "Elasmobranchii", geo = "Peru")
        processid      sampleid recordID    catalognum     fieldnum                          institution_storing collection_code
 1: ANGBF10913-15      KJ146022  5651960                   KJ146022                     Mined from GenBank, NCBI              NA
 2: ANGBF10914-15      KJ146023  5651961                   KJ146023                     Mined from GenBank, NCBI              NA
 3: ANGBF10915-15      KJ146024  5651962                   KJ146024                     Mined from GenBank, NCBI              NA
 4: ANGBF10916-15      KJ146025  5651963                   KJ146025                     Mined from GenBank, NCBI              NA
 5: ANGBF10917-15      KJ146026  5651964                   KJ146026                     Mined from GenBank, NCBI              NA
 6: ANGBF11015-15      KJ146027  5652062                   KJ146027                     Mined from GenBank, NCBI              NA
 7: ANGBF11043-15      KJ146021  5652090                   KJ146021                     Mined from GenBank, NCBI              NA
 8: ANGBF11064-15      KJ146044  5652111                   KJ146044                     Mined from GenBank, NCBI              NA
 9: ANGBF11065-15      KJ146043  5652112                   KJ146043                     Mined from GenBank, NCBI              NA
10: ANGBF11066-15      KJ146042  5652113                   KJ146042                     Mined from GenBank, NCBI              NA
11: ANGBF11633-15      KJ146045  5652680                   KJ146045                     Mined from GenBank, NCBI              NA
12: ANGBF11634-15      KJ146041  5652681                   KJ146041                     Mined from GenBank, NCBI              NA
13: ANGBF11635-15      KJ146040  5652682                   KJ146040                     Mined from GenBank, NCBI              NA
14: ANGBF11636-15      KJ146039  5652683                   KJ146039                     Mined from GenBank, NCBI              NA
         bin_uri phylum_taxID phylum_name class_taxID     class_name order_taxID        order_name family_taxID    family_name
 1: BOLD:AAB0473           18    Chordata       34196 Elasmobranchii         225       Lamniformes         1201      Alopiidae
 2: BOLD:AAB0473           18    Chordata       34196 Elasmobranchii         225       Lamniformes         1201      Alopiidae
 3: BOLD:AAB0473           18    Chordata       34196 Elasmobranchii         225       Lamniformes         1201      Alopiidae
 4: BOLD:AAB0473           18    Chordata       34196 Elasmobranchii         225       Lamniformes         1201      Alopiidae
 5: BOLD:AAB0473           18    Chordata       34196 Elasmobranchii         225       Lamniformes         1201      Alopiidae
 6: BOLD:ACE6212           18    Chordata       34196 Elasmobranchii         214 Carcharhiniformes          944 Carcharhinidae
 7: BOLD:ABZ0850           18    Chordata       34196 Elasmobranchii         214 Carcharhiniformes          944 Carcharhinidae
 8: BOLD:AAA7096           18    Chordata       34196 Elasmobranchii         214 Carcharhiniformes          944 Carcharhinidae
 9: BOLD:AAA7096           18    Chordata       34196 Elasmobranchii         214 Carcharhiniformes          944 Carcharhinidae
10: BOLD:AAA7096           18    Chordata       34196 Elasmobranchii         214 Carcharhiniformes          944 Carcharhinidae
11: BOLD:AAA6753           18    Chordata       34196 Elasmobranchii         214 Carcharhiniformes          938     Sphyrnidae
12: BOLD:AAA3577           18    Chordata       34196 Elasmobranchii         225       Lamniformes         1214       Lamnidae
13: BOLD:AAA3577           18    Chordata       34196 Elasmobranchii         225       Lamniformes         1214       Lamnidae
14: BOLD:AAA3577           18    Chordata       34196 Elasmobranchii         225       Lamniformes         1214       Lamnidae
    subfamily_taxID subfamily_name genus_taxID   genus_name species_taxID            species_name subspecies_taxID subspecies_name
 1:              NA             NA        3007      Alopias         72350       Alopias pelagicus               NA              NA
 2:              NA             NA        3007      Alopias         72350       Alopias pelagicus               NA              NA
 3:              NA             NA        3007      Alopias         72350       Alopias pelagicus               NA              NA
 4:              NA             NA        3007      Alopias         72350       Alopias pelagicus               NA              NA
 5:              NA             NA        3007      Alopias         72350       Alopias pelagicus               NA              NA
 6:              NA             NA        2960 Carcharhinus         12779 Carcharhinus brachyurus               NA              NA
 7:              NA             NA        2960 Carcharhinus         12769   Carcharhinus obscurus               NA              NA
 8:              NA             NA        2954     Prionace         12783         Prionace glauca               NA              NA
 9:              NA             NA        2954     Prionace         12783         Prionace glauca               NA              NA
10:              NA             NA        2954     Prionace         12783         Prionace glauca               NA              NA
11:              NA             NA        3619      Sphyrna         16977         Sphyrna zygaena               NA              NA
12:              NA             NA       59203        Lamna         89093             Lamna nasus               NA              NA
13:              NA             NA       59203        Lamna         89093             Lamna nasus               NA              NA
14:              NA             NA       59203        Lamna         89093             Lamna nasus               NA              NA
    identification_provided_by identification_method identification_reference tax_note voucher_status tissue_type
 1:                                               NA           Nakamura, 1935       NA             NA          NA
 2:                                               NA           Nakamura, 1935       NA             NA          NA
 3:                                               NA           Nakamura, 1935       NA             NA          NA
 4:                                               NA           Nakamura, 1935       NA             NA          NA
 5:                                               NA           Nakamura, 1935       NA             NA          NA
 6:                                               NA         (GÃ¼nther, 1870)       NA             NA          NA
 7:                                               NA          (Lesueur, 1818)       NA             NA          NA
 8:                                               NA         (Linnaeus, 1758)       NA             NA          NA
 9:                                               NA         (Linnaeus, 1758)       NA             NA          NA
10:                                               NA         (Linnaeus, 1758)       NA             NA          NA
11:                                               NA         (Linnaeus, 1758)       NA             NA          NA
12:                                               NA       (Bonnaterre, 1788)       NA             NA          NA
13:                                               NA       (Bonnaterre, 1788)       NA             NA          NA
14:                                               NA       (Bonnaterre, 1788)       NA             NA          NA
    collection_event_id                                                      collectors collectiondate_start collectiondate_end
 1:                  NA                                                    ProDelphinus                   NA                 NA
 2:                  NA                                                    ProDelphinus                   NA                 NA
 3:                  NA                                                    ProDelphinus                   NA                 NA
 4:                  NA                                                    ProDelphinus                   NA                 NA
 5:                  NA                                                    ProDelphinus                   NA                 NA
 6:                  NA                                                    ProDelphinus                   NA                 NA
 7:                  NA                                                    ProDelphinus                   NA                 NA
 8:                  NA                                                    ProDelphinus                   NA                 NA
 9:                  NA                                                    ProDelphinus                   NA                 NA
10:                  NA                                                    ProDelphinus                   NA                 NA
11:                  NA                                                    ProDelphinus                   NA                 NA
12:                  NA                                                    ProDelphinus                   NA                 NA
13:                  NA                                                    ProDelphinus                   NA                 NA
14:                  NA                                                    ProDelphinus                   NA                 NA
    collectiontime collection_note site_code sampling_protocol lifestage sex reproduction         habitat associated_specimens
 1:             NA              NA        NA                          NA  NA           NA                                   NA
 2:             NA              NA        NA                          NA  NA           NA                                   NA
 3:             NA              NA        NA                          NA  NA           NA                                   NA
 4:             NA              NA        NA                          NA  NA           NA                                   NA
 5:             NA              NA        NA                          NA  NA           NA                                   NA
 6:             NA              NA        NA                          NA  NA           NA                                   NA
 7:             NA              NA        NA                          NA  NA           NA                                   NA
 8:             NA              NA        NA                          NA  NA           NA                                   NA
 9:             NA              NA        NA                          NA  NA           NA                                   NA
10:             NA              NA        NA                          NA  NA           NA                                   NA
11:             NA              NA        NA                          NA  NA           NA                                   NA
12:             NA              NA        NA                          NA  NA           NA                                   NA
13:             NA              NA        NA                          NA  NA           NA                                   NA
14:             NA              NA        NA                          NA  NA           NA                                   NA
    associated_taxa    extrainfo          notes       lat      lon coord_source coord_accuracy elev depth elev_accuracy
 1:              NA  taxon:57979                       NA       NA           NA             NA   NA    NA            NA
 2:              NA  taxon:57979                       NA       NA           NA             NA   NA    NA            NA
 3:              NA  taxon:57979                       NA       NA           NA             NA   NA    NA            NA
 4:              NA  taxon:57979                       NA       NA           NA             NA   NA    NA            NA
 5:              NA  taxon:57979                       NA       NA           NA             NA   NA    NA            NA
 6:              NA taxon:671158                       NA       NA           NA             NA   NA    NA            NA
 7:              NA   taxon:7807                       NA       NA           NA             NA   NA    NA            NA
 8:              NA   taxon:7815                       NA       NA           NA             NA   NA    NA            NA
 9:              NA   taxon:7815                       NA       NA           NA             NA   NA    NA            NA
10:              NA   taxon:7815                       NA       NA           NA             NA   NA    NA            NA
11:              NA taxon:195335                       NA       NA           NA             NA   NA    NA            NA
12:              NA   taxon:7849                       NA       NA           NA             NA   NA    NA            NA
13:              NA   taxon:7849                       NA       NA           NA             NA   NA    NA            NA
14:              NA   taxon:7849                       NA       NA           NA             NA   NA    NA            NA
    depth_accuracy country province_state                    region         sector exactsite image_ids image_urls media_descriptors
 1:             NA    Peru                                                                NA        NA         NA                NA
 2:             NA    Peru                                                                NA        NA         NA                NA
 3:             NA    Peru                                                                NA        NA         NA                NA
 4:             NA    Peru                                                                NA        NA         NA                NA
 5:             NA    Peru                                                                NA        NA         NA                NA
 6:             NA    Peru                                                                NA        NA         NA                NA
 7:             NA    Peru                                                                NA        NA         NA                NA
 8:             NA    Peru                                                                NA        NA         NA                NA
 9:             NA    Peru                                                                NA        NA         NA                NA
10:             NA    Peru                                                                NA        NA         NA                NA
11:             NA    Peru                                                                NA        NA         NA                NA
12:             NA    Peru                                                                NA        NA         NA                NA
13:             NA    Peru                                                                NA        NA         NA                NA
14:             NA    Peru                                                                NA        NA         NA                NA
    captions copyright_holders copyright_years copyright_licenses copyright_institutions photographers
 1:       NA                NA              NA                 NA                     NA            NA
 2:       NA                NA              NA                 NA                     NA            NA
 3:       NA                NA              NA                 NA                     NA            NA
 4:       NA                NA              NA                 NA                     NA            NA
 5:       NA                NA              NA                 NA                     NA            NA
 6:       NA                NA              NA                 NA                     NA            NA
 7:       NA                NA              NA                 NA                     NA            NA
 8:       NA                NA              NA                 NA                     NA            NA
 9:       NA                NA              NA                 NA                     NA            NA
10:       NA                NA              NA                 NA                     NA            NA
11:       NA                NA              NA                 NA                     NA            NA
12:       NA                NA              NA                 NA                     NA            NA
13:       NA                NA              NA                 NA                     NA            NA
14:       NA                NA              NA                 NA                     NA            NA
 [ reached getOption("max.print") -- omitted 34 rows ]
```

