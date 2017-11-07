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

file.rename(from = oldnamesF, to = newnamesF)
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

file.rename(from = oldnamesR, to = newnamesR)
```
Si en ambas corridas aparece el siguiente mensaje:

```
Error in newnamesF[i] <- gsub(paste("FISH_P1-", pocillo[i], sep = ""),  : 
  replacement has length zero
```
Es porque no todos los nombres del vector `meta2$pocillo` han encontrado un emparejamiento dentro del loop. 


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
Si queremos tener la información, por ejemplo, de especímenes de todos los elasmobranquios en el Peru depositados en BOLD, podemos usar el siguiente código:
```R
specimendata <- SpecimenData(taxon = "Elasmobranchii", geo = "Peru")
```
Luego, solo para evaluar las dimensiones de la tabla obtenida usamos el paquete _tibble_
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
Podemos también incluir en las últimas 13 columnas información de secuencias con el argumento `seq = "combined"`:

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
Si solo se desea las secuencias de la anterior tabla, se debe modificar el argumento `seq = combined` a `seq = only`:

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
Se pueden guardar las secuencias anteriores usando la funcion `write.dna()` del paquete _ape_:
```R
seqs <- SpecimenData(taxon = "Elasmobranchii", geo = "Peru", seqs = "only")[1:5] ##seleccionamos 5 secuencias
write.dna(seqs, 'secuencias.txt', format = 'fasta', nbcol = 1, colw = 90)
```
## ID_engine

La función **ID_engine** nos permite encontrar los especímenes del repositorio BOLD que generan los mejores porcentajes de similitud dada una secuencia cualquiera a través de algoritmos inspirados en BLASTn. Los argumentos de esta función son `query` y `db`. El primer argumento son las secuencias problema o consulta y el segundo argumento es una de las cuatro base de datos disponibles en BOLD (i.e. `COX1`, `COX1_SPECIES`, `COX1_SPECIES_PUBLIC` y `COX1_L640bp`).

La función es la siguiente:
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

Si queremos tener la identificación a nivel de especie, por ejemplo, de las secuencias guardadas en [secuencias.txt](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/secuencias.txt), podemos usar el siguiente código:

```R
out <- ID_engine(query = read.FASTA('secuencias.txt'), db = "COX1_SPECIES")
```
Vemos las 7 primeras filas y las columnas 1, 5 y 6 de cada elemento de la lista de resultados `out`:
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
