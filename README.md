# BOLD-mineR

DNA barcode insights no only is used by researchers, but also by decision-makers for their corresponding projects such as those for food fraud or illegal species commercialization. Conversely to its big-scaled demand of both online services and information, there are few ways to automatize either species identification or assessment of barcode quality per species directly from the web interface. 

[BOLD system](http://www.boldsystems.org/) is the main database of DNA barcode worldwide. This database has been stepply  growing through time since its release ([Ratnasingham and Hebert 2007](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1471-8286.2007.01678.x)) and its accessibility is pivotal for projects regards DNA barcodes. Nowdays APIs, to some extant, are offering access to well-know databases such as [FishBase](https://fishbase.ropensci.org/), [Worms](http://www.marinespecies.org/rest/) or [BOLD](http://www.boldsystems.org/index.php/api_home). Depite of BOLD's API mostly involve public data only, this leverages its data retrieving for wider purposes. The API's applicability, however,  seems to be wholly held up by its own needs of having either standalone softwares or functions which could wrap up blocks of information. The main objective of these functions (i.e. BOLD-mineR's functions) is justly circumscribe the BOLD's API performance with R-based scripts to get insights about DNA barcodes by using public information.


## SpecimenData 

This function **SpecimenData** let us mine associated metadata from any specimen according to following arguments:

* `taxon` (e.g. Aves|Elasmobranchii).
* `ids` (e.g. ANGBF12704-15).
* `bin` (e.g. BOLD:AAA4689).
* `container` (e.g. FIPP).
* `institution` (e.g. Smithsonian Institution).
* `researchers` (including identifiers and collectors).
* `geo` (e.g. Peru).

You can find this function here: [SpecimenData](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/SpecimenData.R)

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

**ID_engine** finds best matches between a query sequence and a database of BOLD by using BLASTn-based algorithms. Arguments of this function are `query` and `db`. The first one are query sequences and the second one are one of avilable databases in BOLD:

* `COX1`
* `COX1_SPECIES`
* `COX1_SPECIES_PUBLIC`
* `COX1_L640bp`

This script also take account for those sequences which are not, by mistake, at the right sense and sends them to perform a BLAST search through its [API](https://ncbi.github.io/blast-cloud/dev/api.html)

You can find this function here: [ID_engine](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/ID_engine.R)

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

This function adds an **audition** step ([Oliveira _et al._ 2016](https://onlinelibrary.wiley.com/doi/full/10.1111/jfb.13169)) to each selected specimen by `ID_engine()` (see above), given a certain threshold. This function, in turn, uses another function called `AuditionBarcodes()`. This last one has two version. The first one is coupled with `addAudition()` and the second one is `addAudition()`-independent and also normalizes species names by taking accepted names from [Worms database](http://www.marinespecies.org/).


Both versions of `AuditionBarcodes()` can be found here:

* [AuditionBarcodes](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/AuditionBarcodes.R)
* [AuditionBarcodes.v.2](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/AuditionBarcode.v.2.R) which in turn is coupled with [worms.py](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/worms.py)

`addAudition()` function can be found here: 

* [addAudition](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/ID_engine.R)

In order to test its efficiency, species-level identification of samples stored in [secuencias.txt](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/secuencias.txt) is conducted by using a `threshold = 0.99` (i.e. 99% of similarity):

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

#### AuditionBarcodes

Despite `AuditionBarcodes()` function is coupled with `addAudition()` function, its independent version can work with just a list of names. Furthermore, there is an argument which enables to chose if sequences from GenBank are considered. It is pending, however, assess whether these sequences used to assess barcode's quality come from either a published article or direct submission:

```R
species <- c("Ceratostoma foliatum" , "Ocinebrina aciculata",
              "Morula margariticola", "Concholepas concholepas",
              "Nucella emarginata" ,"Thais luteostoma", "Thais bufo")
              
quality_assessment <- AuditionBarcodes(species, include_ncbi = F)

data.frame(species, quality_assessment)
#                  species Grades                                                         Observations                                                                                                               BIN_structure
#1    Ceratostoma foliatum      A                                 Matched BIN with external congruence                                                                                  'BOLD:ABV5037':{'Ceratostoma foliatum':32}
#2    Ocinebrina aciculata      B                            Matched BIN with internal congruence only                                                                                  'BOLD:ACG0291':{'Ocinebrina aciculata':18}
#3    Morula margariticola      C                                                         Splitted BIN                                    'BOLD:AAD8263':{'Drupella margariticola':5}, 'BOLD:AAE7335':{'Drupella margariticola':1}
#4 Concholepas concholepas      D Insufficient data. Institution storing: 1. Total specimen records: 1                                                                                                                            
#5      Nucella emarginata     E*                                                           Merged BIN                                                                   'BOLD:ABA3756':{'Nucella emarginata':16,'Nucella lima':7}
#6        Thais luteostoma    E**                                                         Mixtured BIN 'BOLD:ACB7390':{'Reishia bronni':3,'Reishia luteostoma':16}, 'BOLD:AAW6905':{'Reishia clavigera':11,'Reishia luteostoma':1}
#7              Thais bufo      F                           Barcodes mined from GenBank or unvouchered                                                                                                                            
```

Please notice that grades are obtained with accepted names of species according to [Worms database](http://www.marinespecies.org/) by using the [worms.py](https://github.com/Ulises-Rosas/BOLD-mineR/blob/master/worms.py) script. Hence, since currently accepted names within `species` vector has not been figured out, unevenness between the column `BIN_structure` and `species` could pop up. For this reason, on above example, _Morula margariticola_ has another name into BIN structure column (i.e. _Drupella margariticola_) which actually is its currently accepted name.