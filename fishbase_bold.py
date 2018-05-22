import re
import urllib
import pandas as pd
import argparse
##please delete me later on
parser = argparse.ArgumentParser(description='Short module for getting insights about DNA barcode given a list of species')

parser.add_argument('fastq', metavar='SpeciesList',
                    # nargs = "+",
                    type=argparse.FileType('r'),
                    help='Target species in text format.')

args = parser.parse_args()

all_species = [i.replace('\n', '') for i in args.fastq]

class fishbase:
    '''Validate names against to the fishbase database
    '''
    def __init__(self, species):
        '''It is just a single entry.
        Here is where the species name should enter
        as initial value.
        '''
        self.species = species.replace('\n', '')
        self.synonyms = []


        self.summary_url = 'https://www.fishbase.de/summary/'
        self.api_url = 'https://fishbase.ropensci.org/'
        self.synonym_url = 'http://www.fishbase.se/Nomenclature/SynonymsList.php?'


    def validate_name_html(self):
        '''The initial value is confronted to the Fishbase repository
        to validate whether it really exists
        complete_url = 'https://www.fishbase.de/summary/Sarda-chilnss.html'
        '''
        complete_url = self.summary_url  + self.species.replace(' ', '-')  + '.html'
        page = urllib.request.urlopen(complete_url).read().decode('utf-8')

        f = open(self.species + '.html', 'w')
        f.write(page)
        f.close()

    def validate_name(self):
        '''The initial value is confronted to the Fishbase repository
        to validate whether it really exists
        complete_url = 'https://www.fishbase.de/summary/Sarda-chilnss.html'
        '''
        complete_url = self.summary_url + self.species.replace(' ', '-') + '.html'

        print("Accessing to: {}".format(complete_url))

        page = urllib.request.urlopen(complete_url).read().decode('utf-8')

        empty_html_body = "Species name is not in the public version of FishBase."

        if page.find(empty_html_body) != -1:
            return empty_html_body

        else:

            IDs = [i.replace('SynonymsList.php?ID=', '') for i in re.findall("SynonymsList\.php\?ID=[0-9]+", page)]

            if len(set(IDs)) == 1:
                return r"Species ID : {}".format(IDs[0])

            else:
                return r"Species IDs: {}".format(
                    str(set(IDs)). \
                        replace("'", ""). \
                        replace("{", ""). \
                        replace("}", "")
                )


valid_names = []

all_ids = []

for spps in all_species:

    fish_tmp = fishbase(spps)
    ID = fish_tmp.validate_name()

    all_ids.append(ID)

    if len(re.findall("Species ID", ID)) == 1:
        valid_names.append(spps)


table0 = pd.concat([
    pd.Series(all_species, name='Species'),
    pd.Series(all_ids, name='Species IDs')],
    axis=1)

table0.to_csv("Names_and_IDs.txt", header=True, index=False, sep= "\t") ##possible bug from the sep argument

print(f"\n {table0} \n\nNames with their corresponding were stored at Names_and_IDs.txt")



class bold:

    def __init__(self, all_species):
        '''It gets just a list of entries.
        Here is where the species names should enter
        as initial value. These values are taken from arg.parser objects
        '''
        self.names = all_species

        self.speciesName = []
        self.Nseqs = []

    def get_Meta_and_Seqs(self):


        file = open('sequences.fasta', 'w')

        def export_lines(byte_list):
            for i in byte_list:
                file.write(i.decode('utf-8'))

        for query in self.names:

            base_url = 'http://www.boldsystems.org/index.php/API_Public/sequence?taxon='

            if len(re.findall('\s', query)) == 1:
                self.speciesName.append(query)

                complete_url = base_url + query.replace(' ', '%20')
                seqs_object = urllib.request.urlopen(complete_url)

                print( "Accessing to: {}".format(complete_url) )

                seqs = seqs_object.readlines()

                export_lines(seqs)

                self.Nseqs.append(
                    sum( [ len(re.findall('>', i.decode('utf-8'))) for i in seqs ] )
                    )

        file.close()

        table = pd.concat([
            pd.Series(self.speciesName, name='Species'),
            pd.Series(self.Nseqs, name='Number of sequences')],
            axis=1)

        table.to_csv("table.csv", header=True, index=False)

        print(f"\n {table} \n\nMetadata and sequences were stored in table.csv and sequences.fasta respectively")


bold_tmp = bold(all_species)

bold_tmp.get_Meta_and_Seqs()
