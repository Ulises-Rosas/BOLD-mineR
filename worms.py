#!/usr/bin/python3

import re
import urllib

class worms:

    def __init__(self, taxon):

        self.taxonomic_ranges = ["Phylum",
                                 "Subphylum",
                                 "Superclass",
                                 "Class",
                                 "Infraclass",
                                 "Subclass",
                                 "Superorder",
                                 "Order",
                                 "Suborder",
                                 "Superfamily",
                                 "Family",
                                 "Subfamily",
                                 "Genus",
                                 "Species"]

        self.taxon = taxon.replace(" ", "%20")

        aphiaID_url = "http://www.marinespecies.org/rest/AphiaIDByName/" + \
                      self.taxon + \
                      "?marine_only=true"

        self.aphiaID = urllib.request.urlopen(aphiaID_url).read().decode('utf-8')

        self.records_url = "http://www.marinespecies.org/rest/AphiaChildrenByAphiaID/" + \
                      self.aphiaID + \
                      "?marine_only=true&offset=1"
        self.accepted_name = ""

    def get_children_names(self, till = "Species"):

        records_url = 'http://www.marinespecies.org/rest/AphiaChildrenByAphiaID/205965?marine_only=true&offset=1'

        page = urllib.equest.urlopen(records_url).read().decode('utf-8')

        names = [names.replace('\"','').replace('valid_name:', '') for names in re.findall('"valid_name":"[A-Z][a-z]+[ a-z]+"', page)]

        ## in progress
        pass

    def get_accepted_name(self):

        if len(self.aphiaID) == 0:

            species_binary = self.taxon.split("%20")
            #species_binary = worms("Thais luteostoma").taxon.split("%20")

            genus_id_url = "http://www.marinespecies.org/rest/AphiaIDByName/" + species_binary[0] + "?marine_only=true"

            genus_id = urllib.request.urlopen(genus_id_url).read().decode('utf-8')

            complete_url = "http://www.marinespecies.org/aphia.php?p=taxdetails&id=" + genus_id

            page = urllib.request.urlopen(complete_url).read().decode('utf-8')


            lines = re.findall(".*>Thais[\(A-Za-z\) ]{0,} [a-z]+<.*", page)
            epitopes = [re.findall("<i>[A-Z][a-z]+[\(\)A-Za-z ]{0,} [a-z]+</i>", i)[0].\
                            split(" ")[-1].\
                            replace("</i>", "") for i in lines]


            def get_pieces(string, amplitude):

                pieces = [string[i:i + amplitude] for i in range(len(list(string)))]

                trimmed_pieces = [i for i in pieces if len(i) > amplitude - 1]

                return trimmed_pieces

            for index in range(len(list(species_binary[1])) - 1):

                a = get_pieces(species_binary[1], index + 1)

                lengths = []

                for string in epitopes:

                    matches = [re.findall(i, string) for i in set(a)]

                    if len(get_pieces(string, index + 1)) == 0:
                        d1 = 1
                    else:
                        d1 = len(get_pieces(string, index + 1))

                    n1 = sum([len(c) for c in matches])

                    d2 = len(set(a))
                    n2 = len(set(a) - set(["".join(set(b)) for b in matches if len(b) > 0]))

                    lengths.append(n1/d1 + 1 - n2/d2)

                if len([d for d in lengths if d == max(lengths)]) == 1:
                    page_line = lines[lengths.index(max(lengths))]

                    self.accepted_name = re.findall("<i>[A-Z][a-z]+ [a-z]+</i>", page_line)[-1].replace("<i>", "").replace("</i>", "")

                    self.aphiaID = re.findall("aphia.php\?p=taxdetails&id=[0-9]+", page_line)[0].replace("aphia.php?p=taxdetails&id=","")

                    break

            return self.accepted_name

        else:
            complete_url = "http://www.marinespecies.org/aphia.php?p=taxdetails&id=" + self.aphiaID

            page = urllib.request.urlopen(complete_url).read().decode('utf-8')

            if len(re.findall(">unaccepted<", page)) == 1:

                line = re.findall("p=taxdetails&id=(?!" + self.aphiaID + ").*<i>[A-Z][a-z]+ [a-z]+</i>", page)[0]
                self.accepted_name = re.sub(".*</i><i>(.*)</i>", "\\1", line)

                return self.accepted_name

            else:
                self.accepted_name = self.taxon.replace("%20", " ")

                return self.accepted_name

    def taxamatch(self):

        complete_url = "http://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames%5B%5D=" + \
                       self.taxon + \
                       "&marine_only=true"

        page = urllib.request.urlopen(complete_url).read().decode('utf-8')

        self.accepted_name = re.sub('.*"valid_name":"([A-Z][a-z]+ [a-z]+)".*',"\\1", page )

        return self.accepted_name
