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

        page = urllib.request.urlopen(records_url).read().decode('utf-8')

        names = [names.replace('\"','').replace('valid_name:', '') for names in re.findall('"valid_name":"[A-Z][a-z]+[ a-z]+"', page)]

        ## in progress
        pass

    def get_accepted_name(self):

        complete_url = "http://www.marinespecies.org/aphia.php?p=taxdetails&id=" + self.aphiaID

        page = urllib.request.urlopen(complete_url).read().decode('utf-8')

        if len(re.findall(">unaccepted<", page)) == 1:

            line = re.findall("p=taxdetails&id=(?!"+ self.aphiaID +").*<i>[A-Z][a-z]+ [a-z]+</i>", page)[0]
            self.accepted_name = re.sub(".*</i><i>(.*)</i>", "\\1", line)

            return self.accepted_name

        else:
            self.accepted_name = self.taxon.replace("%20", " ")

            return self.accepted_name

worms("Manta birostis").get_accepted_name()
