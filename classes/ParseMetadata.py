#   This class takes a metadata json file from the TCGA
#   and returns a comprenhensive dataframe to match cases and
#   alliquots
#    
#   Author: Santiago García Martín
#   Date:   6-6-2019


import json
import pandas as pd

class MetaParser:

    def __init__(self, sourceMetaData):
        print("Reading metadata...")
        self.ReadJSON(sourceMetaData)

    
    def ReadJSON(self, sourceMetaData):
        
        with open(sourceMetaData) as json_file:
            rawData = json.load(json_file)
        
        self.PacientCases = dict()

	#Fetch cases and aliquots ids together
        for pacients in rawData[0]['associated_entities']:
            self.PacientCases[pacients['entity_id']] = pacients['case_id']


    

