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
        print("Metadata loaded")

    
    def ReadJSON(self, sourceMetaData):
        
        with open(sourceMetaData) as json_file:
            rawData = json.load(json_file)
        
        self.PacientBarCodes = dict()

	#Fetch cases and aliquots ids together Entity submitter id is TCGA bar code
        for pacients in rawData[0]['associated_entities']:
            self.PacientBarCodes[pacients['entity_id']] = pacients['entity_submitter_id']


    

