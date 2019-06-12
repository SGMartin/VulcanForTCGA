#   This class summarises the data based on gene alterations 
#   Gene alterations are clasified on six categories: truncated,missense_activating,missense_unknown,
#   missense_loss,copy_gain,copy_loss
#   Role in cancer is taken from the Cancer Gene Census
#   Author: Santiago García Martín
#   Date:   10-6-2019

import numpy  as np
import pandas as pd
import requests

class VulcanGeneReport:

    def __init__(self, mafFilteredData, cnvFilteredData):
        
        print("Vulcan query started")

        self.VulcanResponses = {}
        self.QueryVulcan(mafFilteredData, cnvFilteredData)

    def QueryVulcan(self, mutationData, cnvData):

        #query only unique gene entries, we'll sort after that
        self.MafGenes = list(mutationData['Gene'])
        self.CNVGenes = list(cnvData['Gene'])
        
        self.AllGenes = self.MafGenes + self.CNVGenes

        #final gene list
        self.GenesToQuery = np.unique(self.AllGenes)
        print("Querying Vulcan for " + str(len(self.GenesToQuery)) + " genes")

        for gene in self.GenesToQuery:
            print(gene, "...")
            response = requests.get(("http://vulcanspot.org/api/genes/" + str(gene) + "/treatments"))
            self.VulcanResponses[str(gene)] = response.json()

        print("Vulcan query COMPLETED")

            


