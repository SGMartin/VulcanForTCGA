#   This class summarises the data based on gene alterations 
#   Gene alterations are clasified on six categories: truncated,missense_activating,missense_unknown,
#   missense_loss,copy_gain,copy_loss
#   Role in cancer is taken from the Cancer Gene Census
#   Author: Santiago García Martín
#   Date:   10-6-2019

import Utils
import numpy  as np
import pandas as pd
import requests

class VulcanGeneReport:

    def __init__(self, mutationAndCNV, context):

        self.VulcanResponses = {}
        self.TypeOfCancer    = context
        self.Alterations     = mutationAndCNV
        self.DrugTable       = pd.DataFrame(columns=['Drug', 'GeneA', 'GeneB', 'Rank'])

        self.QueryVulcan()
        self.DrugTable.to_csv(('/home/sagarcia/Desktop/Report/DrugTable.csv'), index=None, header=True)


    def QueryVulcan(self):
        
        print("Vulcan query started")

        #Do not query for unknown impact genes 
        unknownMut =self.Alterations['Impact'] == 'Unknown'
        self.GenesToQuery = self.Alterations[~unknownMut]

        #get gene unique entries 
        self.GenesToQuery = list(self.GenesToQuery['Gene'])
        self.GenesToQuery = np.unique(self.GenesToQuery)

       
        print("Querying Vulcan for " + str(len(self.GenesToQuery)) + " genes")

        for gene in self.GenesToQuery:
            print(" Querying ", gene, "...")
            response = requests.get(("http://vulcanspot.org/api/genes/" + str(gene) + "/treatments"))
            self.VulcanResponses[str(gene)] = response.json()
        
        for gene,response in self.VulcanResponses.items():
            self.GetAllDrugsForGeneA(gene,response)

        

    #let's organize the data
    def GetAllDrugsForGeneA(self, geneA, vulcanData):

        print("Getting drugs for ", geneA)

        #even If most of the time there should be only ONE value (GoF/LoF), let's not assume that
        #What if context == 'Pancancer'? let's handle both cases
        genImpact = self.Alterations['Impact'][self.Alterations['Gene'] == geneA].unique()

        for pImpact in genImpact:
            if pImpact in vulcanData['data'][geneA][self.TypeOfCancer]['alterations']:                  
                #let's rename the response for convenience and readability
                self.GeneticDependencies = vulcanData['data'][geneA][self.TypeOfCancer]['alterations'][pImpact]

                for geneB in self.GeneticDependencies:
                    if 'drugs' in self.GeneticDependencies[geneB]:
                        for drug,scores in self.GeneticDependencies[geneB]['drugs'].items():
                            drugRank = Utils.GetDrugRank(scores)
                            self.DrugTable = self.DrugTable.append({'Drug':drug,'GeneA':geneA,'GeneB':geneB,'Rank':drugRank}, ignore_index=True)
                            
        


       

            


