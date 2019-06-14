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
        self.AlternativeDrugTable = pd.DataFrame(columns=['MutatedGene','Druggable', 'AlternateDrug', 'Target', 'Rank'])
        self.DirectDrugsTable = pd.DataFrame(columns=['MutatedGene', 'Drug', 'Dscore'])

        self.QueryVulcan()
        self.AlternativeDrugTable.to_csv(('/home/sagarcia/Desktop/Report/DrugTable.csv'), index=None, header=True)
        self.DirectDrugsTable.to_csv(('/home/sagarcia/Desktop/Report/DirectDrugTable.csv'), index=None, header=True)


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
            self.GetAllDrugsForGene(gene, self.TypeOfCancer, response)


    def GetAllDrugsForGene(self, inputGene, typeOfCancer, vulcanData):

        self.IsDruggable = False
        self.GeneStatus  = self.Alterations['Impact'][self.Alterations['Gene'] == inputGene].unique()

        self.TumorDrugData = vulcanData['data'][inputGene][typeOfCancer]

        if 'drugs' in self.TumorDrugData:
            self.IsDruggable = True
            
            #gene A drugs
            for drug,score in self.TumorDrugData['drugs'].items():
                self.DirectDrugsTable = self.DirectDrugsTable.append({'MutatedGene': inputGene,
                'Drug': drug, 'Dscore' : score}, ignore_index=True)


        
        #Let's search for indirect targets nonetheless

        for status in self.GeneStatus:  #check if this gene is LoF/GoF and if there are Genetic Dep. for it
            if status in self.TumorDrugData['alterations']:
                #loop through genetic dependencies
                    for GeneticDependency in self.TumorDrugData['alterations'][status]:
                        #check if there are drugs for this genetic dependency
                            if 'drugs' in self.TumorDrugData['alterations'][status][GeneticDependency]:
                                #found drugs, grab them
                                    for drug, drugData in self.TumorDrugData['alterations'][status][GeneticDependency]['drugs'].items():
                                         self.DrugRank = Utils.GetDrugRank(drugData)
                                         self.AlternativeDrugTable = self.AlternativeDrugTable.append({
                                             'MutatedGene' : inputGene,
                                             'Druggable' : self.IsDruggable,
                                             'AlternateDrug': drug,
                                             'Target': GeneticDependency,
                                             'Rank' : self.DrugRank},
                                         ignore_index=True)


        






        


       

            


