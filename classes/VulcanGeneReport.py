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
        
        print("Generating genetic report...")

        self.MutSummary = pd.DataFrame
        self.CNVSummary = pd.DataFrame

        self.VulcanResponses = {}

        self.PrepareMAFData(mafFilteredData)
        self.PrepareCNVData(cnvFilteredData)
        self.QueryVulcan()



    #Further cleans the MAF data and summarises genetic alterations in one column
    def PrepareMAFData(self, mafFilteredData):

        self.MutSummary = mafFilteredData.MutFilteredData.filter(['Gene', 'Consequence', 'Role in Cancer', 'OncoDrive',
         'Truncated', 'Missense'], axis=1)

        #Gets rid of the _mutation suffix and lowercases everything
        self.MutSummary['Consequence'] = self.MutSummary['Consequence'].str.split('_').str[0].str.strip().str.lower()
        
        #Replace  NaN with unknown
        self.MutSummary['OncoDrive'] = self.MutSummary['OncoDrive'].replace(np.nan, 'unknown', regex=True)

        #Those wich the tag Truncated == True, will be called truncated now
        #Remember that if truncated == false, missense == true atm
        #missense + _activating, missense + _loss, missense +_NaN


        self.MutSummary['Alteration'] = np.where(self.MutSummary['Truncated'] == True, 'truncated', 
        self.MutSummary['Consequence'] + '_' + self.MutSummary['OncoDrive'])

        self.MutSummary.to_csv(('/home/sagarcia/Desktop/Report/' + 'MutSummary.csv'), index=None, header=True)

    #Counts CNV gain and loss across all samples and returns a dataframe 
    def PrepareCNVData(self, cnvFilteredData):
        
        self.CNVSummary = cnvFilteredData.CNVFilteredData.apply(pd.Series.value_counts, axis=1)[[-1,0,1]].fillna(0)

        self.CNVSummary.columns = ['Loss', 'Neutral', 'Gain']
        self.CNVSummary['Gene'] = cnvFilteredData.CNVFilteredData['Gene Symbol']
        self.CNVSummary.reset_index()
        self.CNVSummary.to_csv(('/home/sagarcia/Desktop/Report/' + 'CNVSummary.csv'), index=None, header=True)


    #let's rock
    def QueryVulcan(self):

        #query only unique gene entries, we'll sort after that
        self.MafGenes = list(self.MutSummary['Gene'])
        self.CNVGenes = list(self.CNVSummary['Gene'])
        self.GenesToQuery = self.MafGenes + self.CNVGenes

        #final gene list
        self.GenesToQuery = np.unique(self.GenesToQuery)
        print("Querying Vulcan for " + str(len(self.GenesToQuery)) + " genes")

        for gene in self.GenesToQuery:
            print(gene, "...")
            response = requests.get(("http://vulcanspot.org/api/genes/" + str(gene) + "/treatments"))
            self.VulcanResponses[str(gene)] = response.json()

            


