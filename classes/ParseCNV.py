#   This class takes a GISTIC-2 scores CNV file from the TCGA
#   and returns a pd dataframe ready to further analysis
#    
#   Author: Santiago García Martín
#   Date:   6-6-2019

import json
import mygene
import pandas as pd

class ParseCNV:
    
    def __init__(self, sourceFile, vulcanGeneList):

        print("Starting CNV analysis...")

        self.CNVFilteredData = pd.DataFrame
        self.RawCNV = pd.read_csv(sourceFile, sep='\t')

        #drop gene id and the cytoband, we don't need it anymore
        self.CNVFilteredData = self.RawCNV.drop(['Gene ID','Cytoband'],axis=1)

        self.TranslateEnsemblToHugo()
        self.FilterVulcanGenes(vulcanGeneList)

    #Translates ID from Ensembl to Hugo.
    #More info: https://docs.mygene.info/en/latest/doc/data.html

    def TranslateEnsemblToHugo(self):
        
        self.mg = mygene.MyGeneInfo()
        
        #DO NOT query version of EnsemBL IDs, trim them first
        self.CNVFilteredData['Gene Symbol'] = self.CNVFilteredData['Gene Symbol'].str.split('.').str[0].str.strip()

        self.geneQuery = self.CNVFilteredData['Gene Symbol']

        #pandas dataframe with the annotated genes
        print("Annotating CNV...")
        self.AnnotationDT = self.mg.getgenes(self.geneQuery, fields='symbol', as_dataframe=True)

        #map to ensembl... indexes are the same so this ought to work
        self.CNVFilteredData['Gene Symbol'] = self.AnnotationDT['symbol'].values

        #dispose of the dataframe
        del self.AnnotationDT
    
    #TODO: THis is shitty code and you know it. Do smth about it. Do not allow it to depend on above method execution
    #Keeping genes that can be input to vulcan only.
    def FilterVulcanGenes(self, vulcanGenes):

        self.IsInVulcan      = self.CNVFilteredData.loc[:,'Gene Symbol'].isin(vulcanGenes)
        self.CNVFilteredData = self.CNVFilteredData[self.IsInVulcan]

        #redo the indexes
        self.CNVFilteredData.reset_index(drop=True)
    

    #alliquot ID is not needed, pacient ID is
    #TODO: handle odd rename failing...
    def AnnotateCases(self, testDict):
        self.CNVFilteredData.rename(testDict, axis = 'columns', inplace=True)



        

        





