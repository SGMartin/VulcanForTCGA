#   This class uses filtered MAF and CNV information to
#   rebuild the pacient genomic landscape. Then the information
#   retrieved by Vulcan is filtered based on this pacient genetic events

#   NOTE that there is no need to call Vulcan once per pacient, as the 
#   cohort vulcan query contains all eligible genes already.

#   Author: Santiago García Martín
#   Date:   03-6-2019

import Utils
import pandas as pd
import numpy  as np

class Pacients:

    def __init__(self, uniqueID, mutDataFrame, cnvDataFrame, vulcanResults):

        self.UID              = uniqueID
        self.MutData          = mutDataFrame
        self.CNVData          = cnvDataFrame
        self.GeneticLandscape = pd.DataFrame

        self.GetCNVForPacient()
        self.GetMutationTable()
        self.GeneticLandscape.to_csv(('/home/sagarcia/Desktop/Report/' + 'GLand.csv'), index=None, header=True)
    
    def GetCNVForPacient(self):
        
        #don't worry about multiple bar codes, the participant should be the same
        #https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

        self.BarCode = self.MutData['TumorBarCode']
        self.BarCode = list(self.BarCode)
        self.BarCode = self.BarCode[0]

        self.PacientBarCode = self.BarCode.split('-')[2]

        #get all CNV columns with matching pacient bar code
        self.pacientSamples = [col for col in self.CNVData if self.PacientBarCode in col ]
        self.columnsOfInterest = ['Gene', 'Role in Cancer', 'OncoDrive', 'DriveValue']
        self.columnsOfInterest = self.columnsOfInterest + self.pacientSamples

        self.CNVData = self.CNVData.loc[:,self.columnsOfInterest]

        #Now, add sample type column and melt it. Just in case there are more samples 
        #than one per case id

        self.CNVData = pd.melt(self.CNVData, id_vars=['Gene', 'Role in Cancer', 'OncoDrive', 'DriveValue'],
        var_name='BarCode', value_name='Alteration')

        
        self.AltDict = {0 : 'cnv_neutral', 1: 'cnv_gain', -1:'cnv_loss'}
        self.CNVData = self.CNVData.replace({'Alteration': self.AltDict})

        #Rename Sample from TCGA-XXXXXX- to TP/TM etc...
        self.CNVData['Sample'] = self.CNVData['BarCode'].apply(Utils.TranslateSampleType)
        self.CNVData           = self.CNVData.drop(['BarCode', 'OncoDrive', 'DriveValue'], axis=1)

        


    
    def GetMutationTable(self):

        self.MutData = self.MutData.filter(['Gene', 'Consequence', 'Role in Cancer',
       'OncoDrive', 'Truncated', 'Missense', 'Sample'], axis=1)

       #Gets rid of the _mutation suffix and lowercases everything
        self.MutData['Consequence'] = self.MutData['Consequence'].str.split('_').str[0].str.strip().str.lower()
        
        #Replace  NaN with unknown
        self.MutData['OncoDrive'] = self.MutData['OncoDrive'].replace(np.nan, 'unknown', regex=True)

        #Those wich the tag Truncated == True, will be called truncated now
        #Remember that if truncated == false, missense == true atm
        #missense + _activating, missense + _loss, missense +_NaN


        self.MutData['Alteration'] = np.where(self.MutData['Truncated'] == True, 'truncated', 
        self.MutData['Consequence'] + '_' + self.MutData['OncoDrive'])

        self.MutData = self.MutData.drop(['Missense', 'Truncated', 'Consequence', 'OncoDrive'], axis=1)

        self.GeneticLandscape = self.MutData.append(self.CNVData)
       