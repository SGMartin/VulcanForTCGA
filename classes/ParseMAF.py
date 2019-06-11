#   This class takes a single nucleotide variation MAF file from the TCGA
#   and returns a pd dataframe ready to further analysis
#    
#   Author: Santiago García Martín
#   Date:   03-6-2019

import pandas as  pd

class ParseMAF:
 
    #TODO: relevantColumns class-wide manip. through class methods ?
    #TODO: Is annotation files loading instance scoped? class scoped?
    #TODO: Should we kill RawMaF to save of memory?

    #Trim the dataset, store relevant columns. This is a class-wide array.
    relevantColumns = [0,8,15,39,41,42,115]

    TruncatedStatus = ["De_novo_Start_OutOfFrame","Frame_Shift_Del","Frame_Shift_Ins",
    "In_Frame_Del","In_Frame_Ins","Nonsense_Mutation","Nonstop_Mutation","Splice_Site",
    "Start_Codon_Del","Start_Codon_Ins","Stop_Codon_Del", "Stop_Codon_Ins"]
        
    def __init__(self, sourceMAF, vulcanGenes):

        #Load empty, instance-scope dataframe
        self.MutFilteredData = pd.DataFrame()

        print("Starting MAF analysis...")

        self.RawMaF = pd.read_csv(sourceMAF, sep='\t', header=5, low_memory=False) #TODO: handle the '#' headers instead of skipping 4 lines
       
       #Filter the Dataframe
        self.FilterForVulcan(vulcanGenes) #TODO: this argument passing could be improved
        self.AnnotateCancerRole()
        self.AnnotateOncoDrive()
        self.AnnotateTruncatedMissense()
        self.CalculateVariantFrequency()

        print("MAF analysis COMPLETED!")

    def FilterForVulcan(self, vulcanGeneList):

        #TODO: This method could be split in two and reused between ParseMAF and CNV-cleaning
        self.RawMaF = self.RawMaF.iloc[:, self.relevantColumns]

        #Trim the datased based on Vulcan genes available
        self.VulcanVector = self.RawMaF.iloc[:,0].isin(vulcanGeneList)
        self.MutFilteredData = self.RawMaF[self.VulcanVector]

        #Rename columns
        self.MutFilteredData.columns = ['Gene','Consequence','TumorBarCode','TotalTumorDepth','VariantDepth','NormalDepth', 'CaseID']

        #Drop entries with low depth (<30 depth)
        minimumTumoralDepth =  self.MutFilteredData['TotalTumorDepth'] >= 30
        minimumNormalDepth  =  self.MutFilteredData['NormalDepth'] >= 30
        self.MutFilteredData = self.MutFilteredData[(minimumNormalDepth & minimumTumoralDepth)]

        #Drop silent mutations
        SilentMutations         = self.MutFilteredData['Consequence'] == "Silent"
        self.MutFilteredData = self.MutFilteredData[~SilentMutations]

        #Clean unnecessary rows
        self.MutFilteredData = self.MutFilteredData.drop('NormalDepth', axis=1)
        #Reset the index
        self.MutFilteredData = self.MutFilteredData.reset_index(drop=True)

    
    def AnnotateCancerRole(self):

        #Load CGC .TSV file
        self.CancerGeneCensus = pd.read_csv("/home/sagarcia/CancerGeneCensus.tsv", sep='\t')
        
        #Annotate gene role from CGC
        self.CancerGeneCensus   = self.CancerGeneCensus.loc[:,['Gene Symbol', 'Role in Cancer']]
        self.MutFilteredData = self.MutFilteredData.merge(self.CancerGeneCensus, how='left', left_on='Gene', right_on='Gene Symbol')
        self.MutFilteredData = self.MutFilteredData.drop('Gene Symbol', axis=1)
    
            
    #Further annotation based on https://bbglab.irbbarcelona.org/oncodriverole/
    #Data downloaded 6-6-19
    def AnnotateOncoDrive(self):

        #Load oncodrive table        
        self.OncoDrive = pd.read_csv("/home/sagarcia/oncodrive.txt", sep='\t')

        #Annotate inferred LoF/GoF
        self.OncoDrive = self.OncoDrive.loc[:,['SYM', 'oncodriveROLE', 'Value']]
        self.OncoDrive.columns = ['GenSYM', 'OncoDrive','DriveValue']
        self.MutFilteredData = self.MutFilteredData.merge(self.OncoDrive, how='left', left_on='Gene', right_on='GenSYM')
        self.MutFilteredData = self.MutFilteredData.drop('GenSYM', axis=1)

        #I replace Loss of function and Activating with activating and loss for later convenience
        self.MutFilteredData['OncoDrive'] = self.MutFilteredData['OncoDrive'].str.replace('Loss of function', 'loss')
        self.MutFilteredData['OncoDrive'] = self.MutFilteredData['OncoDrive'].str.replace('Activating', 'activating')    
    
    def AnnotateTruncatedMissense(self):
        
        self.MutFilteredData.loc[:,'Truncated'] = self.MutFilteredData['Consequence'].isin(self.TruncatedStatus)
        self.MutFilteredData.loc[:,'Missense']  = self.MutFilteredData['Consequence'] == 'Missense_Mutation'
    

    def CalculateVariantFrequency(self):

        #Add column with VAF calculated
        self.MutFilteredData.loc[:,'VAF']      = self.MutFilteredData['VariantDepth'] / self.MutFilteredData['TotalTumorDepth']
        self.MutFilteredData.loc[:,'Missense'] = self.MutFilteredData['Consequence'] == 'Missense_Mutation'

        #Drop depths columns. Not needed anymore
        self.MutFilteredData = self.MutFilteredData.drop(['TotalTumorDepth', 'VariantDepth'], axis=1)