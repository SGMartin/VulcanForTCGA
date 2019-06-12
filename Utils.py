#   
#   These functions annotates the dataframes using COSMIC for 'Role in cancer'
#   and oncodrive for gene behaviour prediction
#
#   Author: Santiago García Martín
#   Date:   11-6-2019



#Taken from COSMIC cancer gene census
def AnnotateCancerRole(rawData):

    import pandas as pd

    CancerGeneCensus = pd.read_csv("./Data/CancerGeneCensus.tsv", sep='\t')
       
    #Annotate gene role from CGC
    CancerGeneCensus = CancerGeneCensus.loc[:,['Gene Symbol', 'Role in Cancer']]

    AnnotatedData  = rawData.merge(CancerGeneCensus,
    how='left', 
    left_on='Gene', 
    right_on='Gene Symbol')

    AnnotatedData  = AnnotatedData.drop('Gene Symbol', axis=1)

    return AnnotatedData


#Annotate gene role from CGC      
#Further annotation based on https://bbglab.irbbarcelona.org/oncodriverole/
#Data downloaded 6-6-19

def AnnotateOncoDrive(rawData):

    import pandas as pd

    OncoDrive = pd.read_csv("./Data/oncodrive.tsv", sep='\t')
    
    #Annotate inferred role
    OncoDrive = OncoDrive.loc[:,['SYM', 'oncodriveROLE', 'Value']]
    OncoDrive.columns = ['GenSYM', 'OncoDrive','DriveValue']

    AnnotatedData = rawData.merge(OncoDrive, how='left', left_on='Gene', right_on='GenSYM')
    AnnotatedData = AnnotatedData.drop('GenSYM', axis=1)

    #I replace Loss of function and Activating with activating and loss for later convenience
    AnnotatedData['OncoDrive'] = AnnotatedData['OncoDrive'].str.replace('Loss of function', 'loss')
    AnnotatedData['OncoDrive'] = AnnotatedData['OncoDrive'].str.replace('Activating', 'activating')

    return AnnotatedData

def TranslateSampleType(barcode):

    #https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
    CodesDictionary = {
        "01":"TP",
        "02":"TR",
        "03":"TB",
        "04":"TRBM",
        "05":"TAP",
        "06":"TM",
        "07":"TAM",
        "08":"THOC",
        "09":"TBM",
        "10":"NB",
        "11":"NT",
        "12":"NBC",
        "13":"NEBV",
        "14":"NBM",
        "15":"15SH",
        "16":"16SH",
        "20":"CELLC",
        "40":"TRB",
        "50":"CELL",
        "60":"XP",
        "61":"XCL",
        "99":"99SH"
    }
    sampleCode = str(barcode).split('-')
    doubleDigit = sampleCode[3][0:2]

    return CodesDictionary[doubleDigit]
    



   