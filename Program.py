#
# C# >>>>>> Python :P
import sys
import pandas as pd
import classes.VulcanGeneWrapper as VulcanGeneWrapper
import classes.ParseMAF as MAFParser
import classes.ParseCNV as CNVParser
import classes.ParseMetadata as MetaParser
import classes.VulcanGeneReport as Reporter

ReportDirectory = ""
SourceFiles     = ""


def main():
    #TODO: Remember to uncomment this
  #  SourceFiles     = input("Input project directory:\n")
  #  ReportDirectory = input("Input the directory to save the report to:\n")
    
    SourceFiles     = "/home/sagarcia/rawtest/"
    ReportDirectory = "/home/sagarcia/Desktop/Report/"
    VulcanContext   =  ""
    MutationData   = pd.DataFrame
    CopyNumberData = pd.DataFrame

    ParseRawData()
    GenerateReports(ReportDirectory)

    
def ParseRawData():

    global MutationData
    global CopyNumberData

    #Call VulcanSpot to ask for input gene list
    VulcanGenes = VulcanGeneWrapper.VulcanInputGenes("/home/sagarcia/Desktop/Report/", "PANCREAS")

    #Call the wrappers
    MetaData        = MetaParser.MetaParser("/home/sagarcia/rawtest/metadata.json")
    MutationData    = MAFParser.ParseMAF("/home/sagarcia/rawtest/test.maf",VulcanGenes.VulcanGeneList)
    CopyNumberData  = CNVParser.ParseCNV("/home/sagarcia/rawtest/CNV.txt", VulcanGenes.VulcanGeneList)
    
    #TODO: This should not be called here
    CopyNumberData.AnnotateCases(MetaData.PacientCases)

    #TEST
    test = Reporter.VulcanGeneReport(MutationData, CopyNumberData)

    #TODO: move this to classes?
def GenerateReports(reportDirectory):

    print("Generating reports...")

    #Write filtered wrapper data to a comma separated file
    MutationData.MutFilteredData.to_csv((reportDirectory + 'Mutations.csv'), index=None, header=True)
    CopyNumberData.CNVFilteredData.to_csv((reportDirectory + 'CNV.csv'), index=None, header=True)

    print("DONE!!!!")


#If the python source file is imported as module, 
#python interpreter sets the __name__ value to module name, 
#so the if condition will return false and main method will not be executed.    
if __name__ == '__main__':
    main()
