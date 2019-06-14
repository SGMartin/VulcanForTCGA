#   
#   VulcanForTCGA 
#     
#   Author: Santiago García Martín
#   Date:   03-6-2019

import sys
import pandas as pd

import classes.Pacients as Pacients
import classes.ParseMAF as MAFParser
import classes.ParseCNV as CNVParser
import classes.ParseMetadata as MetaParser
import classes.VulcanGeneReport as Reporter
import classes.VulcanGeneWrapper as VulcanGeneWrapper

ReportDirectory = ""
SourceFiles     = ""


def main():
  #TODO: Remember to uncomment this
  #  SourceFiles     = input("Input project directory:\n")
  #  ReportDirectory = input("Input the directory to save the report to:\n")
    
  #Input variables
  SourceFiles     = "/home/sagarcia/rawtest/"
  ReportDirectory = "/home/sagarcia/Desktop/Report/"
  VulcanContext   =  ""
    
#Call VulcanSpot to ask for input gene list
  VulcanGenes = VulcanGeneWrapper.VulcanInputGenes("/home/sagarcia/Desktop/Report/", "PANCREAS")

  MutationData   = ParseMAFData(1, VulcanGenes.GeneList).MutFilteredData #filtered data from maf file
  CopyNumberData = ParseCNVData(1, VulcanGenes.GeneList).CNVFilteredData #filtered data from cnv file  
  MutationAndCNV = AnalyzePacientGeneticData(MutationData, CopyNumberData) #final table summarising CNVs and mutations

  VulcanQuery    = QueryVulcanForGenes(MutationAndCNV, "PANCREAS") #Treatments from VulcanSpot as json response

  GenerateReports(ReportDirectory, MutationData, CopyNumberData, MutationAndCNV)    

def ParseMAFData(rawMafData, vulcanGenes):

  FilteredMutations = MAFParser.ParseMAF("/home/sagarcia/rawtest/test.maf", vulcanGenes)

  return FilteredMutations

def ParseCNVData(rawCNVData, vulcanGenes):

   #Call the wrapper for metadata
    MetaData = MetaParser.MetaParser("/home/sagarcia/rawtest/metadata.json")
    FilteredCNVData = CNVParser.ParseCNV("/home/sagarcia/rawtest/CNV.txt", vulcanGenes)

    FilteredCNVData.AnnotateCases(MetaData.PacientBarCodes)

    return FilteredCNVData

def QueryVulcanForGenes(mafFilteredData, cnvFilteredData):
  
  Results = Reporter.VulcanGeneReport(mafFilteredData, cnvFilteredData)
  return Results


#TODO: ask for analysis method: per tumor/per pacient 
def AnalyzePacientGeneticData(filteredMAF, filteredCNV):

  PacientsData = []
  PacientGroup = filteredMAF.groupby('CaseID')
    
  #split mutation data based on ID cases. Returns a collection of dataframes 
  #with the column as key  
  for caseid, Mutations in PacientGroup:
    pacient = Pacients.Pacients(caseid, Mutations, filteredCNV)
    PacientsData.append(pacient.GeneticLandscape)
    
  MergedMutAndCNV = pd.concat(PacientsData, axis=0, ignore_index=True)
      
  return MergedMutAndCNV

#TODO: move this to classes?
def GenerateReports(reportDirectory, mutData, cnvData, pacientSum):

    print("Generating reports...")

    #Write filtered wrapper data to a comma separated file
    mutData.to_csv((reportDirectory + 'Mutations.csv'), index=None, header=True)
    cnvData.to_csv((reportDirectory + 'CNV.csv'), index=None, header=True)
    pacientSum.to_csv((reportDirectory + 'PacientSummary.csv'), index=None, header=True)

    print("DONE!!!!")


#If the python source file is imported as module, 
#python interpreter sets the __name__ value to module name, 
#so the if condition will return false and main method will not be executed.    
if __name__ == '__main__':
    main()
