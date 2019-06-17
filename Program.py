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
import classes.Figures as FMaker

ReportDirectory = ""
SourceFiles     = ""


def main():

  #TODO: better handling of input etc...
  #Input variables
  SourceFiles     =  input("Input project directory: ")
  ReportDirectory =  input("Input results directory: ")
  VulcanContext   =  input("Input Vulcan context: ")
    
  #Call VulcanSpot to ask for input gene list
  VulcanGenes = VulcanGeneWrapper.VulcanInputGenes(ReportDirectory, VulcanContext)

  MutationData   = ParseMAFData(SourceFiles, VulcanGenes.GeneList).MutFilteredData #filtered data from maf file
  CopyNumberData = ParseCNVData(SourceFiles, VulcanGenes.GeneList).CNVFilteredData #filtered data from cnv file  
  
  MutationAndCNV = AnalyzePacientGeneticData(MutationData, CopyNumberData) #final table summarising CNVs and mutations
  VulcanQuery    = QueryVulcanForGenes(MutationAndCNV, VulcanContext) #Treatments from VulcanSpot as json response

  #Write data
  GenerateReports(ReportDirectory, MutationData, CopyNumberData, MutationAndCNV, VulcanQuery)    

def ParseMAFData(sourceDir, vulcanGenes):

  FilteredMutations = MAFParser.ParseMAF((sourceDir + "/mutation.maf"), vulcanGenes)

  return FilteredMutations

def ParseCNVData(sourceDir, vulcanGenes):

   #Call the wrapper for metadata
    MetaData = MetaParser.MetaParser((sourceDir + "/metadata.json"))
    FilteredCNVData = CNVParser.ParseCNV((sourceDir + "/cnv.txt"), vulcanGenes)

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
def GenerateReports(reportDirectory, mutData, cnvData, pacientSum, vulcanReport):

    print("Generating reports...")

    #Write filtered wrapper data to a comma separated file
    mutData.to_csv((reportDirectory + '/Mutations.csv'), index=None, header=True)
    cnvData.to_csv((reportDirectory + '/CNV.csv'), index=None, header=True)
    pacientSum.to_csv((reportDirectory + '/PacientSummary.csv'), index=None, header=True)

    vulcanReport.AlternativeDrugTable.to_csv((reportDirectory + '/DrugTable.csv'), index=None, header=True)
    vulcanReport.DirectDrugsTable.to_csv((reportDirectory + '/DirectDrugTable.csv'), index=None, header=True)

    #testing
    print("Writing figures")
    FMaker.NewSummaryforCNV(cnvData, reportDirectory)
    FMaker.NewSummaryForMutations(mutData, reportDirectory)
    
    #CNV_figure.savefig((reportDirectory + "/cnv.png"))
    #Mut_figure.savefig((reportDirectory + "/mut.png"))

    print("DONE!!!!")


#If the python source file is imported as module, 
#python interpreter sets the __name__ value to module name, 
#so the if condition will return false and main method will not be executed.    
if __name__ == '__main__':
    main()
