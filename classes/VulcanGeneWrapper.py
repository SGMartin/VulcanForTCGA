# This class sends a GET request to VulcanSpot to retrieve the list of input genes

import os.path
import requests

class VulcanInputGenes:

    def __init__(self, _inputDirectory, geneContext):

        print("Retrieving Vulcan input genes...")
        
        self.ReportDirectory = _inputDirectory

        self.VulcanGeneList = self.__GetGeneList(geneContext)
        self.__WriteGeneReport(self.VulcanGeneList, self.ReportDirectory)
       
        
    #These private methods are not encapsulated, that's the way Python works. But it prevents superclass overriding
    def __GetGeneList(self, typeOfCancer):

        #TODO: make this dynamic
       # self.VulcanResponse = requests.get("http://vulcanspot.org/api/genes?class=A")
        self.VulcanResponse = requests.get("http://vulcanspot.org/api/genes?context=" + typeOfCancer)
        self.SourceJson =  self.VulcanResponse.json()

        self.GeneList = list()
     
        for genes in self.SourceJson['data']:
             self.GeneList.append(genes['key'])
    
        return self.GeneList

    def __WriteGeneReport(self,GeneListToReport, _inputDir):

        self.reportDir = _inputDir + "/EligibleGenes.txt"
        self.GeneReport = open(self.reportDir, 'w')

        for lines in GeneListToReport:
            self.GeneReport.write(lines + '\n')
        self.GeneReport.close()

        print("DONE!")