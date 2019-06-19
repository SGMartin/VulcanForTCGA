#   This class parses a text file of cases IDs
#   and creates a list to be passed to the filters
#   Author: Santiago García Martín
#   Date:   19-6-2019

import os.path

class CaseFilter:

    def __init__(self, inputDirectory):

        self.IgnoreList = list()
        self.IgnoreFile = inputDirectory + '/ignorelist.txt'
        
        if os.path.exists(self.IgnoreFile):
            print("Found ignore list!")
            self.ParseIgnoreList(self.IgnoreFile)
    

    def ParseIgnoreList(self, sourceFile):

        self.fileToParse = open(sourceFile)

        for rawLine in self.fileToParse:
            self.cleanLine = rawLine.rstrip('\n')
            self.IgnoreList.append(self.cleanLine)

        self.fileToParse.close()

        print("Detected ", len(self.IgnoreFile))
          