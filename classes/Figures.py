#using pandas and seaborn
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def NewSummaryforCNV(cnvFilteredData, sourceDir):


    #Count occurrences of amplification and losses
    CNVSummary = cnvFilteredData.apply(pd.Series.value_counts, axis=1)[[-1,0,1]].fillna(0)

    #get pacient count, should be the same across all rows
    caseCount = CNVSummary.iloc[0,:].sum(0)
    
    #Name columns, drop neutrals
    CNVSummary.columns = ['Deletion', 'Neutral', 'Amplification']
    CNVSummary['Gene'] = cnvFilteredData['Gene']   
    CNVSummary.drop('Neutral', axis=1)

    #get relative values
    CNVSummary['Amplification'] = CNVSummary['Amplification'] / caseCount * 100
    CNVSummary['Deletion'] = CNVSummary['Deletion'] / caseCount * 100

    #melting dataframe and setting Gain/Loss as categorical 
    CNVSummary = pd.melt(CNVSummary, id_vars=['Gene'], value_vars=['Amplification', 'Deletion'], var_name='Category',
     value_name='Frequency')
    
    CNVSummary['Category'] = pd.Categorical(CNVSummary['Category'])

    #start new figure so that no figure is saved on top of the last one
    plt.figure(figsize=(15,10))
    #plot the data as horizontal barplot
    CNVplot = sns.barplot(x='Frequency', y='Gene', hue='Category', data=CNVSummary,ci=None,palette=['C3','C0','k'])
    CNVplot.set(xlabel="% of cases", ylabel='CNV affected Genes')
    CNVplot.set_xticklabels(CNVplot.get_xticklabels(),rotation=30)

    plt.savefig((sourceDir + "/cnv.png"))

  

def NewSummaryForMutations(mutationFilteredData, sourceDir):

    #start new figure so that no figure is saved on top of the last one
    plt.figure(figsize=(15,10))

    pacientCount = mutationFilteredData['CaseID'].nunique()

    #group by genes, then count the cases. This method counts affected genes only once per case

    GenesTable = mutationFilteredData.groupby('Gene')['CaseID'].nunique()
    GenesTable = pd.DataFrame(GenesTable)

    GenesTable.columns = ['Hits']
    GenesTable['Gene'] = GenesTable.index
    GenesTable['Frequency'] = GenesTable['Hits'] / pacientCount * 100

    #sort in descending order
    GenesTable = GenesTable.sort_values('Frequency', ascending=False)

    #set the plots
    MutationPlot = sns.barplot(x='Gene', y='Frequency', data=GenesTable)
    MutationPlot.set(xlabel='Affected genes', ylabel='Cases proportion with at least one SPM')
    MutationPlot.set_xticklabels(MutationPlot.get_xticklabels(),rotation=30)


    plt.savefig((sourceDir + "/mutation.png"))


def NewImpactSummary(pacientData, sourceDir):

    plt.figure(figsize=(15,10))

    figure = sns.catplot(x='Impact', kind='count', data=pacientData, palette=['C0','C3','k'])
    figure.set(xlabel='Mutations by predicted impact', ylabel='')
    plt.savefig((sourceDir + "/ImpactSummary.png"))


def NewAlterationSummary(pacientData, sourceDir):

    plt.figure(figsize=(15,10))

    pacientCount = pacientData['Pacient'].nunique()

    AlterationCount    = pacientData.groupby('Gene')['Pacient'].nunique()
    PacientAlterations = pd.DataFrame(AlterationCount)
    PacientAlterations['Gene'] = PacientAlterations.index
    PacientAlterations['Frequency'] = PacientAlterations['Pacient'] / pacientCount * 100

    PacientAlterations = PacientAlterations.sort_values('Frequency', ascending=False)

    #plot alteration frequency (no matter if cnv or somatic point mut)
    figure = sns.barplot(x='Gene',y='Frequency', data=PacientAlterations)
    figure.set(xlabel='Affected genes', ylabel='% of cases')
    figure.set_xticklabels(figure.get_xticklabels(),rotation=30)

    plt.savefig((sourceDir + "/AlterationSummary.png"))

def NewCoverageRatio(pacientData, drugData, sourceDir):

    totalPacientAlterations = pacientData['Gene'].nunique()
    DrugCoverage = drugData.groupby('AlternateDrug')['MutatedGene'].nunique()
    DrugCoverage = pd.DataFrame(DrugCoverage)
    DrugCoverage = DrugCoverage.reset_index()

    DrugCoverage['Coverage'] = DrugCoverage['MutatedGene'] / totalPacientAlterations * 100
    DrugCoverage = DrugCoverage.sort_values('Coverage', ascending=False)

    #get top 10
    DrugCoverage = DrugCoverage.sort_values('Coverage', ascending=False).head(10)

    plt.figure(figsize=(15,10))
    figure = sns.barplot(x='AlternateDrug', y='Coverage', data=DrugCoverage)
    figure.set(xlabel='Top 10 GD drugs', ylabel='% of pacient alterations covered')

    plt.savefig((sourceDir  + "/AlternateDrugCoverage.png"))


def NewGeneticDependencies(pacientData, drugData, sourceDir):

    totalAlterations = pacientData['Gene'].nunique()

    GeneticDep = drugData.groupby('Target')['MutatedGene'].nunique()
    GeneticDep = GeneticDep.reset_index()
    GeneticDep = GeneticDep.sort_values('MutatedGene', ascending=False).head(10)

    GeneticDep['Frequency'] = GeneticDep['MutatedGene'] * 100 / totalAlterations

    plt.figure(figsize=(15,10))
    figure = sns.barplot(x='Target', y='Frequency', data=GeneticDep)
    figure.set(xlabel='Top genetic dependencies', ylabel='% of alterations')

    plt.savefig((sourceDir + "/GeneticDependencies.png"))














    
