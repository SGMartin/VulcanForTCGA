#using pandas and seaborn
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def NewSummaryforCNV(cnvFilteredData, sourceDir):

    #start new figure so that no figure is saved on top of the last one
    plt.clf()

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

    #plot the data as horizontal barplot
    CNVplot = sns.barplot(x='Frequency', y='Gene', hue='Category', data=CNVSummary,ci=None,palette=['C3','C0','k'])
    CNVplot.set(xlabel="% of cases", ylabel='CNV affected Genes')

    fig = CNVplot.get_figure()
    fig.save_figure((sourceDir + '/cnv.png'), dpi=400)


def NewSummaryForMutations(mutationFilteredData, sourceDir):

    #start new figure so that no figure is saved on top of the last one
    plt.clf()

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

    MutationPlot = MutationPlot.get_figure()
    MutationPlot.save_figure((sourceDir + '/mutations.png'))





    
