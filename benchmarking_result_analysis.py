import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_original = pd.read_csv("benchmarkingFile4.tsv", sep="\t") #read TSV file into pandas DataFrame
#df_original = df_original[573:1155] # tot 1156 momenteel
#df_original.columns = ['timeOfJob', 'inputFile', 'repeat', 'errorRateInInference','siteSpecificInference', 'siteSpecificSimulation', 'errorRateInSimulation',  'lRef', 'leaves', '||', 'runtime', 'LK', 'RF', 'normalisedRF', 'foundBranches', 'missedBranches', 'notFoundBranches', 'RFL', 'totalBranchLength', 'totalBranchLengthTrue']

for y_variable in [ 'RF', 'RFL', 'LK', 'runtime']: # used for site specific errors.
    df = df_original
    group_variable = 'errorRateInSimulation'#errorRateInEstimation
    x_variable = 'leaves'
    # df = df[df['errorRateInSimulation'] == trueError]
    df = df.groupby([x_variable, group_variable, 'siteSpecificInference'],as_index=False).agg([np.mean, np.std, np.size])
    df = df.reset_index()

    means = pd.DataFrame({'simulated site specific error rate = 0.0001' :
                              list(df[[x and y for x, y in zip(df[group_variable] == 0.00010, df['siteSpecificInference']==False)]][y_variable]['mean']),
                       'simulated site specific error rate = 0.0001, site specific inference ':
                           list(df[[x and y for x, y in zip(df[group_variable] == 0.00010, df['siteSpecificInference'] == True)]][y_variable]['mean']),
                          'simulated site specific error rate = 0.0005' :
                              list(df[[x and y for x, y in zip(df[group_variable] == 0.00050, df['siteSpecificInference']==False)]][y_variable]['mean']),
                       'simulated site specific error rate = 0.0005, site specific inference ':
                           list(df[[x and y for x, y in zip(df[group_variable] == 0.00050, df['siteSpecificInference'] == True)]][y_variable]['mean']),
                          },index=list(set(df[x_variable])))
    stds = pd.DataFrame({'simulated site specific error rate = 0.0001' :
                              list(df[[x and y for x, y in zip(df[group_variable] == 0.00010, df['siteSpecificInference']==False)]][y_variable]['std']),
                       'simulated site specific error rate = 0.0001, site specific inference ':
                           list(df[[x and y for x, y in zip(df[group_variable] == 0.00010, df['siteSpecificInference'] == True)]][y_variable]['std']),
                          'simulated site specific error rate = 0.0005' :
                              list(df[[x and y for x, y in zip(df[group_variable] == 0.00050, df['siteSpecificInference']==False)]][y_variable]['std']),
                       'simulated site specific error rate = 0.0005, site specific inference ':
                           list(df[[x and y for x, y in zip(df[group_variable] == 0.00050, df['siteSpecificInference'] == True)]][y_variable]['std']),
                          },index=list(set(df[x_variable])))
    ax = means.plot.bar(rot=0, yerr=stds)
    plt.ylabel(y_variable)
    plt.xlabel('Number of samples')
    plt.legend(title=group_variable)
    plt.show()

#comparing the performance for different inference errors
for trueError in set(df_original['errorRateInSimulation']):
    df = df_original
    group_variable = 'errorRateInInference'#errorRateInEstimation
    x_variable = 'leaves'
    df = df[df['errorRateInSimulation'] == trueError]
    df = df.groupby([x_variable, group_variable],as_index=False).agg([np.mean, np.std])
    df = df.reset_index()
    y_variable = 'RFL'

    means = pd.DataFrame({'error rate = 0' + ' (true error rate)' if 0 == trueError else 'error rate = 0': list(df[df[group_variable] == 0][y_variable]['mean']),
                       'error rate = 0.0001' + ' (true error rate)' if 0.0001 == trueError else 'error rate = 0.0001': list(df[df[group_variable] == 0.00010][y_variable]['mean']),
                       'error rate = 0.0005' + ' (true error rate)' if 0.0005 == trueError else 'error rate = 0.0005': list(df[df[group_variable] == 0.00050][y_variable]['mean'])
                           },index=list(df[df[group_variable] == 0][x_variable])) #set(df[group_variable])
    stds = pd.DataFrame({'error rate = 0'+ ' (true error rate)' if 0 == trueError else 'error rate = 0': list(df[df[group_variable] == 0][y_variable]['std']),
                       'error rate = 0.0001'+ ' (true error rate)' if 0.0001 == trueError else 'error rate = 0.0001':  list(df[df[group_variable] == 0.00010][y_variable]['std']),
                       'error rate = 0.0005' + ' (true error rate)' if 0.0005 == trueError else 'error rate = 0.0005': list(df[df[group_variable] == 0.00050][y_variable]['std'])
                           },index=list(df[df[group_variable] == 0][x_variable]))
    ax = means.plot.bar(rot=0, yerr=stds)
    plt.ylabel(y_variable)
    plt.xlabel('Number of samples')
    plt.legend(title=group_variable)
    plt.show()

#comparing the performance for different simulation errors
for y_variable in [ 'RF', 'RFL', 'LK']:
    simulationError = 0
    df = df_original
    group_variable = 'errorRateInSimulation'
    x_variable = 'leaves'
    df = df[df['errorRateInInference'] == simulationError]
    df = df.groupby([x_variable, group_variable],as_index=False).agg([np.mean, np.std, np.size])
    df = df.reset_index()

    means = pd.DataFrame({'error rate = 0' + ' (estimation error rate)' if 0 == simulationError else 'error rate = 0': list(df[df[group_variable] == 0][y_variable]['mean']),
                       'error rate = 0.0001' + ' (estimation error rate)' if 0.0001 == simulationError else 'error rate = 0.0001': list(df[df[group_variable] == 0.00010][y_variable]['mean']),
                       'error rate = 0.0005' + ' (estimation error rate)' if 0.0005 == simulationError else 'error rate = 0.0005': list(df[df[group_variable] == 0.00050][y_variable]['mean'])
                           },index=list(df[df[group_variable] == 0][x_variable])) #set(df[group_variable])
    stds = pd.DataFrame({'error rate = 0'+ ' (estimation rate)' if 0 == simulationError else 'error rate = 0': list(df[df[group_variable] == 0][y_variable]['std']),
                       'error rate = 0.0001'+ ' (estimation error rate)' if 0.0001 == simulationError else 'error rate = 0.0001':  list(df[df[group_variable] == 0.00010][y_variable]['std']),
                       'error rate = 0.0005' + ' (estimation error rate)' if 0.0005 == simulationError else 'error rate = 0.0005': list(df[df[group_variable] == 0.00050][y_variable]['std'])
                           },index=list(df[df[group_variable] == 0][x_variable]))
    ax = means.plot.bar(rot=0, yerr=stds)
    plt.ylabel(y_variable)
    plt.xlabel('Number of samples')
    plt.legend(title=group_variable)
    plt.show()
