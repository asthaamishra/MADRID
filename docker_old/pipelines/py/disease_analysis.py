#!/usr/bin/python3
import re
import sys
import getopt
import os
import json
import pandas as pd
import numpy as np
from project import configs
from GSEpipelineFast import *

from rpy2.robjects import pandas2ri
pandas2ri.activate()

def breakDownEntrezs(Disease_UP):
    Disease_UP['Gene ID'] = Disease_UP['Gene ID'].str.replace('///','//')
    singleGeneNames = Disease_UP[~Disease_UP['Gene ID'].str.contains('//')].reset_index(drop=True)
    multipleGeneNames = Disease_UP[Disease_UP['Gene ID'].str.contains('//')].reset_index(drop=True)
    breaksGeneNames = pd.DataFrame(columns=['Gene ID'])
    print(singleGeneNames.shape)
    print(multipleGeneNames.shape)
    for index, row in multipleGeneNames.iterrows():
        for genename in row['Gene ID'].split('//'):
            breaksGeneNames = breaksGeneNames.append({'Gene ID': genename}, ignore_index=True)
    GeneExpressions = singleGeneNames.append(breaksGeneNames, ignore_index=True)
    return GeneExpressions

def get_entrez_id(regulated, outputFullPath, fullflag=False):
    Disease_UP = fetch_entrez_gene_id(list(regulated.index.values), input_db='Affy ID')
    Disease_UP.drop(columns=['Ensembl Gene ID'],inplace=True)
    Disease_UP.replace(to_replace='-', value=np.nan, inplace=True)
    if not fullflag:
        Disease_UP.dropna(how='any', subset=['Gene ID'], inplace=True)
        GeneExpressions = breakDownEntrezs(Disease_UP)
        # GeneExpressions.set_index('Gene ID', inplace=True)
    else: 
        GeneExpressions = Disease_UP
            
    GeneExpressions['Gene ID'].to_csv(outputFullPath, index=False)
    return GeneExpressions


def pharse_configs(inqueryFullPath):
    xl = pd.ExcelFile(inqueryFullPath)
    sheet_name = xl.sheet_names
    inqueries = pd.read_excel(inqueryFullPath, sheet_name=sheet_name, header=0)
    inqueries['Sheet1'].fillna(method='ffill',inplace=True)
    df = inqueries['Sheet1'].loc[:,['GSE ID','Samples','GPL ID','Instrument']]
    df_target = inqueries['Sheet1'].loc[:,['Samples','Experiment']]
    df_target.rename(columns={'Samples':'FileName','Experiment':'Condition'},inplace=True)
    df_target['FileName'] = df_target['FileName'].astype(str) + '.txt.gz'
    df_target['SampleNumber']= 1 + df_target.index.values
    df_target = df_target[['SampleNumber','FileName','Condition']]
    df_target['Condition'] = df_target['Condition'].str.lower()
    return df, df_target



def main(argv):
    inputfile = 'disease_transcriptomics_data_inputs.xlsx'
    targetfile = 'targets.txt'
    try:
        opts, args = getopt.getopt(argv, "hi:", ["ifile="])
    except getopt.GetoptError:
        print('python3 disease_analysis.py -i <inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 disease_analysis.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
    print('Input file is "', inputfile)
    # print('Target file is "', targetfile)
    # print(configs.rootdir)
    # filename = 'Disease_Gene_Analyzed.xlsx'
    inqueryFullPath = os.path.join(configs.rootdir, 'data', inputfile)
    querytable, df_target = pharse_configs(inqueryFullPath)
    targetdir = os.path.join(configs.datadir, targetfile)
    df_target.to_csv(targetdir, index=False, sep='\t')
    sr = querytable['GSE ID']
    gse_ids = sr[sr.str.match('GSE')].unique()
    GSE_ID = gse_ids[0]
    gseXXX = GSEproject(GSE_ID, querytable, configs.rootdir)
    for key,val in gseXXX.platforms.items():
        rawdir = os.path.join(gseXXX.genedir,key)
        print('{}:{}, {}'.format(key, val, rawdir))
        data2 = affyio.fitaffydir(rawdir, targetdir)
        data2 = ro.conversion.rpy2py(data2)

    data2['abs_logFC'] = data2['logFC'].abs()
    data2.sort_values(by='abs_logFC', ascending=False, inplace=True)
    regulated = data2[data2['abs_logFC']>=1.0]
    down_regulated = regulated[regulated['logFC']<0]
    up_regulated = regulated[regulated['logFC']>0]
    
    print("data2 before:")
    print(data2.head())
    
    print("down:")
    print(down_regulated)
    up_file = os.path.join(configs.datadir,'Disease_UP_{}.txt'.format(GSE_ID))
    down_file = os.path.join(configs.datadir,'Disease_DOWN_{}.txt'.format(GSE_ID))
    Disease_UP = get_entrez_id(up_regulated, up_file)
    Disease_UP.dropna(how='any', subset=['Gene ID'], inplace=True)
    Disease_DOWN = get_entrez_id(down_regulated, down_file)
    Disease_DOWN.dropna(how='any', subset=['Gene ID'], inplace=True)
    # Disease_UP = pd.read_csv(up_file, index_col=False, header=None)
    # Disease_UP.rename(columns={0:'Gene ID'}, inplace=True)
    # Disease_UP['Gene ID'] = Disease_UP['Gene ID'].astype(str)
    # Disease_UP = breakDownEntrezs(Disease_UP)
    # Disease_DOWN = pd.read_csv(down_file, index_col=False, header=None)
    # Disease_DOWN.rename(columns={0:'Gene ID'}, inplace=True)
    # Disease_DOWN['Gene ID'] = Disease_DOWN['Gene ID'].astype(str)
    # Disease_DOWN = breakDownEntrezs(Disease_DOWN)
    print(Disease_DOWN)
    print(Disease_UP)

    all_file = os.path.join(configs.datadir, 'Disease_FULL_{}.txt'.format(GSE_ID))
    Disease_FULL = get_entrez_id(data2, all_file, True)
    print("dis full", Disease_FULL.shape[0])
    print(Disease_FULL)
    data2.index.name = 'Affy ID'
    print("data2", data2.shape[0])
    print(data2)
    data2['ENTREZ_GENE_ID'] = Disease_FULL['Gene ID']
    #data2.dropna(how='any', subset=['ENTREZ_GENE_ID'], inplace=True)
    raw_file = os.path.join(configs.datadir, 'Raw_Fit_{}.csv'.format(GSE_ID))
    print('Raw Data saved to\n{}'.format(raw_file))
    print("data2 after:")
    print(data2.head())
    data2.to_csv(raw_file, index=False)

    files_dict = {'GSE': GSE_ID, 'UP_Reg': up_file, 'DN_Reg': down_file, 'RAW_Data': raw_file}
    files_json = os.path.join(configs.datadir, 'step2_results_files.json')
    with open(files_json, 'w') as fp:
        json.dump(files_dict, fp)
    os.remove(targetdir)


if __name__ == "__main__":
   main(sys.argv[1:])
