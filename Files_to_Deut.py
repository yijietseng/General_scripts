'''
Original Author: Nate Zuniga
Modified by Yi Jie(Josh) Tseng on 11/1/2023

This script is to convert the PEAKS11 output to PEAKS10 format which is compatible to deuterater.

Requirement:
Make sure to have the features.csv, protein.csv and peptide.csv from PEAKS11 output


Usage:

python Files_to_Deut.py -l <csv_file_location>

'''



import os, glob, argparse
import pandas as pd

def mkdir(path):
    folderls = ['DeuteraterReady']
    for folder in folderls:
        os.makedirs(f'{path}/{folder}', exist_ok=True)

#-------------------------Below section is definitions of input options------------------------------------------

parser = argparse.ArgumentParser(description='''This script is developd to make the work flow of autodock_vina more smooth
                                                 Some functions can only be run in PyMOL. ''')
parser.add_argument('-l',
                    '--location',
                    nargs="*",
                    type=str,
                    help='Path to the folder where the csv data files are')
args = parser.parse_args()

#---------------------------------------------Section ends-------------------------------------------------------

def main():
    if args.location:

        where_are_the_files = f'{args.location[0]}'
        if '\\' in args.location[0]:
            where_are_the_files = args.location[0].replace('\\', '/')
    else:
        raise ValueError('Please input a valid path, ex: C:/Desktop/')
    
    
    mkdir(where_are_the_files)

    all_file_names = glob.glob(f'{where_are_the_files}/*.csv')
    
    print('[+] Converting...')
    
    full_features_DF = pd.DataFrame()
    for file_name in all_file_names:
        # concat all feature files
        if 'features' in file_name:
            peaks11_protPep_DF = pd.read_csv(file_name)
            full_features_DF = pd.concat([full_features_DF,peaks11_protPep_DF])
            full_features_DF.to_csv(f'{where_are_the_files}/DeuteraterReady/peaks11_feature_file.csv',index=False)
        # Rename the column in protein file
        elif 'proteins.csv' in file_name:
            peaks11_proteins_DF = pd.read_csv(file_name, sep=',')
            peaks11_proteins_DF = peaks11_proteins_DF.rename(columns={'Average Mass':'Avg. Mass'})
            peaks11_proteins_DF.to_csv(f'{where_are_the_files}/DeuteraterReady/peaks11_proteins_file.csv',index=False)
        # Convert the rest of the files
        else:
            peaks11_prot_pep_DF = pd.read_csv(file_name, sep=',')
            peaks11_prot_pep_DF = peaks11_prot_pep_DF.rename(columns={'Accession':'Protein Accession'})
            peaks11_prot_pep_DF['Peptide'] = "X." + peaks11_prot_pep_DF['Peptide'] + ".X"
            peaks11_prot_pep_DF.to_csv(f'{where_are_the_files}/DeuteraterReady/peaks11_prot_pept_file.csv',index=False)
    
    print('----------------\n[+] Done')
if __name__ == '__main__':
    main()