'''
Original Author: Nate Zuniga
Modified by: Yi Jie (Josh) Tseng on 10/01/2023

Software requirement: Python
Python module requirement: scipy, sklearn, numpy, pandas and matplotlib


Purpose:
To perform normalization based on the Aguilan et al. 2020. It includes the following steps:
1. Clean and filter data.
2. log2 transformation, normalization with average, normalization of slope.
3. Separate the samples into groups for comparison and calculate the stats of those comparison.
4. Correct the -10logPV using Benjamini_Hochberg method.
5. Prepare the data for Ontology.
6. Plot volcano plots.
7. Convert the data format to DAVID compatible.


***************************************** IMPORTANT IF RUNNING RNA DATA *****************************************
If running RNA data, the below softwares and modules are only required when you need to convert the gene ensembl
ID to UniprotID:

Softwares: chrome.exe, chromedriver.exe
Python modules: selenium, requests, bs4

Please note the version of your chrome.exe and use the corresponding version of chromedriver.exe
chromedriver.exe can be downloaded from https://googlechromelabs.github.io/chrome-for-testing/
Once chromedriver.exe is downloaded, place it in the same folder as this script.
*****************************************************************************************************************

Useage:
For proteome data:
    python MultidataSet_analyzer.py -l <path_to_data_file_or_folder>
For RNAseq data:
    python MultidataSet_analyzer.py -l <path_to_data_file_or_folder> --RNA
If no need to convert data format for ontology for DAVID
    python MultidataSet_analyzer.py -l <path_to_data_file_or_folder> --convert False

'''

from scipy import stats as st
from sklearn.impute import KNNImputer
from itertools import combinations
import numpy as np
import pandas as pd
#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import argparse, os, glob, re


# Parameter Settings
significance_cutoff = np.log10(0.05)*-10
FC = np.log2(2)


def mkdir(path):
    folderls = ['Normalized', 'Stats', 'PVFC', 'VP', 'Ontology_input']
    for folder in folderls:
        os.makedirs(f'{path}/{folder}', exist_ok=True)


def wget_UniprotID(df):
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver.chrome.service import Service
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    import numpy as np
    import time
    import pandas as pd


    df = pd.read_csv('./Normalized_RNA_Seq_counts-all_AvsC_Final_DataSet.csv', sep=',')

    # Configure Chrome WebDriver options
    options = webdriver.ChromeOptions()
    options.add_argument("--headless")  # Run Chrome in headless mode
    options.add_argument("--disable-extensions")
    options.add_argument("--disable-gpu")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--disable-software-rasterizer")
    options.add_argument("--disable-browser-side-navigation")
    options.add_argument("--disable-infobars")
    options.add_argument("--no-sandbox")
    options.add_argument("--mute-audio")  # Mute audio to suppress audio-related console messages
    options.add_argument("--log-level=3")  # Set the log level to suppress console messages


    # Specify the path to the ChromeDriver executable
    chrome_driver_path = "./chromedriver.exe"  # Replace with the actual path to chromedriver.exe (Windows) or chromedriver (Linux/Mac)

    # Create a ChromeDriver service with the specified path
    chrome_driver_service = Service(chrome_driver_path)

    print('[+] Obtaining Uniprot IDs...')

    # Create a Chrome WebDriver instance using the service
    driver = webdriver.Chrome(service=chrome_driver_service, options=options)

    #IDls = ['ENSMUSG00000029368', 'ENSMUSG00000002985', 'ENSMUSG00000064370', 'ENSMUSG00000064351', 'ENSMUSG00000064341','ENSMUSG00000064339']
    IDls = df.Accession
    UniprotIDls = []
    click_count = 0
    for ID in IDls:
        driver.get(f'https://www.uniprot.org/uniprotkb?query={ID}&facets=reviewed%3Atrue')
        # Rest of your Selenium code here
        
        if click_count != 0:
            pass
        else:
            # Find and click the <span> element with the text "Table"
            table_button = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.XPATH, "//span[text()='Table']")))
            table_button.click()

            # Find and click the <button> element with class "button primary" and text "View results"
            view_results_button = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.XPATH, "//button[@class='button primary' and text()='View results']")))
            view_results_button.click()
            click_count += 1
        # Find the <a> element with class "BqBnJ" and href attribute "/uniprotkb/P07724/entry"
        try:
            a_elements = WebDriverWait(driver, 3).until(EC.presence_of_all_elements_located((By.XPATH, "//a[contains(@class, 'BqBnJ') and contains(@href, '/uniprotkb/')]")))
            for a_element in a_elements:
                # Get the href attribute
                a_href = a_element.get_attribute("href")

                # Parse the href to extract "P07724"
                p_value = a_href.split("/")[-2]

                # Print or use the extracted value as needed
                
        except:
            p_value = np.nan

        UniprotIDls.append(p_value)  
    driver.quit()

    df = df.set_index(pd.Index(UniprotIDls))
    
    return df


def get_UniprotID(df):
    # prep UniprotID for 
    UniprotID = []
    for ID in df.index.to_list():
        UniprotID.append(ID.split('|')[0])
    df = df.set_index(pd.Index(UniprotID))
    
    return df


def get_sample_name(df):

    column_ls = df.columns.tolist()
    sample_groups = {}
    
    # Finding control group and names 
    for name in column_ls:
        match = re.search(r'Area\s+(\D*)(\d*)', name)
        if match:
            if match.group(1) not in sample_groups:
                sample_groups[match.group(1)] = [match.group(1)+match.group(2)]
            else:
                sample_groups[match.group(1)].append(match.group(1)+match.group(2))   

    return sample_groups


def data_filter(df, RNA=False):  
    print('[+] Cleaning and filtering data...')
    if RNA:
        # Cleaning data   
        # Set the 'Accession' as the index, filter out proteins with <2 unique peptides, and filter out non-area data.
        df = df.set_index('Accession')
        df = df.filter(regex='Area')
        # Sort the columns and rows and remove the values at base line
        df.sort_index(axis=1, inplace=True)
        #new_df = new_df.sort_values(by=list(new_df.columns)[0], ascending=False)
        #new_df = new_df[new_df.ne(new_df.shift())]
        #df = df.clip(lower=0)
        df = df.where(df >= 6, np.nan)

    else:  
        # Cleaning data   
        # Set the 'Accession' as the index, filter out proteins with <2 unique peptides, and filter out non-area data.
        df = df.set_index('Accession')
        df = df.loc[(df['#Unique'] >= 2),]
        df = df.loc[df['Top'] == True,]
        df = df.filter(regex='Area')
        df.sort_index(axis=1, inplace=True)

    # remove protein with >1 missing value per sample group. thresh = the number of NOT NULL values required to keep the row.
    sample_groups = get_sample_name(df)
    new_df = pd.DataFrame()
    if len(sample_groups) == 1:
        for sample_group, _ in sample_groups.items():
            # Remove duplicated sample_name
            groups = list(set(sample_groups[sample_group]))
            groups.sort()
            for group in groups:
                df.filter(regex=group)
                df_temp = df.filter(regex=group)
                df_temp = df_temp.dropna(thresh=len(df_temp.columns)-1, axis=0)
                new_df = pd.concat([new_df, df_temp],axis=1).filter(regex='Area').copy()
        new_df = new_df.dropna(thresh=len(new_df.columns)-len(groups), axis=0)
    elif len(sample_groups) > 1:
        for sample_group, _ in sample_groups.items():
            # Remove duplicated sample_name
            df.filter(regex=f'Area\s{sample_group}')
            df_temp = df.filter(regex=f'Area\s{sample_group}')
            # First get rid of the >1 missing value per sample group
            df_temp = df_temp.dropna(thresh=len(df_temp.columns)-1, axis=0)
            new_df = pd.concat([new_df, df_temp],axis=1).filter(regex='Area').copy()
            # Further clean up the combined dataframe
            new_df = new_df.dropna(thresh=len(new_df.columns)-2, axis=0)

    return new_df  
 

def separate_into_groups(df):

    sample_groups = get_sample_name(df)
    sample_df_dic = {}
    if len(sample_groups) == 1: 
        for _, sample_names in sample_groups.items():
            for sample_name in sample_names:
                sample_df = df.filter(regex=sample_name).copy()
                sample_df_dic[f'{sample_name}_DF'] = sample_df.copy()
    elif len(sample_groups) > 1:
        for exp_name, _ in sample_groups.items():
            sample_df = df.filter(regex=exp_name).copy()
            sample_df_dic[f'{exp_name}_DF'] = sample_df.copy()
    
    comparison_df_dic = {}
    # Generate all possible comparison and combine dataframes
    double_combinations = list(combinations(list(sample_df_dic.keys()),2))           

    for A_B in double_combinations:
        A, B = A_B
        df_A = sample_df_dic[A]
        df_B = sample_df_dic[B]
        df_name = 'vs'.join(A_B).replace('_DF', '').strip('_')
        comparison_df_dic[f'{df_name}'] = pd.concat([df_A, df_B],axis=1).dropna().filter(regex='Area').copy()
   
    return comparison_df_dic 
        
        
def slope_normalize(input_df, log=True, ave=True, slp=True, impute=True):
    print('[+] Normalizing...')
    df_to_normalize = input_df.copy()
    if log:
        # Creat allele comparison dataframes, log2 normalize.
        df_to_normalize = np.log2(df_to_normalize)
    if ave:
        # Normalize Data by the Average
        for sample in df_to_normalize.columns:
            df_to_normalize[sample] = df_to_normalize[sample] - np.mean(df_to_normalize[sample])
        df_to_normalize['protein_avg'] = np.mean(df_to_normalize,axis=1)
    if slp:
        # Normalize Data by the Slope
        for sample in df_to_normalize.columns:
            x= df_to_normalize['protein_avg']
            y= df_to_normalize[sample]
            mask = ~np.isnan(x) & ~np.isnan(y)
            slope = st.linregress(x[mask], y[mask])[0]
            df_to_normalize[sample] = df_to_normalize[sample]/slope
        normalized_DF = df_to_normalize.filter(regex='Area')
    if impute:
        # Impute missing values using KNN imputer
        imputer = KNNImputer(n_neighbors=2)
        normalized_DF = pd.DataFrame(imputer.fit_transform(normalized_DF),columns=normalized_DF.columns,index=normalized_DF.index)
        
    return normalized_DF

   
def stats_calc(input_df, comparison=None):
    
    print('[+] Calculating stats...')
    def f_test(group1, group2):        
        # Calculate the sample variances
        var1 = np.var(group1, ddof=1)
        var2 = np.var(group2, ddof=1)
        
        # Calculate the degrees of freedom
        dof1 = len(group1) - 1
        dof2 = len(group2) - 1
        
        # Calculate the F-statistic
        if var1 > var2:
            f_value = var1 / var2
            # Calculate the p-value
            p_value = 1 - st.f.cdf(f_value, dof1, dof2)
        else:
            f_value = var2 / var1
            # Calculate the p-value
            p_value = st.f.cdf(f_value, dof2, dof1)
            
        return p_value
    
          
    if comparison:
        allele_a, allele_b = comparison.split('vs')
    else:
        raise ValueError('Please specify what needs to be compared') 
      
    comparison_DF_b = input_df.filter(regex=f'Accession|{allele_b}').copy()
    comparison_DF_a = input_df.filter(regex=f'Accession|{allele_a}').copy()
    
    # Create dictionaries to store PVs and FCs
    comparison_p_val_dict = {}
    comparison_FC_dict = {}
    
    for protein in input_df.index:
        #ApoEb vs ApoEa P-value and FC
        a = comparison_DF_b.loc[protein,]
        b = comparison_DF_a.loc[protein,]
        
        f_PV = f_test(a, b)
    
        if f_PV > 0.05:
            comparison_p_val = st.ttest_ind(a, b, equal_var=True, alternative='two-sided',nan_policy='omit')[1]
            comparison_p_val_dict[protein] = comparison_p_val
            
        elif f_PV < 0.05:
            comparison_p_val= st.ttest_ind(a,b, equal_var=False, alternative='two-sided',nan_policy='omit')[1]
            comparison_p_val_dict[protein] = comparison_p_val
        
        # Calculate FC. FC = mean_area_treament/mean_area_control (is a substraction becuase of log2 transformation)
        comparison_FC = (np.mean(b)) - (np.mean(a))
        comparison_FC_dict[protein] = comparison_FC
    comparison_df_PV = pd.DataFrame.from_dict(comparison_p_val_dict.items()).rename(columns={0 : 'Accession',1 : f'P_val'}).set_index('Accession')
    comparison_df_FC = pd.DataFrame.from_dict(comparison_FC_dict.items()).rename(columns={0 : 'Accession',1 : f'Fold_change'}).set_index('Accession')
    
    
    result_df = pd.concat([input_df, comparison_df_PV, comparison_df_FC], axis=1)
    result_df[' -log10(PV)'] = (np.log10(result_df['P_val']) * -10)
        
    return(result_df)
    

def performBH_correction(input_df):
    print('[+] Performing Benjamini_Hochberg correction on p-values...')    
    input_df['Benjamini_Hochberg_pval'] = None
    input_df = input_df.reset_index()
    #sort cleaned_df by pvalue jc mod
    input_df = input_df.sort_values(by='P_val')
    input_df = input_df.reset_index(drop=True) #sort keeps the origional index value so you need to re-index to use it in the BH calc
    #calculate benjamini_hochberg correction as (rank/total numer of tests)*probability of false positive jc mod
    total_rows = len(input_df.index)
    
    for row in input_df.itertuples():
        BH_pval = ((row.Index+1)/total_rows)*0.25
        input_df.at[row.Index, 'Benjamini_Hochberg_pval'] = BH_pval
    input_df[' -10*log10(BH)'] = np.log10(input_df['Benjamini_Hochberg_pval'].astype(float)) * -10
    input_df = input_df.set_index('Accession')
    
    return input_df


def plotVP(input_df, fileName=None):

    df = input_df
    #sns.set_style('white')
    plt.figure(figsize=(10, 10))
    #sns.set(style='white', context='talk')
    plt.axvline(x=FC, ymin=0, ymax=1, linestyle=':',color='gray')
    plt.axvline(x=-FC, ymin=0, ymax=1, linestyle=':',color='gray')
    plt.axhline(y=significance_cutoff, xmin=0, xmax=1, linestyle=':',color='gray')
    plt.xlabel("Fold_change")
    plt.ylabel("-10*log10(BH)")
    plt.title(f'P-value vs. Fold Change',fontsize=30)
    plt.xlim(-6, 6)

    #df['cutoff'] = ''

    #df.loc[((df['Fold_change'] < FC) | (df['Fold_change'] > -FC)) & (df[' -10*log10(BH)'] < significance_cutoff),['cutoff']] = 100
    #df.loc[((df['Fold_change'] < FC) | (df['Fold_change'] > -FC)) & (df[' -10*log10(BH)'] >= significance_cutoff),['cutoff']] = 100
    #df.loc[((df['Fold_change'] >= FC)) & (df[' -10*log10(BH)'] >= significance_cutoff),['cutoff']] = 250
    #df.loc[((df['Fold_change'] <= -FC)) & (df[' -10*log10(BH)'] >= significance_cutoff),['cutoff']] = 15
    
    significant = ((df['Fold_change'] <= -FC) & (df[' -10*log10(BH)'] >= significance_cutoff) | (df['Fold_change'] >= FC) & (df[' -10*log10(BH)'] >= significance_cutoff))
    colors = ['black', 'red', 'red']
    cmap = ListedColormap(colors)
    plt.scatter(df["Fold_change"], df[' -10*log10(BH)'], c=significant, cmap=cmap, s=20)
    plt.savefig(f'{fileName}', format="svg", transparent=True, dpi=300)


#-------------------------Below section is definitions of input options------------------------------------------

parser = argparse.ArgumentParser(description='''This script is developd to make the work flow of autodock_vina more smooth
                                                 Some functions can only be run in PyMOL. ''')
parser.add_argument('-l',
                    '--location',
                    nargs=1,
                    type=str,
                    help='Path to csv data files')
parser.add_argument('--RNA', 
                    action='store_true',
                    default=False,
                    help='normalize data from RNA seq')
parser.add_argument('--convert', 
                    action='store_false',
                    default=True,
                    help='Pass False if you do not want to convert to the format for ontology. For RNA, this is time consuming, becasue it will search Uniprot and convert gene ensembl ID to UniprotID')
args = parser.parse_args()

#---------------------------------------------Section ends----------------------------------------------------

def main():
    if not args.location:
        raise ValueError('Please input a valid path, ex: C:/Desktop/')
    # This is to deal with spaces in the file path
    if len(args.location) > 1:
        where_are_the_files = ' '.join(args.location)
    if '\\' in where_are_the_files:
            where_are_the_files = where_are_the_files.replace('\\', '/')  
    # NORMALIZATION AND STATS
    where_are_the_files = where_are_the_files.rsplit('/', 1)[0]
    mkdir(where_are_the_files)
    
    all_file_names = glob.glob(f'{where_are_the_files}/*.csv')

    for dataSet in all_file_names:
        # Read in 
        print(f'[+] Reading {dataSet}...')
        dataDF = pd.read_csv(dataSet, delimiter=',')
        filename = f'{dataSet.rsplit(".",1)[0]}'.rsplit('/',1)[1]

        # Clean and filter data set and normalize and save a plot for visualization
        if args.RNA:
            cleaned_df = data_filter(dataDF, RNA=True)
            norm_df =  slope_normalize(cleaned_df, log=False)
        else:
            cleaned_df = data_filter(dataDF)
            norm_df =  slope_normalize(cleaned_df)
        #cleaned_df.to_csv(f'{where_are_the_files}/cleaned.csv', sep=',')
        
        norm_df.to_csv(f'{dataSet.rsplit("/",1)[0]}/Normalized/{filename}_NormDF.csv', sep=',')
        print('[+] Plotting normalized distribution...')
        density_plot = norm_df.plot.density(figsize=(10, 10)) # This will not show plot when running under WSL
        density_plot.get_figure().savefig(f'{dataSet.rsplit("/",1)[0]}/Normalized/{filename}_normalized_result.png', dpi = 300)
        
        # Generate all comparison groups
        comparison_dataset_dic = separate_into_groups(norm_df)
        #print(comparison_dataset_dic)
        
        for key, df in comparison_dataset_dic.items():

            # Calculate the significance and fold change
            significant_df = stats_calc(df, key)
            significant_df.to_csv(f'{dataSet.rsplit("/",1)[0]}/Stats/{filename}_{key}_stats.csv', sep=',')
            
            # Perform the Benjamini_Hochberg_pval correction
            corrected_df = performBH_correction(significant_df)
            corrected_df.to_csv(f'{dataSet.rsplit("/",1)[0]}/PVFC/{filename}_{key}_PVFC.csv', sep=',')
            
            if args.convert:
                if args.RNA:
                    # Prep for Ontology
                    df_onto = get_UniprotID(corrected_df)
                    #df_onto = wget_UniprotID(corrected_df)
                else:
                    # Prep for Ontology
                    df_onto = get_UniprotID(corrected_df)
                df_positive = df_onto[(df_onto['Fold_change'] >= FC) & (df_onto[' -10*log10(BH)'] >= significance_cutoff)] 
                df_negative = df_onto[(df_onto['Fold_change'] <= -FC) & (df_onto[' -10*log10(BH)'] >= significance_cutoff)] 
                df_onto.to_csv(f'{dataSet.rsplit("/",1)[0]}/Ontology_input/{filename}_{key}_bg.csv', sep=',')
                df_positive.to_csv(f'{dataSet.rsplit("/",1)[0]}/Ontology_input/{filename}_{key}_pos.csv', sep=',')
                df_negative.to_csv(f'{dataSet.rsplit("/",1)[0]}/Ontology_input/{filename}_{key}_neg.csv', sep=',')
            
            # plot volcano plots
            print(f'[+] Plotting volcano plot for {key}...')
            plotVP(corrected_df, f'{dataSet.rsplit("/",1)[0]}/VP/{filename}_{key}_VP.svg')
        print('-----------------\n[+] Done')
if __name__ == '__main__':
    main()
