'''
Created by Yi Jie(Josh) Tseng on 11/14/2023


This script was created to convert the PD output to deuterater compatible. Because I was in a hurry,
so this script is not fully optimized. 
To use this script, simply pass the PD peptide output to this script:
python PD_2_Duet.py <file_name>

It will output a ID.csv file in the same folder as this script. 

'''


import pandas as pd
import numpy as np
import argparse



def convert(df_in, df_out, col_in:str, col_out:str, split_str=None, index=None, sec_split=None, sec_index=None):
    col_values = []
    for value in df_in[col_in]:
        if sec_split is None and sec_index is None:
            col_values.append(value.split(split_str)[index])
        else:
            col_values.append(value.split(split_str)[index].split(sec_split)[sec_index])
    df_out[col_out] = col_values

    return df_out



#-------------------------Below section is definitions of input options------------------------------------------


parser = argparse.ArgumentParser(description='''This script is developd to make the work flow of autodock_vina more smooth
                                                 Some functions can only be run in PyMOL. ''')
parser.add_argument('--peptide',
                    nargs=1,
                    type=str,
                    help='Path to csv or xlsx file of the peptide')
parser.add_argument('--protein',
                    nargs=1,
                    type=str,
                    help='Path to csv or xlsx file of the protein')
args = parser.parse_args()
#---------------------------------------------Section ends----------------------------------------------------

def main():
    
    
    if args.peptide and args.protein:
        peptide_grp_file = args.peptide[0]  # input a peptide file from PD
        protein_file = args.protein[0]  # input a protein file from PD
    else:
        print('Please input a peptide file and a protein file')

    if peptide_grp_file.rsplit('.', 1)[-1] == 'xlsx':
        df_PD_input = pd.read_excel(peptide_grp_file)
        if protein_file.rsplit('.', 1)[-1] == 'xlsx':
            df_protein = pd.read_excel(protein_file)
        else:
            df_protein = pd.read_csv(protein_file, sep=',')
    else:
        df_PD_input = pd.read_csv(peptide_grp_file, sep=',')
        if protein_file.rsplit('.', 1)[-1] == 'xlsx':
            df_protein = pd.read_excel(protein_file)
        else:
            df_protein = pd.read_csv(protein_file, sep=',')
        

    # Clean up all nan 
    df_PD_input = df_PD_input.dropna(subset=['Top Apex RT [min]'], how='any', inplace=False)
    df_PD_input.reset_index(drop=True, inplace=True)

    # For Deuterater5
    df_ID = pd.DataFrame(columns=['Sequence','Protein ID','Protein Name','Precursor Retention Time (sec)','rt_start','rt_end',
                                'rt_width','Precursor m/z','Peptide Theoretical Mass','Identification Charge','ptm','avg_ppm',
                                'start_loc','end_loc','num_peptides','num_unique','Homologous Proteins','species','gene_name',
                                'protein_existence','sequence_version','cf','neutromers_to_extract','literature_n'])

    print('[+] Converting...')
    
    # Converting sequence
    df_ID = convert(df_PD_input, df_ID, 'Annotated Sequence', 'Sequence', '.', 1)

    # Converting Protein name
    df_ID = convert(df_PD_input, df_ID, 'Master Protein Descriptions', 'Protein Name', ' OS=', 0)

    # Converting Protein ID
    for i in range(len(df_PD_input)):
        if len(df_PD_input['Master Protein Accessions'][i].split(';')) > 1:
            df_ID['Protein ID'][i] = df_PD_input['Master Protein Accessions'][i].split(';')[0]
        else:
            df_ID['Protein ID'][i] = df_PD_input['Master Protein Accessions'][i]

    # Converting Homologous Proteins
    df_ID['Homologous Proteins'] = df_PD_input['Master Protein Accessions'].str.split(';').apply(lambda x: [item.strip() for item in x])

    # Converting precursor RT
    df_ID['Precursor Retention Time (sec)'] = df_PD_input['Top Apex RT [min]'] * 60

    # Converting precursor m/z
    df_ID['Precursor m/z'] = df_PD_input['m/z [Da] (by Search Engine): CHIMERYS Identification']

    # Converting Peptide Theoretical Mass
    df_ID['Peptide Theoretical Mass'] = df_PD_input['Theo. MH+ [Da]'] - 1.01

    # Converting charge
    df_ID['Identification Charge'] = df_PD_input['Charge (by Search Engine): CHIMERYS Identification']

    # Converting # peptides and # uniques
    merged_df = df_ID.merge(df_protein[['Accession', '# Peptides', '# Unique Peptides']], left_on='Protein ID', right_on='Accession', how='left')
    df_ID['num_peptides'] = merged_df['# Peptides']
    df_ID['num_unique'] = merged_df['# Unique Peptides']

    # Converting species
    df_ID = convert(df_PD_input, df_ID, 'Master Protein Descriptions', 'species', ' OS=', 0, ' OX=', 0)

    # Converting gene_name
    df_ID = convert(df_PD_input, df_ID, 'Master Protein Descriptions', 'gene_name', ' GN=', 1, ' PE=', 0)

    # Adding protein existance
    df_ID['protein_existence'] = 1

    # Converting chemical formula
    # Loading element table
    df_element = pd.read_csv('aa_elem_compositions.tsv', sep='\t')
    composition_table = df_element.set_index('amino_acid').to_dict(orient='index')

    # Calculate chemical formulas
    
    cf_ls = []
    for seq in df_ID['Sequence']:
        overall_composition = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0, 'P': 0}
        for AA in seq:
            if AA in composition_table:
                for element, count in composition_table[AA].items():
                    if element != 'amino_acid_name':
                        overall_composition[element] += count
        cf_ls.append(''.join(f"{element}{count}" if count > 0 else '' for element, count in overall_composition.items()))

    df_ID['cf'] = cf_ls

    # Calculate neutromers_to_extract
    mass_cutoff1 = 1500
    mass_cutoff2 = 2400

    df_ID['neutromers_to_extract'] = df_ID['Peptide Theoretical Mass'].apply(lambda mass: 3 if mass < mass_cutoff1 else (4 if mass < mass_cutoff2 else 5))

    # Calculate literature_n
    # Load in labeling site and elem compositions
    df_label_site = pd.read_csv('aa_labeling_sites.tsv', sep='\t')

    # Set the study_type as the index in df_label_site
    label_table = df_label_site.set_index('study_type').to_dict(orient='index')

    n_ls = []
    for seq in df_ID.Sequence:
        n_value = 0
        for aa in seq:
            n_value += label_table['tissue'][aa]
        n_ls.append(n_value)
    df_ID['literature_n'] = n_ls

    # Assume df_ID is your original DataFrame

    # Calculate different charge states using iterative approach
    min_charge = 2
    max_charge = 6

    # Create an empty list to store DataFrames for concatenation
    dfs_to_concat = []

    # Iterate through charge states
    for charge in range(min_charge, max_charge + 1):
        # Skip if Identification Charge is the same as the iteration
        if 'Identification Charge' in df_ID and charge == df_ID['Identification Charge'].iloc[0]:
            continue
        
        # Create a DataFrame for the current charge state
        charge_state_df = df_ID.copy()

        # Assign the charge state to the new DataFrame
        charge_state_df['Identification Charge'] = charge

        # Update Precursor m/z based on Identification Charge
        charge_state_df['Precursor m/z'] = charge_state_df['Peptide Theoretical Mass'] / charge_state_df['Identification Charge']

        # Append the DataFrame to the list
        dfs_to_concat.append(charge_state_df)

    # Concatenate the original df_ID with charge_states_df
    df_ID = pd.concat([df_ID] + dfs_to_concat, ignore_index=True)

    # Save the resulting DataFrame to a CSV file
    df_ID.to_csv('ID.csv', sep=',', index=False)
    
    print('-------------\n[+] Done')

if __name__ == '__main__':
    main()