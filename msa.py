'''
This script is to graph a MSA figure. This script requires a master fasta file with multiple fasta formatted sequences.
It can also combine multiple fasta file into a master file by calling a -c option, and if the -g option is called, it 
will automatically graph the figure. Available options are listed below. 

usage: msa.py [-h] [-c COMBINE] [-f FIXEDRESI [FIXEDRESI ...]] [-g] -i INPUT [-o OUTPUT]

Prep fasta files and generate a MSA figure

optional arguments:
  -h, --help            show this help message and exit
  -c COMBINE, --combine COMBINE
                        combine fasta files into a master file, takes an
                        argument for the output file
  -f FIXEDRESI [FIXEDRESI ...], --fixedresi FIXEDRESI [FIXEDRESI ...]
                        Residue number that were fixed during design
  -g, --graph           Graph the MSA
  -i INPUT, --input INPUT
                        input fasta file path/name
  -o OUTPUT, --output OUTPUT
                        output file path/name

'''
from pymsaviz import MsaViz
import argparse, os, glob


def graph(master_fasta:str ,outputFN:str, fixed_resi):
    msa_file = master_fasta
    mv = MsaViz(msa_file, wrap_length=60)

    # Extract MSA positions less than 50% consensus identity
    probability = []
    ident_list = mv._get_consensus_identity_list()
    for pos, ident in enumerate(ident_list, 1):
        if ident >= 75 :
                probability.append(pos)

    # Add markers
    mv.add_markers(fixed_resi)
    mv.set_highlight_pos(probability)


    mv.savefig(outputFN, dpi=300)

def combine(path, output_name):

    all_fasta = glob.glob(f'{path}/*.fasta')
    for fastafile in all_fasta:
        with open(fastafile, 'r') as fin:
            seq = fin.read()
        with open(output_name, 'a') as fout:
            fout.write('\n'+seq)
    # removing the leading space
    with open(output_name, 'r') as fin:
        seq = fin.read()
    with open(output_name, 'w') as fout:
        fout.write(seq[1:])

#-------------------------Below section is definitions of input options------------------------------------------

parser = argparse.ArgumentParser(description='Prep fasta files and generate a MSA figure')
parser.add_argument('-c',
                    '--combine', 
                    nargs=1,
                    type=str,
                    help='combine fasta files into a master file, takes an argument for the output file')
parser.add_argument('-f',
                    '--fixedresi',
                    nargs='+',
                    type=int,
                    help='Residue number that were fixed during design')
parser.add_argument('-g',
                    '--graph',
                    action='store_true',
                    help='Graph the MSA')
parser.add_argument("-i", 
                    '--input',
                    nargs=1,
                    type=str,
                    required=True,
                    help='input fasta file path/name')
parser.add_argument('-o',
                    '--output', 
                    nargs=1,
                    type=str,
                    help='output file path/name')

args = parser.parse_args()

#---------------------------------------------Section ends----------------------------------------------------

def main():
    # store input path
    fasta_input = args.input[0]
    
    if args.fixedresi:
        fixed_ls = args.fixedresi
    else:
        fixed_ls = []
    
    if args.combine:
        # Check if the input is a file or a directory
        if os.path.isdir(fasta_input):
            # If the last character is a / then remove it.
            if fasta_input[-1] == '/':
                fasta_input = fasta_input[:-1]
            combine(fasta_input, args.combine[0])
            fasta_input = args.combine[0]
        else:
            raise ValueError('Please input a path to the directory where all fasta files are in.')
    else:
        # Check if the input is a file or a directory
        if os.path.isdir(fasta_input):
            raise ValueError('Please input the name for the master fasta file')
        else:
            pass

    if args.graph:
        try:
            outputfile = args.output[0]
        except:
            raise ValueError('Missing output file name.')
        
        # execute the multialignment and graph
        if fixed_ls:
            graph(fasta_input, outputfile, fixed_resi=fixed_ls)
        else:
            graph(fasta_input,outputfile)
            print('test')


if __name__ == '__main__':
    main()