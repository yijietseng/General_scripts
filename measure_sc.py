from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
import sys, glob, argparse, os
from sys import argv
init()

def get_residue_ls(pose, chain='A'):

	total = pose.total_residue()
	ligstr=''
	enzstr=''
	for i in range(1,total+1):
		resi_info = pose.pdb_info().pose2pdb(i)
		resi = resi_info.replace(' ','')
		if resi[-1] == chain:
			ligstr += resi+','
		else:
			enzstr += resi+','
	ligstr = ligstr[:-1]
	enzstr = enzstr[:-1]

	return (ligstr, enzstr)


def measrue_sc(pose, ligstr='', enzstr=''):
	sc = ShapeComplementarityFilter()

	# if a ligand is identified, it measure the ligand's sc. Otherwise, it measures global sc
	if ligstr and enzstr:
		sc.residues1(ligstr)
		sc.residues2(enzstr)
	
	sc_score = sc.score(pose)
	
	return sc_score


def ligand_sc(PDBname, chain='A', logpath=''):
	pose = pose_from_pdb(PDBname)
		
	ligstr, enzstr = get_residue_ls(pose, chain)
	sc_score = measrue_sc(pose, ligstr, enzstr)

	if logpath:
		f = open(logpath,'a+')
		f.write('PDBFN\tSC\n')
		f.write(PDBname+'\t'+str(sc_score)+'\n')
		f. close()
		print(f'log is saved to {logpath}')

	print('SC score for '+PDBname+' is:'+str(sc_score))


def multi_sc(PDBpath, logpath=''):
	'''
	Measures SC for multiple structures.
	'''

	PDB = glob.glob(PDBpath+'*.pdb')
	if logpath:
		f = open(logpath,'a+')
		f.write('PDBFN\tSC\n')

	SCls = []
	for i in PDB:
		pose = pose_from_pdb(i)
		sc_score = measrue_sc(pose)
		if logpath:
			f.write(i+'\t'+str(sc_score)+'\n')
		SCls.append(i+':'+str(sc_score))
	if logpath:
		f.close()
	print(SCls)


def single_sc(PDBname, logpath=''):
	pose = pose_from_pdb(PDBname)
	sc_score = measrue_sc(pose)
	if logpath:
		f = open(logpath,'a+')
		f.write('PDBFN\tSC\n')
		f.write(PDBname+'\t'+str(sc_score)+'\n')
		f.close()

	print('SC score for '+PDBname+' is:'+str(sc_score))

#-------------------------Below section is definitions of input options------------------------------------------

parser = argparse.ArgumentParser(description='Measure the SC score for 1 or multiple structures, it can also measure the SC for just the Ligand')
parser.add_argument("-i", 
                    '--input',
                    nargs=1,
                    type=str,
                    required=True,
                    help='input PDB file path/name')
parser.add_argument('-m',
                    '--multi', 
					action='store_true',
                    help='Measure SC for multiple structures')
parser.add_argument('-l',
                    '--ligand',
                    nargs='?',
		    		const='A',
                    type=str,
                    help='Measure SC of ligand, argument is the ChainID of the ligand, defualt = A')
parser.add_argument('-r',
                    '--record',
                    nargs=1,
		    		type=str,
                    help='Record a log file or not')
parser.add_argument('-s',
                    '--single', 
                    action='store_true',
                    help='Measure SC for a single structure')

args = parser.parse_args()

#---------------------------------------------Section ends----------------------------------------------------


def main():
	PDB = args.input[0]
	
	if args.record:
		logpath = args.record[0]
	else:
		logpath = ''

	if args.single:
		if os.path.isdir(PDB):
			raise ValueError('Please input a path to a single PDB file')
		else:
			single_sc(PDB, logpath)
	elif args.multi:
		if os.path.isdir(PDB):
			multi_sc(PDB, logpath)
		else:
			raise ValueError('Please input a path to multiple PDB file')
	
	elif args.ligand:
		if os.path.isdir(PDB):
			raise ValueError('Please input a path to a single PDB file with ligand')
		else:
			ligand_sc(PDB, args.ligand[0],logpath)
	


if __name__ == '__main__':
	main()
