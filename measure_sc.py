from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import PerResidueSasaMetric
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.core.pose import get_chain_id_from_chain
import sys, glob, argparse, os, operator
from sys import argv
init()

def get_chain_resi_ls(pose, chains=['A']):
	# getting rosetta chainIDs and total residues
	chainIDls = []
	for i in chains:
		chainIDs = get_chain_id_from_chain(i, pose)
		chainIDls.append(chainIDs)

	res_str_ls = []
	for chainID in chainIDls:
		res_str = ''
		chain_begin = pose.chain_begin(chainID)
		chain_end = pose.chain_end(chainID)
		for i in range(chain_begin, chain_end+1):
			resID = pose.pdb_info().pose2pdb(i).replace(' ', '')
			res_str += resID+','
		res_str = res_str[:-1]
		res_str_ls.append(res_str)


	return res_str_ls
	
	







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


def measure_sc(pose, ligstr='', enzstr=''):
	sc = ShapeComplementarityFilter()

	# if a ligand is identified, it measure the ligand's sc. Otherwise, it measures global sc
	if ligstr and enzstr:
		sc.residues1(ligstr)
		sc.residues2(enzstr)
	
	sc_score = sc.score(pose)
	
	return sc_score


def ligand_sc(PDBname, chain='A', logpath='', buried=False):
	pose = pose_from_pdb(PDBname)
		
	ligstr, enzstr = get_residue_ls(pose, chain)
	sc_score = measrue_sc(pose, ligstr, enzstr)
	if buried:
		buried_sa = measure_buried_SA(pose)
		if logpath:
			f = open(logpath,'a+')
			f.write('PDBFN\tSC\tBuried_sa\n')
			f.write(f'{PDBname}\t{sc_score}\t{buried_sa}\n')
			f. close()
			print(f'log is saved to {logpath}')

		print(f'SC score for {PDBname} is: {sc_score}')
		print(f'Buried surface area of the ligand is: {buried_sa}')
	else:
		if logpath:
			f = open(logpath,'a+')
			f.write('PDBFN\tSC\n')
			f.write(f'{PDBname}\t{sc_score}\n')
			f. close()
			print(f'log is saved to {logpath}')

		print(f'SC score for {PDBname} is: {sc_score}')


def multi_sc(PDBpath, logpath='', buried=False, focus=[]):
	'''
	Measures SC for multiple structures.
	'''

	PDB = glob.glob(PDBpath+'*.pdb')
	if logpath:
		if buried:
			f = open(logpath,'a+')
			f.write('PDBFN\tSC\tBuried_SA\n')
		else:
			f = open(logpath,'a+')
			f.write('PDBFN\tSC\n')

	SCls = []
	for i in PDB:
		pose = pose_from_pdb(i)
		sc_score = measure_sc(pose)
		if buried:
			buried_sa = measure_buried_SA(pose)
			if logpath:
				f.write(f'{i}\t{sc_score}\t{buried_sa}\n')

			SCls.append(f'{i} : {sc_score}, {buried_sa}')
			if logpath:
				f.close()
		else:
			if logpath:
				f.write(f'{i}\t{sc_score}\n')

			SCls.append(f'{i} : {sc_score}')
			if logpath:
				f.close()
		
	print(SCls)


def single_sc(PDBname, logpath='', buried=False, focus=[]):
	pose = pose_from_pdb(PDBname)
	sc_score = measure_sc(pose)
	
	if buried:
		buried_sa = measure_buried_SA(pose)
		if logpath:
			f = open(logpath,'a+')
			f.write('PDBFN\tSC\tBuried_SA\n')
			f.write(f'{PDBname}\t{str(sc_score)}\t{str(buried_sa)}\n')
			f.close()
		print(f'SC score for {PDBname} is: {sc_score}')
		print(f'Calculated buried surface area is: {buried_sa}')
	else:
		if logpath:
			f = open(logpath,'a+')
			f.write('PDBFN\tSC\n')
			f.write(f'{PDBname}\t{str(sc_score)}\n')
			f.close()
		print('SC score for '+PDBname+' is:'+str(sc_score))


def measure_buried_SA(pose):
	sasa_metric = PerResidueSasaMetric()
	# Creating a translated pose
	pose_t = Pose()
	pose_t.assign(pose)
	pose_t = RBTranslate(pose_t)

	# Getting the SASA for each residue
	resi_sasa1 = sasa_metric.calculate(pose)
	resi_sasa2 = sasa_metric.calculate(pose_t)
	# Formatting the SASA
	resi_sasa1 = sorted(resi_sasa1.items(), key=operator.itemgetter(1), reverse=False)
	resi_sasa2 = sorted(resi_sasa2.items(), key=operator.itemgetter(1), reverse=False)
	relist1, salist1 = zip(*resi_sasa1)
	relist1, salist2 = zip(*resi_sasa2)

	# Calculating buried surface area
	SA1 = sum(salist1)
	SA2 = sum(salist2)
	return SA2 - SA1


def RBTranslate(pose, step_size=500):
	# perform rigid body translation
	trans_mover = RigidBodyTransMover(pose, 1) 
	trans_mover.step_size(step_size) # 500 Ang is the default
	trans_mover.apply(pose)
	return pose


#-------------------------Below section is definitions of input options------------------------------------------

parser = argparse.ArgumentParser(description='Measure the SC score for 1 or multiple structures, it can also measure the SC for just the Ligand')
parser.add_argument('-b',
		    		'--buried', 
                    action='store_true',
                    help='Calculate buried surface area')
parser.add_argument('-f',
		    		'--focus',
					nargs='+',
				    type=str,
                    help='Arguments=chain ID(s). Focus mode is to calculate SC at specific interface(s).')
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
			if args.buried:
				single_sc(PDB, logpath, buried=True)
			else:
				single_sc(PDB, logpath)
	elif args.multi:
		if os.path.isdir(PDB):
			if args.buried:
				multi_sc(PDB, logpath, buried=True)
			else:
				multi_sc(PDB, logpath)
		else:
			raise ValueError('Please input a path to multiple PDB file')
	
	elif args.ligand:
		if os.path.isdir(PDB):
			raise ValueError('Please input a path to a single PDB file with ligand')
		else:
			if args.buried:
				ligand_sc(PDB, args.ligand[0],logpath, buried=True)
			else:
				ligand_sc(PDB, args.ligand[0],logpath)


if __name__ == '__main__':
	main()
