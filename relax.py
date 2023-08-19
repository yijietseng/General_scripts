#-------------------------Set up Pyrosetta------------------------------
import time
from sys import argv
from pyrosetta.rosetta import *
from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.rosetta.core.scoring import atom_pair_constraint
from pyrosetta.rosetta.protocols.relax import FastRelax
import itertools, os, argparse, glob
import multiprocessing as mp
import numpy as np

#init('-beta_nov16 True')
init()

def format_outputFN(pdbname):
	if '_R' in pdbname:
		pdbname = pdbname.split('.pdb')[0][:-2]
		R_prefix = ''
	else:
		pdbname = pdbname.split('.pdb')[0]
		R_prefix = '_R'

	for i in itertools.count(start=1):
		if not os.path.exists(pdbname+R_prefix+"{:02d}".format(i)+'.pdb'):
			outputFN = pdbname+R_prefix+"{:02d}".format(i)+'.pdb'
			print(f'\noutput file is set to {outputFN}\n')
			break

	return outputFN


def set_cst(sf, pose, cstfile):
	constraint_file = cstfile
	constraint = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
	constraint.constraint_file(constraint_file)
	constraint.add_constraints(True)
	constraint.apply(pose)
	sf.set_weight(atom_pair_constraint, 1)


def gen_cst(poseName, filename, weight=0.25):

	pose = pose_from_pdb(poseName)
	
	# Bind the CA of the backbone
	# List out the residue indices
	bkbone = np.arange(1,pose.total_residue()+1,dtype=int)

	filepath = f'{filename}.cst'

	c1 = open(filepath, "w")
	# Loop through and only bind together residues that are at least 7 away in index, which typically follows the sequence
	c1.write('\n# Connect the CA of the backbone together\n\n')
	try:
		for index in bkbone:
			for ind in bkbone:
				if (index-ind > 7):
					dVec= pose.residue(index).xyz("CA")-pose.residue(ind).xyz("CA")
					# We will see if it allows me to remove the spaces from the output of pose2pdb this easily
					pdbind = pose.pdb_info().pose2pdb(ind).replace(' ','')
					pdbindex = pose.pdb_info().pose2pdb(index).replace(' ','')
					c1.write(f'AtomPair CA {pdbind} CA {pdbindex} FLAT_HARMONIC {dVec.norm()} {weight} {weight}\n')
	except:
		print('some residue has no CA')		
	print(f'Constraint file: {filepath} written')
		

def set_movemap(set_bb, set_chi, set_jump):
	movemap = MoveMap()
	movemap.set_bb(set_bb) 
	movemap.set_chi(set_chi)
	movemap.set_jump(set_jump)
	print(movemap)
	return movemap

def RunRelax(sf, mmp, pose, apply=True):
	relax = FastRelax()
	relax.set_scorefxn(sf)

	if apply:
		relax.apply(pose)
		return pose
	else:
		print('\npassed FastRelax Mover\n')
	

def cal_wall_time(start, end):
	d = end-start
	if d >= 3600:
		print ('\n'*2, 'The time spent on relaxing with constraints was:', d/3600, 'hours')
	elif 60 <= d < 3600:
		print ('\n'*2, 'The time spent on relaxing with constraints was:', d/60, 'minutes')
	elif d < 60:
		print ('\n'*2, 'The time spent on relaxing with constraints was:', d, 'seconds')
	return d


def RunProcess(PDB, bb, chi, jump, cst=''):

	pose = pose_from_pdb(PDB)
	outputFN = format_outputFN(PDB)
	#scorefxn = create_score_function('beta_nov16')
	scorefxn = get_fa_scorefxn()

	if cst:
		set_cst(scorefxn, pose, cst)
	# Set movemap
	mmp = set_movemap(set_bb=bb, set_chi=chi, set_jump=jump)
	
	# Run relax
	score_before = scorefxn(pose)
	pose_relaxed = RunRelax(scorefxn, mmp, pose, apply=True)
	score_after = scorefxn(pose_relaxed)
	
	# output
	print ('\n'*2, f'The initial score of {outputFN} structure is:', '\n', score_before)
	print ('The score after relaxing is:', '\n', score_after)
	print ('\n','The distribution of scores', '\n')
	scorefxn.show(pose_relaxed)

	pose.dump_pdb(outputFN)

	return f'{outputFN}\t{score_after}'

#-------------------------Below section is definitions of input options------------------------------------------

parser = argparse.ArgumentParser(description='Relax inputed structure')
parser.add_argument('-c',
					'--constraint',
					nargs=1,
					type=str,
					help='Input a constraint file')
parser.add_argument('-g',
					'--generate_constraint',
					nargs='+',
					type=str,
					help='''generate constraint file. Two or three arguments are taken, 
							first = PDB name, second = cst file name, third = weight,
							if the third argument is not inputted, it will assume the weight to be 0.25''')
parser.add_argument("-i", 
					'--input',
					nargs=1,
					type=str,
					help='input PDB file path/name')
parser.add_argument('-m',
					'--multi', 
					action='store_true',
					help='Relax multiple structures by multithreading')
parser.add_argument('-r',
					'--record',
					nargs=1,
					type=str,
					help='Record a log file or not')
parser.add_argument('--bb_off',
					action='store_false',
					help='''Turn backbone movemap off. No argument, the default is "on"
							when this option is called on the commandline, backbone movemap is turned off''')
parser.add_argument('--chi_off',
					action='store_false',
					help='''Turn sidechain movemap off. No argument, the default is "on"
							when this option is called on the commandline, sidechain movemap is turned off''')
parser.add_argument('--jump_off',
					action='store_false',
					help='''Turn jump movemap off. No argument, the default is "on"
							when this option is called on the commandline, jump movemap is turned off''')						


args = parser.parse_args()

#---------------------------------------------Section ends----------------------------------------------------


def main():
	#----------------------------Start timer--------------------------------
	start = time.time()
	#----------------------------Main function------------------------------
	bb = args.bb_off
	chi = args.chi_off
	jump = args.jump_off

	print(bb, chi, jump)
	print(type(bb), type(chi), type(jump))

	if args.generate_constraint:
		if len(args.generate_constraint) == 3:
			gen_cst(args.generate_constraint[0], args.generate_constraint[1], float(args.generate_constraint[2]))
		else:
			gen_cst(args.generate_constraint[0], args.generate_constraint[1])

	if args.input:
		if os.path.isfile(args.input[0]):
			if args.constraint:
				score = RunProcess(args.input[0], bb, chi, jump, cst=args.constraint[0])
			else:
				score = RunProcess(args.input[0], bb, chi, jump)

			if args.record:
				with open(args.record[0], 'a') as fin:
					fin.write(score)
			
			end = time.time()
			# Calculate wall time
			cal_wall_time(start, end)

		elif os.path.isdir(args.input[0]):

			PDBls = glob.glob(args.input[0]+'*.pdb')

			pool = mp.Pool(int(mp.cpu_count()*2/3))

			mmplist = [bb, chi, jump]
			if args.constraint:
				cstls = []
				for _ in range(len(PDBls)):
					cstls.append(args.constraint[0]) 
				combine_ls = []
				for i in PDBls:
					combine_ls.append([i]+mmplist)
				final_ls = []
				for i in range(len(combine_ls)):
					final_ls.append(combine_ls[i]+[cstls[i]])
				scores = pool.starmap(RunProcess, final_ls)
			else:
				final_ls = []
				for i in PDBls:
					final_ls.append([i]+mmplist)
				scores = pool.starmap(RunProcess, final_ls)

			if args.record:
				with open(args.record[0], 'a') as fin:
					for i in scores:
						fin.write(f'{i}\n')
		

		end = time.time()
		# Calculate wall time
		cal_wall_time(start, end)



if __name__ == '__main__':
	main()
