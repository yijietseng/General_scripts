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

def set_cst(sf, pose):
    constraint_file = input('Please enter a path for the constraint file:').strip()
    constraint = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    constraint.constraint_file(constraint_file)
    constraint.add_constraints(True)
    constraint.apply(pose)
    sf.set_weight(atom_pair_constraint, 1)


def set_movemap(set_bb=True, set_chi=True, set_jump=False):
    movemap = MoveMap()
    movemap.set_bb(set_bb) 
    movemap.set_chi(set_chi)
    movemap.set_jump(set_jump)
    return movemap

def RunRelax(sf, mmp, pose, apply=True):
    relax = FastRelax()
    relax.set_scorefxn(sf)
    relax.set_movemap(mmp)
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


def RunProcess(PDB, cst=''):

    pose = pose_from_pdb(PDB)
    outputFN = format_outputFN(PDB)
    #scorefxn = create_score_function('beta_nov16')
    scorefxn = get_fa_scorefxn()

    if cst:
        set_cst(scorefxn, pose)
    # Set movemap
    mmp = set_movemap()
    
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
parser.add_argument("-i", 
                    '--input',
                    nargs=1,
                    type=str,
                    required=True,
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


args = parser.parse_args()

#---------------------------------------------Section ends----------------------------------------------------


def main():
    #----------------------------Start timer--------------------------------
    start = time.time()
    #----------------------------Main function------------------------------
    if os.path.isfile(args.input[0]):
        if args.constraint:
            score = RunProcess(args.input[0], cst=args.constrain[0])
        else:
            score = RunProcess(args.input[0])

        if args.record:
            with open(args.record[0], 'a') as fin:
                fin.write(score)
        
        end = time.time()
        # Calculate wall time
        cal_wall_time(start, end)

    elif os.path.isdir(args.input[0]):

        PDBls = glob.glob(args.input[0]+'*.pdb')

        pool = mp.Pool(int(mp.cpu_count()*2/3))
        if args.constraint:
            cstls = []
            for _ in range(len(PDBls)):
                cstls.append(args.constraint[0])
            combine_ls = list(zip(PDBls,cstls))
            scores = pool.map(RunProcess, combine_ls)
        else:
            scores = pool.map(RunProcess, PDBls)

        if args.record:
            with open(args.record[0], 'a') as fin:
                fin.write(scores)
        

        end = time.time()
        # Calculate wall time
        cal_wall_time(start, end)



if __name__ == '__main__':
    main()
