#----------------------------Start timer--------------------------------
import time
from sys import argv
start = time.time()

#-------------------------Set up Pyrosetta------------------------------
from pyrosetta.rosetta import *
from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.rosetta.core.scoring import fa_rep, atom_pair_constraint
from pyrosetta.rosetta.core.pack.task.operation import ReadResfile, RestrictToRepacking
from pyrosetta.rosetta.protocols.relax import FastRelax
import itertools
init('-beta_nov16 True')

PDB = argv[1]


for i in itertools.count(start=1):
    if not os.path.exists(PDB[:-4]+'_R'+"{:02d}".format(i)+'.pdb'):
        if i > 1:
            file_incre = i-1
        else:
            file_incre = i
        break

outputFN = PDB[:-4]+'_R'+"{:02d}".format(file_incre)+'.pdb'


No_list = ['n','no']


pose = pose_from_pdb(PDB)
scorefxn = create_score_function('beta_nov16')

# Check if constraint file is needed
cst_check = input('\nDo you wish to use a contraint file?(y/n):')
cst_check = cst_check.lower()

if cst_check in No_list:
    print('No constraint file applied!!!!')
else:
    constraint_file = input('Please enter a path for the constraint file:')
    constraint = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    constraint.constraint_file(constraint_file)
    constraint.add_constraints(True)
    constraint.apply(pose)
    scorefxn.set_weight(atom_pair_constraint, 1)


score_before = scorefxn(pose)

movemap = MoveMap()
movemap.set_bb(True) 
movemap.set_chi(True)

relax = FastRelax()
relax.set_scorefxn(scorefxn)
relax.set_movemap(movemap)
relax.apply(pose)


#-----------------------------Out put-----------------------------------
score_after = scorefxn(pose)

print ('\n'*2, 'The initial score of this structure is:', '\n', score_before)
print ('The score after relaxing is:', '\n', score_after)
print ('\n','The distribution of scores', '\n')
scorefxn.show(pose)

pose.dump_pdb(outputFN)

#----------------------------Wall time----------------------------------
end = time.time()
d = end-start
if d >= 3600:
	print ('\n'*2, 'The time spent on relaxing with constraints was:', d/3600, 'hours')
if 60 <= d < 3600:
	print ('\n'*2, 'The time spent on relaxing with constraints was:', d/60, 'minutes')
if d < 60:
	print ('\n'*2, 'The time spent on relaxing with constraints was:', d, 'seconds')
