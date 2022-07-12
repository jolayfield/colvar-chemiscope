import numpy as np 
import ase
import ase.io
from chemiscope import write_input
import plumed
import sys
import argparse
import warnings


warnings.filterwarnings("ignore")

# Create the parser
my_parser = argparse.ArgumentParser(description='Generate a chemiscope json.gz file from a trajectory')

# Add the arguments
my_parser.add_argument('Trajectory',
                       metavar='traj',
                       type=str,
                       help='the pdb trajectory file.')
my_parser.add_argument('COLVAR', metavar='colvar file', type=str, help='The colvar file to analyze.')
my_parser.add_argument('Output', metavar='output_file', type=str, help='The full name of the chemiscope file (it should end in `.json.gz`)')
my_parser.add_argument('-p','--plumed', required=False)

# Execute the parse_args() method
args = my_parser.parse_args()

traj_file = args.Trajectory
colvar_file=args.COLVAR
out_file = args.Output

with open(colvar_file) as f:
    colvar_lines = f.readlines()

print('colvar')
CVS = colvar_lines[0].split()[3:]
print('CVS',CVS)

np_colvar = np.loadtxt(colvar_file, skiprows=1)

# Read in trajectory using ase
traj = ase.io.read(traj_file,':')

# Setup plumed object to do calculation
p = plumed.Plumed()
p.cmd("setMDEngine","python")
p.cmd("setTimestep", 1.)
p.cmd("setKbT", 1.)
natoms = len(traj[0].positions)
p.cmd("setNatoms",natoms)
p.cmd("setLogFile","test.log")
p.cmd("init")

# Read plumed input 
if args.plumed:
    p.cmd("readInputLine",f'INCLUDE FILE={args.plumed}')
else:
    p.cmd("readInputLine",f"MOLINFO STRUCTURE={traj_file}") 

# # Loop over trajectory and get data from plumed
nfram, tt, box = 0, [], np.array([[100.,0,0],[0,100.,0],[0,0,100]])
charges, forces, virial = np.zeros(natoms,dtype=np.float64), np.zeros([natoms,3]), np.zeros((3,3),dtype=np.float64)

for ts in traj :
    p.cmd("setStep",nfram)
    p.cmd("setBox",box )
    p.cmd("setMasses", ts.get_masses() )
    p.cmd("setCharges", charges )
    pos = np.array(ts.get_positions(), dtype=np.float64 )
    p.cmd("setPositions", pos )
    p.cmd("setForces", forces )
    p.cmd("setVirial", virial )
    p.cmd("calc")
    tt.append(nfram)
    nfram = nfram + 1
    
atom_dict= {'H':1, 'C':6, 'N':7, 'O':8, 'S':16, '1':1, '2':1,'3':1}


for frame in traj:
    frame.numbers = np.array(
        [
            atom_dict[am[0]] 
            for am in frame.arrays["atomtypes"]
        ]
    )

    
# This constructs the dictionary of properties for chemiscope
properties = {}  
properties["time"]={"target": "structure","values": tt,"description": "Simulation step number"}

for i,CV in enumerate(CVS):
    properties[CV]= {'target':'structure', "values":np_colvar[:,i+1], "description":'fix this later'} #CVs[i]}
# # This generates our chemiscope output
write_input(out_file, frames=traj, properties=properties )
