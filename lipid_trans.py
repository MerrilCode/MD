#!/usr/bin/python
import numpy as np
import re
import csv 
import sys,getopt
# Directory and input files for comparison
data_dir = '/Users/merrilmathew/Dropbox/md/dppc_opo80/'

solvated_file_name = 'solvated_ions.gro'
md_file_name = 'md.gro'
solvated_f = open(data_dir + solvated_file_name,'rb') 
md_f = open(data_dir + md_file_name,'rb') 
solvated_lines = solvated_f.readlines()[2:-1] 
md_lines = md_f.readlines()[2:-1]
assert len(solvated_lines)==len(md_lines)
# Declaring variables and methods to store user input values from command line
lower_z_atom_plane = '' 
upper_z_atom_plane = '' 
outfile = ''
dimension = '' 
ref_residue = ''
ref_atom = ''

# Declaring lists to store specific lines of the residue and atom from the two files to compare. Declaring lists to store the coordinates of the XYZ dimension of interest from the files.
solvated_list = []
md_list = []
solvated_cord = [] 
md_cord = []

# Command line user input argument fetching and storing
def main(argv): 
	try:
		opts,args = getopt.getopt(argv,"hr:a:l:u:d:o:", ["help","res","atom","lower_plane","upper_plane","dimension","ofile="])
		except getopt.GetoptError: 
			print 'test.py -o <outputfile>' sys.exit(2)
		for o, a in opts:

			if o == '-h':
				print 'test.py -o <outputfile>'
			elif o in ("-r","--res"): 
				global ref_residue 
				ref_residue = a
			elif o in ("-a","--atom"): 
				global ref_atom 
				ref_atom = a
			elif o in("-l","--lower_plane"): 
				global lower_z_atom_plane 
				lower_z_atom_plane = float(a)
			elif o in("-u","--upper_plane"): 
				global upper_z_atom_plane
			elif o in("-d","--dimension"): 
				global dimension
				dimension = a
			elif o in ("-o", "--ofile"): 
				global outfile
				outfile = a

if __name__ == "__main__": 
	main(sys.argv[3:])
# Regular expression to match all the lines contains the reference residue and atom
def reg_search(str):
	my_regex = r'^.*' + re.escape(ref_residue) + r'.*\b' + re.escape(ref_atom) + r'\b.*$'
	match = re.search(my_regex,str,re.IGNORECASE) 
	return match;

for atom_idx_sol, atom_idx_md in zip(solvated_lines,md_lines):
	matchObj_sol = reg_search(atom_idx_sol) matchObj_md = reg_search(atom_idx_md) 

	if matchObj_sol:
		solvated_list.append(atom_idx_sol) 

	if matchObj_md:
		md_list.append(atom_idx_md)

# Store all the xyz coordinates of solvated_list into solvated_cord
for atom_idx in range(len(solvated_list)):
	solvated_line = solvated_list[atom_idx].split() 
	solvated_cord.append(map(lambda x:float(x),solvated_line[3:6]))
# Store all the xyz coordinates of md_list into md_cord
for atom in range(len(md_list)):
	md_line = md_list[atom].split() 
	md_cord.append(map(lambda x: float(x),md_line[3:6]))

# If the user input for the dimension is X then compute the difference between the coordinates of the same atom and store into x_diffs.
x_diffs =[]

if dimension == 'X':
for atom_idx in range(len(md_cord)):
	x_diffs.append(md_cord[atom_idx][0]-solvated_cord[atom_idx][0]) 
f = open(str(outfile),'a')
wtr = csv.writer(f,delimiter= '\t')
for atom_idx in range(len(x_diffs)):
	atom_num = atom_idx
	dist = x_diffs[atom_idx] 
	wtr.writerow((atom_num, "%.3f" % dist))

# If the user input for the dimension is Y then compute the difference between the coordinates of the same atom and store into x_diffs
elif dimension == 'Y':
for atom_idx in range(len(md_cord)): 
	x_diffs.append(md_cord[atom_idx][1]-solvated_cord[atom_idx][1])
f = open(str(outfile),'a')
wtr = csv.writer(f,delimiter= '\t') 
for atom_idx in range(len(x_diffs)):
	atom_num = atom_idx
	dist = x_diffs[atom_idx] 
	wtr.writerow((atom_num, "%.3f" % dist))
# If the user input for dimension (-d) is Z then split the coordinates into lower monolayer/plane (-l) and upper monolayer/plane(-u) and store into lists.
elif dimension == 'Z':
	lower_idxs = []
	for atom_idx in range(len(solvated_cord)):
	if np.abs(solvated_cord[atom_idx][2]-lower_z_atom_plane)<1: 
		lower_idxs.append(atom_idx)
# Checking if the z-coordinates of the atoms is within a range of the upper plane, i.e the second plane of the slab and storing into list.
upper_idxs = []
for atom_idx in range(len(solvated_cord)):
	if np.abs(solvated_cord[atom_idx][2]-upper_z_atom_plane)<1: 
		upper_idxs.append(atom_idx)
# Checking if all atoms are stored into upper and lower plane lists from the previous two steps
assert len(lower_idxs) + len(upper_idxs) == len(solvated_cord)
lower_solvs = [solvated_cord[idx] for idx in lower_idxs] 
upper_solvs = [solvated_cord[idx] for idx in upper_idxs] 
lower_mds = [md_cord[idx] for idx in lower_idxs] 
upper_mds = [md_cord[idx] for idx in upper_idxs]
lower_dists = [] 
lower_trans = []

# Computing and storing the distance between the atom coordinates from the two files in the lower monolayer/plane and storing the sign of change to detect the direction of movement.
for atom_idx in range(len(lower_solvs)):
	z_solv = lower_solvs[atom_idx][2]
	z_md = lower_mds[atom_idx][2]
	lower_dists.append(np.abs(z_solv-z_md)) 
	lower_trans.append([(z_solv-lower_z_atom_plane),(z_md-lower_z_atom_plane)])

# Printing the average distance moved by the reference residue and atom and the direction of the movement to standard output.
print"\nAverage movement in lower monolayer: ", sum(lower_dists)/len(lower_dists)
x = np.asarray(lower_trans) initial_pos = sum(x[:,0]) 
final_pos = sum(x[:,1])
if initial_pos < 0:
	print " Atom started from Air/Water interface" 
elif initial_pos > 0:
	print " Atom started from Air/Water interface" 
if final_pos < 0:
	print "Atom moved towards Air" 
elif final_pos > 0:
	print "Atom moved towards Water"

# Computing and storing the distance between the atom coordinates from the two files in the upper monolayer/plane and storing the sign of change to detect the direction of movement.
upper_dists = []
upper_trans = []
for atom_idx in range(len(upper_solvs)):
z_solv = upper_solvs[atom_idx][2]
z_md = upper_mds[atom_idx][2]
upper_dists.append(np.abs(z_solv-z_md)) 
upper_trans.append([(z_solv-upper_z_atom_plane),(z_md-upper_z_atom_plane)])

print"\nAverage movement in upper monolayer: ", sum(upper_dists)/len(upper_dists) 
y = np.asarray(upper_trans)
initial_pos1 = sum(y[:,0])
final_pos1 = sum(y[:,1])
if initial_pos1 < 0:
	print " Atom started from Air/Water interface"
elif initial_pos1 > 0:
	print " Atom started from Air/Water interface"
if final_pos1 < 0:
	print "Atom moved towards Water"
elif final_pos1 > 0:
	print "Atom moved towards the Air"

# writing the atom number, distance moved, intial position and final position to a CSV file.
f = open(str(outfile),'a')
wtr = csv.writer(f,delimiter= '\t')
for atom_idx in range(len(lower_solvs)):
	atom_num = atom_idx
	dist_lower = lower_dists[atom_idx]
	initial = x[atom_idx][0]
	final = x[atom_idx][1]
	wtr.writerow((atom_num, "%.3f" % dist_lower, "%.3f" % initial, "%.3f" % final))

for atom_idx in range(len(upper_solvs)):
	dist_upper = upper_dists[atom_idx]
	initial1 = y[atom_idx][0]
	final1 = y[atom_idx][1]
	wtr.writerow((atom_idx, "%.3f" % dist_upper,"%.3f" % initial1, "%.3f" % final1))