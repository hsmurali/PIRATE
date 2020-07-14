import os

Mash_Sketch_Paths = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash_Sketch/'
Prog_Path = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/run_mash_sketch_distance.py'

files =os.listdir(Mash_Sketch_Paths)
Mash_Commands = []

for f in files:
	cmd = Prog_Path + ' '  + f + '\n'
	Mash_Commands.append(cmd)

f = open('/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash_Commands.txt','w')
f.writelines(Mash_Commands)
