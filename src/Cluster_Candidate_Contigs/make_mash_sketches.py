import os

input_path = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash/'
Mash_Path = '/fs/cbcb-scratch/hsmurali/Mash-2.1.1/'
op_path = '/fs/cbcb-scratch/hsmurali/Potential_Phages_Updated/Mash_Sketch/'

try:
	os.mkdir(op_path)
except Exception:
	pass
	
files = os.listdir(input_path)
fasta_files = []
for f in files:
        if ".msh" not  in f:
                fasta_files.append(input_path+f)

batch_size = 100

counter = 1

os.chdir(Mash_Path)

for i in range(0, len(fasta_files), batch_size):
	print(counter)
	mash_command = './mash sketch -o '+op_path+'reference_'+str(counter)
	for j in range(0, batch_size):
		index = i+j
		if i+j < len(fasta_files):
			mash_command += ' ' + fasta_files[index]
	mash_command += ' -k 11 -s 100000'
	counter += 1
	result = os.popen(mash_command).read()
		
