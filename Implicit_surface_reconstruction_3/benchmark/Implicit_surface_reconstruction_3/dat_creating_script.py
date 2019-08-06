import subprocess
import os

#-------------------------------------------------------------------------------------------------

path = './data'
param = ["0","1"]

output_file_paths = [ ("dat_files/memory_peak.dat","memory peak"),("dat_files/time.dat", "time") ]

#creates dat files
for o in output_file_paths :
	(file_name, info) = o
	curr_file = open(file_name,"w+")
	curr_file.write("# filename\tnb of pts\t")
	for p in param :
		curr_file.write(info + " for p" + p + "\t")
	curr_file.write("\n")
	curr_file.close()


#for each file calls cpp exec
for r, d, f in os.walk(path):
	for file in f:
		file_path = os.path.join(r, file)

		for o in output_file_paths :
			(file_name, info) = o
			curr_file = open(file_name,"a")
			curr_file.write(file_path + "\t")
			curr_file.close()

		for p in param:
			cp = subprocess.run(["./build/main", file_path, p])

		for o in output_file_paths :
			(file_name, info) = o
			curr_file = open(file_name,"a")
			curr_file.write("\n") 
			curr_file.close()