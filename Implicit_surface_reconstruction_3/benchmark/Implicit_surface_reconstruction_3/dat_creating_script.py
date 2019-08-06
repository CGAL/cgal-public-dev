import subprocess
import os

#-------------------------------------------------------------------------------------------------

path = './data'
param = ["0","1"]

#creates memory peak dat file
mem_chart_file = open("dat_files/mem_chart.dat","w+")
mem_chart_file.write("# filename\tnb of pts\t")
for p in param :
	mem_chart_file.write("memory peak for p" + p + "\t") 
mem_chart_file.write("\n")
mem_chart_file.close();

#creates time dat file
time_chart_file = open("dat_files/time_chart.dat","w+")
time_chart_file.write("# filename\tnb of pts\t")
for p in param :
	time_chart_file.write("time taken for p" + p + "\t") 
time_chart_file.write("\n")
time_chart_file.close();

#for each file calls cpp exec
for r, d, f in os.walk(path):
	for file in f:
		file_path = os.path.join(r, file)
		mem_chart_file = open("dat_files/mem_chart.dat","a")
		mem_chart_file.write(file_path + "\t")
		mem_chart_file.close()
		time_chart_file = open("dat_files/time_chart.dat","a")
		time_chart_file.write(file_path + "\t")
		time_chart_file.close()
		for p in param:
			cp = subprocess.run(["./build/mem_chart", file_path, p])
		mem_chart_file = open("dat_files/mem_chart.dat","a")
		mem_chart_file.write("\n") 
		mem_chart_file.close()
		time_chart_file = open("dat_files/time_chart.dat","a")
		time_chart_file.write("\n") 
		time_chart_file.close()
