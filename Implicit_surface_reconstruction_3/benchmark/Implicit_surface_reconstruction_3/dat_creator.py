import subprocess
import os

#-------------------------------------------------------------------------------------------------

input_path = './data/regular_data'
dat_path = './dat_files/'
param_array = ["0","1"]
artfct_array = ["position_noise"]

output_file_paths_mem_time = [ ("dat_files/0_mem_peak.dat","memory peak"),\
								("dat_files/1_mem_peak.dat","memory peak"),\
								("dat_files/0_time.dat","time"),\
								("dat_files/1_time.dat", "time") ]
# output_file_paths_position_noise = [ ("dat_files/position_noise_mean_dist_ptm.dat","mean distance ptm"),\
# 									("dat_files/position_noise_mean_dist_mtp.dat","mean distance mtp"),\
# 									("dat_files/position_noise_hausdorff_dist_ptm.dat","hausdorff distance ptm"),\
# 									("dat_files/position_noise_hausdorff_dist_mtp.dat","hausdorff distance mtp") ]

#creates dat files
for o in output_file_paths_mem_time :
	(file_name, info) = o
	curr_file = open(file_name,"w+")
	curr_file.write("# filename\tnb of pts\t" + info + "\n")
	curr_file.close()

# for o in output_file_paths_position_noise :
# 	(file_name, info) = o
# 	curr_file = open(file_name,"w+")
# 	curr_file.write("# filename\tnoise lvl\t")
# 	for p in param :
# 		curr_file.write(info + " for p" + p + "\t")
# 	curr_file.write("\n")
# 	curr_file.close()

for curr_param in param_array :
	curr_file = open(("./dat_files/" + curr_param + "_pos_noise_mean_dist_ptm.dat" ),"w+")
	curr_file.write("# filename\tlvl\t" + "mean_dist_ptm" + "\n")
	curr_file.close()

#for each file
for r, d, f in os.walk(input_path) :
	for file in f :
		input_file_path = os.path.join(r, file)
		input_xyz_file_path = input_file_path[:input_file_path.rfind(".")] + ".xyz"
		cp = subprocess.run(["./build/mesh_to_pwnl", input_file_path, input_xyz_file_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		#for each param
		for p in param_array :

			#calls reconstruction.cpp
			cp = subprocess.run(["./build/reconstruction", input_xyz_file_path, p], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

			output = cp.stdout.decode('utf-8')
			for output_line in output.split("\n") :
				#get mem + time + nb of pts
				if ("DAT" in output_line and "_nb_points_" in output_line) :
					begin = output_line.find('_nb_points_') + len('_nb_points_')
					end = output_line.find('_time_')
					nb_points = output_line[begin:end]
					
					begin = end + len('_time_')
					end = output_line.rfind('_mem_peak_')
					time = output_line[begin:end]

					begin = end + len('_mem_peak_')
					mem_peak = output_line[begin:]

					#stores mem + time in dat file
					curr_file = open(dat_path + p + "_time.dat", "a")
					curr_file.write(input_file_path + "\t" + nb_points + "\t" + time + "\n")
					curr_file.close()

					curr_file = open(dat_path + p + "_mem_peak.dat", "a")
					curr_file.write(input_file_path + "\t" + nb_points + "\t" + mem_peak + "\n")
					curr_file.close()
				

		# for p in param_array:
		# 	cp = subprocess.run(["./build/main", input_file_path, p], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		# 	output = cp.stdout.decode('utf-8')
		# 	for output_line in output.split("\n"):
		# 		if "DAT" in output_line :
		# 			begin_fp = output_line.find('_file_') + len("_file_")
		# 			end_fp = output_line.rfind('_xy_values_')
		# 			filepath = output_line[begin_fp:end_fp] #change name

		# 			begin_v = output_line.rfind('_xy_values_') + len("_xy_values_")
		# 			value = output_line[begin_v:]

		# 			curr_file = open(filepath,"a")
		# 			curr_file.write(input_file_path + "\t" + value + "\n")
		# 			curr_file.close()

		# for o in (output_file_paths_mem_time ) :
		# 	(file_name, info) = o
		# 	curr_file = open(file_name,"a")
		# 	curr_file.write("\n") 
		# 	curr_file.close()