import subprocess
import os

#-------------------------------------------------------------------------------------------------

input_path = './data/regular_data'
dat_path = './dat_files/'
param_array = ["0","1"]


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

#for each file calls cpp exec
for r, d, f in os.walk(input_path):
	for file in f:
		file_path = os.path.join(r, file)

		# for o in (output_file_paths_mem_time ) :
		# 	(file_name, info) = o
		# 	curr_file = open(file_name,"a")
		# 	curr_file.write(file_path + "\t")
		# 	curr_file.close()

		for p in param_array:
			cp = subprocess.run(["./build/main", file_path, p], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			output = cp.stdout.decode('utf-8')
			for output_line in output.split("\n"):
				if "DAT" in output_line :
					begin_fp = output_line.find('_file_') + len("_file_")
					end_fp = output_line.rfind('_xy_values_')
					filepath = output_line[begin_fp:end_fp] #change name

					begin_v = output_line.rfind('_xy_values_') + len("_xy_values_")
					value = output_line[begin_v:]

					curr_file = open(filepath,"a")
					curr_file.write(file_path + "\t" + value + "\n")
					curr_file.close()

		# for o in (output_file_paths_mem_time ) :
		# 	(file_name, info) = o
		# 	curr_file = open(file_name,"a")
		# 	curr_file.write("\n") 
		# 	curr_file.close()