import subprocess
import os

#-------------------------------------------------------------------------------------------------

input_path = './data/regular_data'
dat_path = './dat_files/'
ext = '.dat'
param_array = ["0","1"]
artefact_array = [("_pos_noise", "./build/position_noise_generator", "Position noise")]
lvl_array = ["1" , "2"]
mt_array = [("_mean_dist_ptm","Mean distance from pts to mesh")]
time_mem_array = [("_mem_peak","Memory peak (bytes)"),	("_time", "Time (s)")]


#creates dat files
for p in param_array :
	for tm in time_mem_array :
		(file_name, info) = tm
		curr_file = open(dat_path + p + file_name + ext,"w+")
		curr_file.write("# filename\t\t\t\t\t\tnb of pts\t" + info + "\n")
		curr_file.close()
	for a in artefact_array :
		(art_file_name, art_exec, art_info) = a
		for mt in mt_array :
			(mt_file_name,info) = mt
			curr_file = open(dat_path + p + art_file_name + mt_file_name + ext,"w+")
			curr_file.write("# filename\t\t\t\t\t\tlvl\t" + info + "\n")
			curr_file.close()


#fills dat files
for r, d, f in os.walk(input_path) :
	for file in f :

		input_file_path = os.path.join(r, file)
		input_xyz_file_path = input_file_path[:input_file_path.rfind(".")] + ".xyz"
		cp = subprocess.run(["./build/mesh_to_pwnl", input_file_path, input_xyz_file_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		#for each param
		for p in param_array :

			lvl0_values = []
			#calls reconstruction.cpp (no artefact)
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
				#get lvl 0
				elif ("DAT" in output_line and "_measure_type_" in output_line) :
					begin = output_line.find('_measure_type_') + len('_measure_type_')
					end = output_line.rfind('_value_')
					m_type = output_line[begin:end]

					begin = end + len('_value_')
					value = output_line[begin:]
					lvl0_values.append((m_type,value))

			for a in artefact_array :

				(art_file_name, art_exec, art_info) = a
				#stores value into file
				for c in lvl0_values :
					(m_type, value) = c
					curr_file = open(dat_path + p + art_file_name + m_type, "a")
					curr_file.write(input_file_path + "\t" + "0" + "\t"+ value + "\n")
					curr_file.close()

				for lvl in lvl_array :

					output_xyz_file = "output.xyz"
					# output_xyz_file = file[:file.rfind(".")]+p+art_file_name+lvl+".xyz"
					cp = subprocess.run([art_exec, input_xyz_file_path, output_xyz_file, lvl],\
										 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					cp = subprocess.run(["./build/reconstruction", output_xyz_file, p],\
										 stdout=subprocess.PIPE, stderr=subprocess.PIPE)

					output = cp.stdout.decode('utf-8')
					for output_line in output.split("\n") :
						#get value
						if ("DAT" in output_line and "_measure_type_" in output_line) :
							begin = output_line.find('_measure_type_') + len('_measure_type_')
							end = output_line.rfind('_value_')
							m_type = output_line[begin:end]

							begin = end + len('_value_')
							value = output_line[begin:]

							#stores value into file
							curr_file = open(dat_path + p + art_file_name + m_type, "a")
							curr_file.write(input_file_path + "\t" + lvl + "\t"+ value + "\n")
							curr_file.close()

#creates .ps gnuplot script
charts_path = './charts/'
ps_file_name = 'plot_script.ps'
ps_file = open(ps_file_name,"w+")
ps_file.write("set terminal pdf\n")

for tm in time_mem_array :
	(tm_file_name, tm_info) = tm
	chart_file_name = "'" + charts_path + "chart" + tm_file_name + ".pdf" + "'"
	# chart_title = 
	chart_xlabel = '"Number of points"'
	chart_ylabel = '"' + tm_info + '"'
	ps_file.write("set output " + chart_file_name + "\n")
	# ps_file.write("set title " + chart_title + "\n") 
	ps_file.write("set xlabel " + chart_xlabel + "\n")
	ps_file.write("set ylabel " + chart_ylabel + "\n")
	ps_file.write("set key" + "\n")
	plot_line = "plot "
	for p in param_array :
		plot_line = plot_line + "'" + dat_path + p + tm_file_name + ".dat' "
		plot_line = plot_line + "using 2:3 smooth unique with linespoints "
		plot_line = plot_line + "title 'param " + p + "', "
	ps_file.write(plot_line + "\n")

for a in artefact_array :
	(art_file_name, art_exec, art_info) = a
	for mt in mt_array :
		(mt_file_name, mt_info) = mt
		chart_file_name = "'" + charts_path + "chart" + art_file_name + mt_file_name + ".pdf" + "'"
		# chart_title = 
		chart_xlabel = '"Level of ' + art_info + '"'
		chart_ylabel = '"' + mt_info + '"'
		ps_file.write("set output " + chart_file_name + "\n")
		# ps_file.write("set title " + chart_title + "\n") 
		ps_file.write("set xlabel " + chart_xlabel + "\n")
		ps_file.write("set ylabel " + chart_ylabel + "\n")
		ps_file.write("set key" + "\n")
		plot_line = "plot "
		for p in param_array :
			plot_line = plot_line + "'" + dat_path + p + art_file_name + mt_file_name + ".dat' "
			plot_line = plot_line + "using 2:3 smooth unique with linespoints "
			plot_line = plot_line + "title 'param " + p + "', "
		ps_file.write(plot_line + "\n")
	ps_file.write("\n")

ps_file.close()


#creates bash commands
cmds = open("cmds.sh", "w+")
cmds.write('gnuplot -e "load ' + "'"  + ps_file_name + "'" + '"' + "\n") #gnuplot cmd = gnuplot -e "load '<ps_file_name>'"
cmds.write("rm output.xyz" + "\n")
cmds.write("for f in " + input_path + "/*.stl; do\n")
cmds.write('  [ -e "${f%.*}.xyz" ] && rm -- "${f%.*}.xyz"\ndone \n')
cmds.write("rm " + ps_file_name + "\n")
cmds.write("rm cmds.sh")
cmds.close()
os.system("chmod +x cmds.sh")
os.system('./cmds.sh')