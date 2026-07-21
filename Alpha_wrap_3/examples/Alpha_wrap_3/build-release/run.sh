# Cleanup
rm *off
rm steps.csv
rm -r results
mkdir results
root=$(pwd)

# Run
for mesh_path in ~/projects/GSoC-CGAL/dev/models/convavities/*off; do
    for ((alpha=20; alpha<=200; alpha+=20)); do
		for ((offset=20; offset<=200; offset+=20)); do
			mesh_name=$(echo $mesh_path | awk -F"/" '{ print $NF }')
			base_name=$(echo $mesh_name | awk -F"." '{ print $1 }')
			output_name="$base_name"_"$alpha"_"$offset"
			output_dir=$root/results/$output_name
			mkdir $output_dir
			cd $output_dir
			RELAXED=1 SPHERE_MARCHING_OLD=1 DICHOTOMY=1 SPHERE_MARCHING=2\
					  "$root"/triangle_mesh_wrap $mesh_path $alpha $offset
		done        
    done
done
