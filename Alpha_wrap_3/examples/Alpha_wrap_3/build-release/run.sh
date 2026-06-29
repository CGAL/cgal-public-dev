for mesh in ~/projects/GSoC-CGAL/dev/models/convavities/*off; do
    for ((alpha=50; alpha<=300; alpha+=50)); do
		for ((offset=50; offset<=300; offset+=50)); do
			./triangle_mesh_wrap $mesh $alpha $offset
			output=$(echo $mesh | awk -F"/" '{ print $NF }')
			output+=_
			output+=$alpha
			output+=_
			output+=$offset
			./stats.sh $output 1.2
		done        
    done
done
