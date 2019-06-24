for filename in $1/*.off; do
	../build/surface_mesh_statistics $filename;
done
