all_voronoi: voronoi \
	voronoi_no_fone

voronoi:
	$(MAKEF) "BENCH_FONE=$(FILTERED_ZONE)"

voronoi_no_fone:
	$(MAKEF)

# Install:
voronoi_inst:
	$(MAKEF) "BENCH_FONE=$(FILTERED_ZONE)" install

voronoi_no_fone_inst:
	$(MAKEF) install

all_voronoi_inst: voronoi_inst \
	voronoi_no_fone_inst
