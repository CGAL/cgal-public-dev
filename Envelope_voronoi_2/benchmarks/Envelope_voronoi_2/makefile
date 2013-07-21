BASEDIR =.

USE_QT=1
USE_CORE=1

include $(ROOT)/include/make/cgaldef.mak

FILTERED_ZONE = 1

BASENAME = bench
INSTALLDIR0 = $(BINDIR)
CPPSOURCES = envelope_voronoi.cpp
CPPSOURCES+= Option_parser.cpp
CPPSOURCES+= Bench_apollonius.cpp
CPPSOURCES+= Bench_envelope_apollonius.cpp
CPPSOURCES+= Bench_envelope_voronoi.cpp
CPPSOURCES+= Bench_triangulation_voronoi.cpp
CPPSOURCES+= Bench_sphere.cpp

TARGET0 = $(BASENAME)

LOBJDIR =

ifeq ($(DEBUG),1)
GCPPOPTS = -g
# LCPPDEFS+= -DCGAL_DEBUG_SWEEP_LINE
# LCPPDEFS+= -DCGAL_PROFILE
# GCPPOPTS+= -DVERBOSE
else
GCPPOPTS = -O3

ifeq ($(BENCH_TRAITS), $(CORE_CONIC_TRAITS))
# GCPPOPTS = -O2 -fno-strict-aliasing
GCPPOPTS = -O3
endif

ifeq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
ifeq ($(BENCH_NT), $(NIX_CORE_FIELD_WITH_SQRT_NT))
# GCPPOPTS = -O2 -fno-strict-aliasing
GCPPOPTS = -O3
endif
endif

ifeq ($(BENCH_TRAITS), $(CK_CIRCLE_TRAITS))
ifeq ($(BENCH_NT), $(CORE_EXPR_NT))
# GCPPOPTS = -O2 -fno-strict-aliasing
GCPPOPTS = -O3
endif
endif

endif

# Planar map:
TARGET0 := $(TARGET0)Voronoi
LOBJDIR :=$(LOBJDIR)_voronoi

ifeq ($(BENCH_FONE), $(FILTERED_ZONE))
TARGET0 := $(TARGET0)Fone
LOBJDIR :=$(LOBJDIR)_fone
LCPPDEFS+= -DFILTERED_ZONE -DCONVEX_ZONE
endif

# Window system
ifeq ($(USE_CGAL_WINDOW), 1)
LCPPDEFS+= -DUSE_CGAL_WINDOW
TARGET0 := $(TARGET0)CgalWindow
LOBJDIR :=$(LOBJDIR)_cgal_window
endif

# Put is all together:
TARGET0 := $(TARGET0)$(EXEFILESUFFIX)

LCPPINCS = -I./
LCPPINCS+= -I$(BASEDIR)/../../include
LCPPINCS+= -I$(BASEDIR)
LCPPINCS+= -I$(CGAL_WORKDIR)/trunk/Benchmark/include
LCPPINCS+= -I$(CGAL_WORKDIR)/trunk/Arrangement_on_surface_2/include
LCPPINCS+= -I$(CGAL_WORKDIR)/trunk/Envelope_3/include

# LCPPINCS+= -I$(HOME)/trees/old/tcgl/trunk/src/apps/Voronoi/include
# LCPPINCS+= -I$(CGAL_WORKDIR)/trunk/Benchmark/include
# LCPPINCS+= -I$(CGAL_WORKDIR)/trunk/Arrangement_2/include
# LCPPINCS+= -I$(CGAL_WORKDIR)/trunk/Envelope_3/include

# ifneq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
LCPPINCS+= -I$(COREROOT)/inc
# endif

LCPPINCS+= $(CGALCPPINCS)

LLDLIBS = -lCGALBenchmark
LLDLIBS+= -lCGALcore++
ifeq ($(USE_QT), 1)
LLDLIBS+= $(CGALQTLIB) $(QTLIB)
endif
LLDLIBS+= $(CGALLIB)
ifeq ($(USE_LEDA), 1)
LLDLIBS+= $(LEDALIBS)
endif
LLDLIBS+= $(GMPLIBS)
LLDLIBS+= -lm
LLDLIBS+= -lboost_program_options
LLDLIBS+= -lboost_filesystem
LLDLIBS+= -lboost_thread

ifeq ($(USE_EXACUS), 1)
LCPPDEFS+= -DUSE_EXACUS

CPPSOURCES+= Bench_EXACUS_apollonius.cpp
CPPSOURCES+= Bench_EXACUS_Qdx_apollonius.cpp

LCPPINCS+= -I$(EXACUS_DIR)/include

LLDLIBS+= -lNiX
LIBDIR+= -L$(EXACUS_DIR)/lib
endif

include $(ROOT)/include/make/cgalrul.mak

# Dependencies:
$(BASENAME).moc: $(BASENAME).C
	${QT_MOC} -o $@ $<

$(BASENAME).o: $(BASENAME).moc

include $(BASEDIR)/voronoi.mak
