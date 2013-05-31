##
##  Default Makefile options from Configure script
##  (Begin of automatically generated section)
##

CXXFLAGS= -DNDEBUG -O2 -Wall -Wno-unused -pthread 
LIBS=
LDFLAGS= -lpthread -lrt -lstdc++
INCDIR=$(SRCDIR)/include
# CC=$(HOME)/thirdparty_libs/gcc-4.1/gcc-4.1.2/bin/gcc-4.1
OBJS=
LIB=
# INCLUDES=-I$HOME/thirdparty_libs/gcc-4.1/gcc-4.1.2/include/c++/4.1.2 -pthread

PROGRAMS = xtri.cgi

all: xtri.cgi

xtri_cgi_SOURCES = xtri_client.C cgic.C
xtri_cgi_OBJECTS = xtri_client.o cgic.o

.SUFFIXES: .C .S .c .lo .o .obj .s

.C.o:	
	$(CC) -c $(CXXFLAGS) $(INCLUDES) $<

xtri.cgi: $(xtri_cgi_OBJECTS)
	$(CC) $(xtri_cgi_OBJECTS) $(LDFLAGS) -o xtri.cgi

clean:
	rm -f $(LIB) $(OBJS) $(SHLIBS) $(OBJS_PIC)

distclean: clean
	rm -f Makefile



