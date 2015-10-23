all:

EXTRA_CXXFLAGS = -Wno-sign-compare -O3
include ../kaldi.mk

TESTFILES = 

OBJFILES = lattice-simple-decoder.o lattice-faster-decoder.o \
    faster-decoder.o  decoder-wrappers.o

LIBNAME = kaldi-decoder

ADDLIBS = ../transform/kaldi-transform.a ../tree/kaldi-tree.a ../lat/kaldi-lat.a \
     ../gmm/kaldi-gmm.a ../hmm/kaldi-hmm.a ../util/kaldi-util.a \
     ../base/kaldi-base.a ../matrix/kaldi-matrix.a 

include ../makefiles/default_rules.mk

