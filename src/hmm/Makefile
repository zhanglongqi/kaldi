all:


include ../kaldi.mk

TESTFILES = hmm-topology-test hmm-utils-test transition-model-test posterior-test

OBJFILES = hmm-topology.o transition-model.o hmm-utils.o posterior.o

LIBNAME = kaldi-hmm
ADDLIBS = ../tree/kaldi-tree.a ../matrix/kaldi-matrix.a ../util/kaldi-util.a \
          ../base/kaldi-base.a 

include ../makefiles/default_rules.mk

