

all:

include ../kaldi.mk

TESTFILES = feature-mfcc-test feature-plp-test feature-fbank-test \
         feature-functions-test pitch-functions-test feature-sdc-test \
         resample-test online-feature-test sinusoid-detection-test \
         signal-test

OBJFILES = feature-functions.o  \
           mel-computations.o  \
          

LIBNAME = kaldi-feat

ADDLIBS = ../transform/kaldi-transform.a ../gmm/kaldi-gmm.a ../tree/kaldi-tree.a \
	../util/kaldi-util.a ../matrix/kaldi-matrix.a ../base/kaldi-base.a ../thread/kaldi-thread.a

include ../makefiles/default_rules.mk

