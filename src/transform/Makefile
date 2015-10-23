all:

include ../kaldi.mk

TESTFILES = regtree-fmllr-diag-gmm-test lda-estimate-test \
      regression-tree-test fmllr-diag-gmm-test \
      regtree-mllr-diag-gmm-test fmpe-test fmllr-raw-test

OBJFILES = regression-tree.o regtree-mllr-diag-gmm.o  \
    regtree-fmllr-diag-gmm.o  transform-common.o fmllr-diag-gmm.o


LIBNAME = kaldi-transform

ADDLIBS = ../gmm/kaldi-gmm.a ../thread/kaldi-thread.a ../tree/kaldi-tree.a \
   ../util/kaldi-util.a ../matrix/kaldi-matrix.a ../base/kaldi-base.a

include ../makefiles/default_rules.mk
