# This file was generated using the following command:
# ./configure 

# Rules that enable valgrind debugging ("make valgrind")

valgrind: .valgrind

.valgrind:
	echo -n > valgrind.out
	for x in $(TESTFILES); do echo $$x>>valgrind.out; valgrind ./$$x >/dev/null 2>> valgrind.out; done
	! ( grep 'ERROR SUMMARY' valgrind.out | grep -v '0 errors' )
	! ( grep 'definitely lost' valgrind.out | grep -v -w 0 )
	rm valgrind.out
	touch .valgrind


ALTAS = $(KALDI_ROOT)/tools/ATLAS
CONFIGURE_VERSION := 2
FSTROOT = $(KALDI_ROOT)/tools/openfst
OPENFST_VER = 1.3.4
OPENFST_GE_10400 = 0
OPENFSTLIBS = -L$(KALDI_ROOT)/tools/openfst/lib -lfst
OPENFSTLDFLAGS = -Wl,-rpath=$(KALDI_ROOT)/tools/openfst/lib
ATLASINC = $(ALTAS)/include
ATLASLIBS = $(ALTAS)/lib/libatlas.so.3 $(ALTAS)/lib/libf77blas.so.3 $(ALTAS)/lib/libcblas.so.3 $(ALTAS)/lib/liblapack_atlas.so.3 \
	    $(ALTAS)/lib/libgfortran.so.3 
# You have to make sure ATLASLIBS is set...

ifndef FSTROOT
$(error FSTROOT not defined.)
endif

ifndef ATLASINC
$(error ATLASINC not defined.)
endif

ifndef ATLASLIBS
$(error ATLASLIBS not defined.)
endif


CXXFLAGS = -finstrument-functions  -I.. \
      -DKALDI_DOUBLEPRECISION=0 -DHAVE_POSIX_MEMALIGN \
      -Wno-sign-compare -Wno-unused-local-typedefs -Winit-self \
      -DHAVE_EXECINFO_H=1 -DHAVE_CXXABI_H \
      -DHAVE_ATLAS -I$(ATLASINC) \
      -I$(FSTROOT)/include \
      $(EXTRA_CXXFLAGS) \
      -g # -O0 -DKALDI_PARANOID 

ifeq ($(KALDI_FLAVOR), dynamic)
CXXFLAGS += -fPIC
endif

LDFLAGS = $(OPENFSTLDFLAGS)
LDLIBS = $(EXTRA_LDLIBS) $(OPENFSTLIBS) $(ATLASLIBS) -lm -ldl
CC = arm-linux-gnueabihf-g++
CXX = $(CC)
AR = arm-linux-gnueabihf-ar
AS = arm-linux-gnueabihf-as
RANLIB = ranlib

CXXFLAGS += -DHAVE_SPEEX -I$(KALDI_ROOT)/src/../tools/speex/include 
LDLIBS += -L$(KALDI_ROOT)/src/../tools/speex/lib -lspeex
LDFLAGS += -Wl,-rpath=$(KALDI_ROOT)/src/../tools/speex/lib
