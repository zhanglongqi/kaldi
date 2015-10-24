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


CONFIGURE_VERSION := 2
FSTROOT = /home/longqi/kaldi-trunk/tools/openfst
OPENFST_VER = 1.3.4
OPENFST_GE_10400 = 0
OPENFSTLIBS = -L/home/longqi/kaldi-trunk/tools/openfst/lib -lfst
OPENFSTLDFLAGS = -Wl,-rpath=/home/longqi/kaldi-trunk/tools/openfst/lib
ATLASINC = /home/longqi/kaldi-trunk/tools/ATLAS/include
ATLASLIBS = /usr/lib/libatlas.so.3 /usr/lib/libf77blas.so.3 /usr/lib/libcblas.so.3 /usr/lib/liblapack_atlas.so.3
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


CXXFLAGS = -finstrument-functions  -Wall -I.. \
      -DKALDI_DOUBLEPRECISION=0 -DHAVE_POSIX_MEMALIGN \
      -Wno-sign-compare -Wno-unused-local-typedefs -Winit-self \
      -DHAVE_EXECINFO_H=1 -shared -DHAVE_CXXABI_H \
      -DHAVE_ATLAS -I$(ATLASINC) \
      -I$(FSTROOT)/include \
      $(EXTRA_CXXFLAGS) \
      -g # -O0 -DKALDI_PARANOID 

ifeq ($(KALDI_FLAVOR), dynamic)
CXXFLAGS += -fPIC
endif

LDFLAGS = -shared $(OPENFSTLDFLAGS)
LDLIBS = $(EXTRA_LDLIBS) $(OPENFSTLIBS) $(ATLASLIBS) -lm -ldl
CC = /usr/bin/arm-none-eabi-c++
CXX = /usr/bin/arm-none-eabi-c++
AR = /usr/bin/arm-none-eabi-ar
AS = /usr/bin/arm-none-eabi-as
RANLIB = ranlib

CXXFLAGS += -DHAVE_SPEEX -I/home/longqi/kaldi-trunk/src/../tools/speex/include
LDLIBS += -L/home/longqi/kaldi-trunk/src/../tools/speex/lib -lspeex
LDFLAGS += -Wl,-rpath=/home/longqi/kaldi-trunk/src/../tools/speex/lib
