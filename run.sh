#!/bin/bash
KALDI_ROOT=`pwd`
echo Kaldi root is $KALDI_ROOT
cd src/
make KALDI_ROOT=$KALDI_ROOT -j4
cd ../egs/yesno/s5
./decode.sh

