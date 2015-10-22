#!/bin/bash
cd src/
make clean
make -j4
cd ../egs/yesno/s5
./decode.sh
less exp/mono0a/decode_test_yesno/log/decode.1.log 

