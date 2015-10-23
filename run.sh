#!/bin/bash
cd src/
make -j4
cd ../egs/yesno/s5
./decode.sh

