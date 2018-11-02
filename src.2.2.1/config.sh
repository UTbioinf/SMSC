#!/bin/bash

# configure mummer root directory, full path is a must
mummer_path=$MUMMERDIR

# configure mhaps, the executable jar is included. Full path is a must
mhap_path=$LOONLOCAL/openbiosrc/mhap/2.1/mhap-2.1.jar
mhap_mem=32g

# configure muscle, the executable is included. Full path is a must
MUSCLE_EXEC=${LOONLOCAL}/bin/muscle
max_iters=5
