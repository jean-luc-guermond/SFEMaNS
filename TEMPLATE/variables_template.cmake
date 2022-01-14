# Variables to set for Jean-Zay
#set(ENV{CC} "icc")
#set(ENV{CXX} "icpc")
#set(ENV{FC} "ifort")
#set(SFEMaNS_DIR "/gpfswork/rech/nor/commun/SFEMaNS_v5.4/SFEMaNS")
#set(debug_bounds "$ENV(FFLAGS) -O2 -g -traceback -heap-arrays")
#set(release_bounds "-O3")

# Variables to set for Whistler
set(SFEMaNS_DIR "/home/guermond/SFEMaNS/GIT_SFEMaNS/SFEMaNS")
set(ADDITIONAL_LINKS "-lmetis -lz -L /usr/lib/x86_64-linux-gnu/hdf5/serial")
set(debug_bounds "-Wall -fimplicit-none -fbounds-check")
set(release_bounds "-O3")
set(native_bounds "-march=native -mtune=native -Ofast")
