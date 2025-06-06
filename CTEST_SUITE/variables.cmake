# Variables to set for Jean-Zay
#set(ENV{CC} "icc")
#set(ENV{CXX} "icpc")
#set(ENV{FC} "ifort")
#set(SFEMaNS_DIR "/gpfswork/rech/nor/commun/SFEMaNS_v5.4/SFEMaNS")
#set(debug_bounds "$ENV(FFLAGS) -O2 -g -traceback -heap-arrays")
#set(release_bounds "-O3")

# Variables to set for Whistler
set(SFEMaNS_DIR "/mnt/c/Users/MELVIN/PycharmProjects/SFEMaNS_GIT")
set(ADDITIONAL_LINKS "-L /home/creff/petsc/lib")
set(debug_bounds "-Wall -fimplicit-none -fbounds-check -fallow-argument-mismatch")
set(release_bounds "-O3 -fallow-argument-mismatch")
set(native_bounds "-march=native -mtune=native -Ofast")
