# Variables to set for Jean-Zay

set(ENV{CC} "icc")
set(ENV{CXX} "icpc")
set(ENV{FC} "ifort")
set(SFEMaNS_DIR "/gpfswork/rech/nor/commun/SFEMaNS_v5.4/SFEMaNS")
#set(COMPILE_OPTIONS "-fpp")
# Variables to set for Whistler
=======
>>>>>>> Stashed changes
#set(ENV{CC} "icc")
#set(ENV{CXX} "icpc")
#set(ENV{FC} "ifort")
#set(SFEMaNS_DIR "$ENV{SFEMaNS_DIR}")
#set(COMPILE_OPTIONS "-x f95-cpp-input")


# Variables to set for Whistler
set(SFEMaNS_DIR "/home/guermond/SFEMaNS/GIT_SFEMaNS/SFEMaNS")
set(COMPILE_OPTIONS "-x f95-cpp-input")
set(ADDITIONAL_LINKS "-lmetis -lz -L /usr/lib/x86_64-linux-gnu/hdf5/serial")
