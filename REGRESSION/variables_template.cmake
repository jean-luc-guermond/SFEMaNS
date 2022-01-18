# Variables to set for Jean-Zay
#set(ENV{CC} "icc")
#set(ENV{CXX} "icpc")
#set(ENV{FC} "ifort")
#set(debug_bounds "$(FFLAGS) -O2 -g -traceback -heap-arrays")
#set(release_bounds "-O3")
#set(SFEMaNS_DIR "/gpfswork/rech/nor/commun/SFEMaNS_v5.4/SFEMaNS")
#set(RUN_PRE_PROC "srun")
#set(PROC_CALL "--ntasks=")
#set(RUN_POST_PROC "--hint=nomultithread --job-name=regression_SFEMaNS --time=00:02:00 --qos=qos_cpu-dev -A nor@cpu")

# Variables to set for Whistler
set(SFEMaNS_DIR "/home/guermond/SFEMaNS/GIT_SFEMaNS/SFEMaNS")
set(ADDITIONAL_LINKS "-lmetis -lz -L /usr/lib/x86_64-linux-gnu/hdf5/serial")
set(debug_bounds "-Wall -fimplicit-none -fbounds-check")
set(release_bounds "-O3")
set(native_bounds "-march=native -mtune=native -Ofast")
set(RUN_PRE_PROC "srun")
set(PROC_CALL "-n")
set(RUN_POST_PROC "")
