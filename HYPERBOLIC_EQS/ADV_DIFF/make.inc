ROOT           = ${SFEMaNS_DIR}

#######################################################
########Quarck family
#######################################################
RACINE_FFTW    = /usr
PETSC_DIR      = $(ROOT)/../petsc-3.10.3
PETSC_ARCH     = linux-gnu-c
I_PETSC        = -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include 
INC_FFTW       =-I$(RACINE_FFTW)/include
PREPROC        = -x f95-cpp-input
OPT            =  $(FFLAGS) -O3
#OPT            =  $(FFLAGS) -Wall -fimplicit-none -fbounds-check #-O3
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
F90	       = $(FC) $(PREPROC) $(OPT) $(I_PETSC) $(INC_FFTW)
LNK	       = $(FLINKER) 
#######################################################

#######################################################
########Faraday
#######################################################
#PREPROC        = -Mpreprocess
#OPT            = $(FFLAGS) -O3
#PETSC_DIR      = $(ROOT)/../petsc-3.4.0
#PETSC_ARCH     = linux-pgi
#include ${PETSC_DIR}/conf/variables
#F90            = $(FC) $(PREPROC) $(OPT) $(I_PETSC) 
#LNK            = $(FLINKER)
#I_PETSC        = -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include
#######################################################

#######################################################
########Stampede
#######################################################
#PETSC_DIR      = /opt/apps/intel13/mvapich2_1_9/petsc/3.4
#PETSC_ARCH     = 
#PREPROC        = -cpp 
#OPT            = -O3 
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules
#F90            = $(FC) $(PREPROC) $(OPT) $(I_PETSC) $(INC_FFTW)
#LNK            = $(FLINKER)
#I_PETSC        = -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include
#######################################################

#######################################################
########Ada
#######################################################
#RACINE_FFTW    =
#PREPROC        = -fpp
#LIB_FFTW       =
#PA_LIB         =
#OPT            = $(FFLAGS) -O2 -g -traceback -heap-arrays # -O3
#include ${PETSC_DIR}/${PETSC_ARCH}/conf/variables
#include ${PETSC_DIR}/${PETSC_ARCH}/conf/rules
#F90            = $(FC) $(PREPROC) $(OPT) $(I_PETSC) $(INC_FFTW)
#LNK            = $(FLINKER)
#I_PETSC        = -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include
#INC_FFTW       = -I$(RACINE_FFTW)/include
#######################################################
