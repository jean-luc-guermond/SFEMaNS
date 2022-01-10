// ---------------------------------------------------------------------
// $Id$
//
//    This file is part of SFEMaNS.
//
//    SFEMaNS is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    SFEMaNS is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with SFEMaNS.  If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------



/**
 * @page doc_installation Installation Instructions



@section doc_install_sofware Installation of required software

@subsection doc_install_petsc PETSc for parallel computing (required)

SFEMaNS uses the sofware PETSc to do the linear algebra in parallel. The following documentation describes the procedure to install PETSc. We refer to  <a href='https://www.mcs.anl.gov/petsc/download/index.html'><code>https://www.mcs.anl.gov/petsc/download/index.html</code></a> for more information about PETSc. Since it is more efficient to use the native PETSc of your cluster, you may want to ask your root administrator to install it for you. If you do so, you need to specify that you want the version petsc-3.15.1 and that it needs to include blaslapack, hypre, mumps, scalapack, metis, parmetis and blacs.

Here is the procedure to install PETSc by yourself:
<ol>
<li>First go to the directory where you want to install petsc and dowload it with the command:
\code
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.15.1.tar.gz
\endcode
<li>Untar the file with the following command:
\code
tar -xvf petsc-3.15.1.tar.gz
\endcode
This creates a directory named <tt>petsc-3.15.1</tt>. 
<li>Go into the directory petsc-3.15.1 and defined the following environnement variable:
\code
cd petsc-3.15.1
export PETSC_DIR=$PWD
export PETSC_ARCH=linux-gnu-c
\endcode
Note that depending of the work environment, the command <tt>export</tt> can be replaced by <tt>setenv</tt> or equivalent.
<li>Configure petsc with the command:
\code
./configure --configModules=PETSc.Configure --optionsModule=config.compilerOptions --download-f-blaslapack=1 --with-shared-libraries=1 --download-hypre=1 --download-mumps=1 --download-scalapack=1 --download-metis=1 --download-parmetis=1 --download-blacs=1 --with-debugging=0 --with-x=0
\endcode
<li>Then type the command:
\code
make all
make check
\endcode
</ol>

The installation of petsc-3.15.1 is now complete.

@subsection doc_install_fftw FFTW for Fast Fourier Transform (required)

To compute nonlinear terms efficiently, SFEMaNS uses the package FFTW (Fast Fourier Transform in the West).
If it is not yet installed in your work environment, either ask your root administrator to install it or follow the instructions in <a href='http://www.fftw.org/download.html'><code>http://www.fftw.org/download.html</code> </a>.

@subsection doc_install_arpack ARPACK and PARPACK for eigen value problem (optional)

SFEMaNS can solve eigenvalue problem using the software ARPACK. If you do not plan to solve such problems, you don't need to install ARPACK. 
Notice that to solve eigenvalue problems, you will be required to uncomment a few lines in the file <tt>main.f90</tt> as describeb in the next section.

The following describes the procedure to install ARPACK by yourself. We refer to  <a href='http://www.caam.rice.edu/software/ARPACK/'><code>http://www.caam.rice.edu/software/ARPACK/</code> </a> for more information on ARPACK and its installation.

<ol>
<li>First download ARPACK and PARPACK:
\code
wget http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz
wget http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz
wget http://www.caam.rice.edu/software/ARPACK/SRC/parpack96.tar.gz
wget http://www.caam.rice.edu/software/ARPACK/SRC/ppatch.tar.gz
\endcode
<li>Extract these files with the following commands:
\code
tar -xvf arpack96.tar.gz
tar -xvf patch.tar.gz
tar -xvf parpack96.tar.gz
tar -xvf ppatch.tar.gz
\endcode
This creates a directory called ARPACK.
<li>Go into this directory:
\code
cd ARPACK
\endcode
<li>Edit the file <tt>ARmake.inc</tt>. Here are the main steps to follow.
<ol>
<li>To start with a suitable <tt>ARmake.inc</tt>, you need to know which environment you are using (SUN4, SP2, ...).
<li> Copy the suitable ARmake.inc from the directory <tt>ARMAKES</tt>. For example for sun4, type the command:
\code
cp ARMAKES/ARMAKES/ARmake.MPI-SUN4 ARmake.inc
\endcode
<li>Open the file <tt>ARmake.inc</tt>.
<li>Edit the path to the ARPACK directory (variable called home).
<li>Check that the compilation program and flags are the correct one. For example, set the compilers FC and PFC to mpif90, set the FFLAGS and PFFLAGS to -O, set the correct path to the command make (type "which make" in your terminal to know its path). The above list is not exhaustive, more information are provided by ARPACK's developers in the file README of the directory ARPACK.
</ol>
<li>Type the command:
\code
make lib
make plib
\endcode
</ol>

The installation of arpack and parpack is now complete.

@section doc_install_sfemans Dowload and installation of SFEMaNS

In this section, we describe how to download SFEMaNS and how to get an executable called <tt>a.exe</tt>. We also provide information on how to check that SFEMaNS and the above sofwares are installed correctly.

@subsection doc_install_sfemans_download Download

Download SFEMaNS_v5.3 with the following command:
\code
wget  http://math.tamu.edu/~guermond/DOWNLOADS/SFEMaNS_v5.3.tar.gz
\endcode

Then move the file <tt>SFEMaNS_v5.3.tar.gz</tt> where you want to install the code SFEMaNS. Untar it with the command:
\code
tar -xvf SFEMaNS_v5.3.tar.gz
\endcode

Define the following environment variable in your <tt>.bashrc</tt>:
\code
export SFEMaNS_DIR= the_path_you_chose_for_SFEMaNS
\endcode
Note that depending of your shell environment, the command <tt>export</tt> can be replaced by <tt>setenv</tt> or equivalent.

IMPORTANT: this last step is required as the code SFEMaNS needs to know this variable when compiling the executable. As a consequence, we ask every user to add this line in his <tt>.bashrc</tt> or equivalent.

@subsection doc_install_sfemans_compilation  Compilation of the executable a.exe

To generate an executable of the code SFEMaNS, the following files are needed:
<ol>
<li><tt>make.inc</tt> that defines environment variable (path to PETSc, FFTW, ARPACK) and compilation options.
<li><tt>my_make</tt> that allows to generate a file called <tt>makefile</tt>. This file contains the tree of the code (which module need which modules, etc.).
<li><tt>main.f90</tt> that is used to set the desired outputs with the subroutine <tt>my_post_processing</tt>.
<li><tt>condlim.f90</tt> that is used to set the initial conditions, boundary conditions and source term of the problem considered.
<li><tt>read_user_data.f90</tt> that is used to add new variable in the code. For instance, such variable can be used in the condlim to change the amplitude of a forcing term. 
</ol>

Templates of each of the above files are available in the following directory: <tt>$SFEMaNS_DIR/TEMPLATE</tt>. We refer to the section \ref doc_computation_with_SFEMaNS for more details on how to generate your own executable and set properly your <tt>data</tt> file.


@subsection doc_install_sfemans_check_install  Check Installation

The code SFEMaNS presents more than 30 tests that mainly involve manufactured solutions. They have been implemented so one can check the correct installation of the code. It also allows developpers to check the consistency of the new features they implement. To run these test, follows the intructions below:
<ol>
<li>Create a directory for your applications with SFEMaNS. For instance, type the command:
\code
mkdir MY_APPLICATIONS_SFEMaNS_v5.3
\endcode
<li>Go into this directory and type the following command:
\code
cp -rf $SFEMaNS_DIR/TEMPLATE CHECK_INSTALLATION
\endcode
<li>Go in the directory CHECK_INSTALLATION.
\code
cd CHECK_INSTALLATION
\endcode
<li>If ARPACK is installed, you can uncomment the following sections in the file <tt>main.f90</tt>:
<ol>
<li> Line 6
\code
    USE arpack_mhd
\endcode
<li> The section "EIGENVALUE PROBLEMS/ARPACK" (lines 72 to 91)
</ol>
<li>Rename the file <tt>make.in_template</tt> and <tt>my_make_template</tt> as follows:
\code
mv make.inc_template make.inc
mv my_make_template my_make
\endcode
<li>Edit the file <tt>make.inc</tt>.
<ol>
<li>Set RACINE_FFTW. It is the directory that contains include/fftw3.h. To know this path, type the following command:
\code
 locate fftw3.h
\endcode
The result is RACINE_FFTW/include/fftw3.h.
<li>Set LIB_FFTW. It is the result of the command:
\code
locate libfftw3.so
\endcode
<li>Set PETSC_DIR (see compilation of PETSc).
<li>Set PETSC_ARCH (see compilation of PETSc).
<li>Set the preprocessor option: PREPROC.
<li>Set the compilation option OPT (-O3, -C, etc.).
<li>If ARPACK is installed, set PA_LIB as follows:
\code
PA_LIB= path_to_arpack/parpack_MPI-SUN4.a path_to_arpack/libarpack_SUN4.a
\endcode
Note that the the name of the files <tt>parpack_MPI-SUN4.a</tt> and <tt>libarpack_SUN4.a</tt> needs to be change for system other than sun4.
<li>All the other variables should not require any modification.
</ol>
<li>Generate a makefile with the command:
\code
./my_make
\endcode
<li>Compile the executable a.exe with the command:
\code
make a.exe -j 8
\endcode
The option -j 8 means 8 processors are used to do the compilation. It can be removed or used with a different number of processors.
<li>Edit the shell <tt>debug_SFEMaNS_template</tt> as follows:
<ol>
<li>Set the first and last test numbers to check with the variables iter_beg and iter_end. As SFEMaNS_v5.3 presents 40 tests, we have \f$1 \leq \text{iter_beg} \leq \text{iter_end} \leq 40\f$. In order to check all the tests, set iter_beg to 1 and iter_end to 40.
<li>Define the program used to run MPI program (mpirun, mpiexec). It is done by setting the variable MY_MPIRUN as follows:
\code
export MY_MPIRUN=mpirun
\endcode
</ol>
<li>Run the shell <tt>debug_SFEMaNS_template</tt> with the following command:
\code
./debug_SFEMaNS_template
\endcode
 For each test, a "OK" should be displayed. It means the test is performed correctly. 
 We note that the test 14 is disabled by default because it requires ARPACK.
 The ouput of each test, which can be found in the files called fort.10_T*, are compared with reference values. These reference values are at the end of the files $SFEMaNS_DIR/MHD_DATA_TEST_CONV_PETSC/debug_data_* where * is the number of the test considered. 
</ol>

 */
