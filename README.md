# HS_HSSW_quench_NESCGLE
A code to compute the mobility function <img src="https://latex.codecogs.com/svg.latex?\;b(t;\phi,T_f,T_i=\infty)" title="b(t,phi,Tf,Ti)" /> of a quench for the Hard Sphere (HS) + Square Well system starting from an infinite initial temperature (HS limit) to a finite final temperature through the "Non-equilibrium self-consistent generalized Langevin equation" formalism.

<h2>Installation guide:</h2>

This program is written in Fortran 2003, make sure to have a fortran compiler and GNU make utility.
The make file does utilize gfortran as main command to compile (line 4 of Makefile), make sure to modify it to
accordingly to your compiler. gcc options -ffast-math -o3 are used but not necessary (lines 6 and 7 of Makefile).

Tested with gcc version 7.1.0

To install simply execute the command "make" while being in this folder.

<h3>Tools website:</h3>

make utility gnu website: https://www.gnu.org/software/make/
gcc compiler website: https://gcc.gnu.org/

<h3>Tools installation:</h3>

For GNU OS systems make sure to have the apt utility installed in your distribution. If you do not have an utility, then simply execute
in your terminal for the ones you are missing: 

"apt-get make" 
"apt-get gcc"
"apt-get gfortran"

internet connection is needed. If you do not have access to the apt utility then look at the official websites of the tools for installation guides.

***********************************

For Windows user there are several alternatives. For beginner users I recommend mingw-w64 installer option which can be downloaded at: 
https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe/download
which is a friendly windows installer for gcc and make. Please note that the commands will only be loaded in the final mingw-64 terminal 
before installation. The "make" utility is usually name "mingw32-make" which can be found inside the binaries folder "bin" inside the mingw-64
installation folder along "gfortran" binary. If you wish to export all the binaries to the windows "command box" cmd, make sure to add the path of
the binaries to the user or system "path" of environment variables (search for "add path to environment  variables in windows" at any search engine).

mingw-w64 official website: http://mingw-w64.org/doku.php

*********************************

For Mac beginner users I recommend installing the "Homebrew" installer utility. Please make sure to have Xcode, which can be found at the apple 
store (it is free) and Xcode Command Line Tools. Then follow the instructions located at https://brew.sh (Hombrew installer utility main website).
Once the installer utility has been installed simply execute

"brew install gcc"
"brew install --cask gfortran"

Please make sure that the command gfortran is indeed a gfortran compiler (one can check the version of gfortran with "gfortran -v"). There are some cases 
in which OS X have a preinstalled fortran compiled which is also executed with "gfortran" command but it is not gfortran and may have issues with the
programming lines used in the fortran code lines if it does not support fortran 2003 syntaxis. If this is the case, one can make sure to call the correct
"gfortran" by specifying the version i.e. "gfortran-8", with which you will need to modify the line 4 of the Makefile file.



<h2>Execution guide:</h2>

To execute the program simply run the file "program.out" generated by "make" inside this folder. 

Note: GNU OS users do have to specify the folder in which "program.out" is located. In the terminal, 
if you are currently within this folder, the current path can be added by "./", hence you can execute 
the program by running "./program.out" in your terminal. Windows users may simply execute "program.out"
in the cmd when inside this folder. For OSX (Mac) users, follow the instruction of GNU OS users.
