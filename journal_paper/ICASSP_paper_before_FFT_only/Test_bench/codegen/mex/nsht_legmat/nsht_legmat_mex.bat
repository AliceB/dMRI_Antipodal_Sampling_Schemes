@echo off
set MATLAB=C:\PROGRA~1\matlab\r2013b
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\matlab\r2013b\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=nsht_legmat_mex
set MEX_NAME=nsht_legmat_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for nsht_legmat > nsht_legmat_mex.mki
echo COMPILER=%COMPILER%>> nsht_legmat_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> nsht_legmat_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> nsht_legmat_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> nsht_legmat_mex.mki
echo LINKER=%LINKER%>> nsht_legmat_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> nsht_legmat_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> nsht_legmat_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> nsht_legmat_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> nsht_legmat_mex.mki
echo BORLAND=%BORLAND%>> nsht_legmat_mex.mki
echo OMPFLAGS= >> nsht_legmat_mex.mki
echo OMPLINKFLAGS= >> nsht_legmat_mex.mki
echo EMC_COMPILER=msvcsdk>> nsht_legmat_mex.mki
echo EMC_CONFIG=optim>> nsht_legmat_mex.mki
"C:\Program Files\matlab\r2013b\bin\win64\gmake" -B -f nsht_legmat_mex.mk
