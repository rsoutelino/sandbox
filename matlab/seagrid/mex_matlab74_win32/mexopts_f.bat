@echo off
rem INTELF91MSVS2005OPTS.BAT
rem
rem    Compile and link options used for building MEX-files using the
rem    Intel� Fortran Compiler 9.1 with the Microsoft� Visual Studio�
rem    2005 linker.
rem
rem    $Revision: 1.1.6.1 $  $Date: 2006/11/19 18:53:36 $
rem
rem ********************************************************************
rem General parameters
rem ********************************************************************
set MATLAB=%MATLAB%
set IFORT_COMPILER91=C:\Program Files\Intel\Compiler\Fortran\9.1
set VS80COMNTOOLS=%VS80COMNTOOLS%
set LINKERDIR=C:\Program Files\Microsoft Visual Studio 8
set PATH=%IFORT_COMPILER91%\IA32\Bin;%LINKERDIR%\VC\BIN;%LINKERDIR%\Common7\Tools;%LINKERDIR%\Common7\Tools\bin;%LINKERDIR%\Common7\IDE;%LINKERDIR%\SDK\v2.0\bin;%PATH%
set INCLUDE=%IFORT_COMPILER91%\IA32\Include;%LINKERDIR%\VC\ATLMFC\INCLUDE;%LINKERDIR%\VC\INCLUDE;%LINKERDIR%\VC\PlatformSDK\include;%LINKERDIR%\SDK\v2.0\include;%INCLUDE%
set LIB=%IFORT_COMPILER91%\IA32\Lib;%LINKERDIR%\VC\ATLMFC\LIB;%LINKERDIR%\VC\LIB;%LINKERDIR%\VC\PlatformSDK\lib;%LINKERDIR%\SDK\v2.0\lib;%MATLAB%\extern\lib\win32;%LIB%
set MW_TARGET_ARCH=win32

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=ifort
set COMPFLAGS=/fpp /Qprec "/I%MATLAB%/extern/include" -c -nologo -DMATLAB_MEX_FILE /fixed
set OPTIMFLAGS=/MD -Ox -DNDEBUG /real_size:64
set DEBUGFLAGS=/MD -Zi /PDB:"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win32\microsoft
set LINKER=link
set LINKFLAGS=/DLL /EXPORT:MEXFUNCTION /MAP /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /implib:%LIB_NAME%.lib /NOLOGO
set LINKOPTIMFLAGS=
set LINKDEBUGFLAGS=/debug /PDB:"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT="/out:%OUTDIR%%MEX_NAME%%MEX_EXT%"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=

set POSTLINK_CMDS=del "%OUTDIR%%MEX_NAME%.map"
set POSTLINK_CMDS1=mt -outputresource:"%OUTDIR%%MEX_NAME%%MEX_EXT%";2 -manifest "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
set POSTLINK_CMDS2=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest" 
