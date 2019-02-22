rem
rem This script creates the mysolver.def...
rem
rem
cl /c /Od /MD /EHsc /Z7 /Fo"motsub.obj" "motsub.c"
rem
set OPATH=%PATH%
PATH=D:\MSC.Software\Adams\2017_2\utils\win64;%OPATH%
echo EXPORTS > mysolver.def
%PYTHON_EXE% %SHORT_TOPDIR%utils\bld_gendef.py -e -o tmp.def "motsub.obj"
type tmp.def >> mysolver.def
del tmp.def
set PATH=%OPATH%
set OPATH=
