@ECHO OFF
SETLOCAL
REM PLEASE CHANGE THE VALUE OF JROOT TO THE CORRECT PATH ON YOUR SYSTEM
SET JROOT="c:\users\frug\j602"

REM DO NOT CHANGE ANYTHING BELOW THIS LINE!
SET jcon="%JROOT%\bin\jconsole.exe
SET PATH=%JROOT%\bin;%PATH%
REM If no clustal alignment file supplied show usage:
if "%1"=="" goto usage
SET CLUSTAL_FILE=%1
if "%2"=="" goto usage
SET PROPERTIES_FILE=%2
if "%3"=="" goto norefseq
SET REF_SEQ=%3
goto readalignment
:norefseq
SET REF_SEQ=longest
:readalignment
REM Remove any previous intermediate file
if exist alignment.im del /Q alignment.im
if exist run.log del /Q run.log
REM Create the intermediate file
python read_alignment.py %CLUSTAL_FILE%
jconsole  snplist.ijs alignment.im %PROPERTIES_FILE%  %REF_SEQ%
goto finish
:usage
echo USAGE:
echo snplist clustalfilename propertiesfilename referencesequencename
echo or ( use longest sequence as reference sequence ):
echo snplist clustalfilename propertiesfilename 
goto finish
:finish
