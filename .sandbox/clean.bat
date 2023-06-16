@echo off
set folder="build\"
if not exist %folder% (
   echo "no build folder. creating..."
   mkdir build
)

cd /d %folder%
for /F "delims=" %%i in ('dir /b') do (rmdir "%%i" /s/q || del "%%i" /s/q) >nul 2>&1
copy NUL .keep
cd ..