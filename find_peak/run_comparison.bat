@echo off
echo Running validation tests between Python and C versions

echo.
echo 0. Installing Python dependencies...
call install_dependencies.bat
if %ERRORLEVEL% NEQ 0 (
    echo Failed to install Python dependencies. Exiting.
    exit /b 1
)

echo.
echo 1. Compiling C version...
call compile.bat
if %ERRORLEVEL% NEQ 0 (
    echo Failed to compile C program. Exiting.
    exit /b 1
)

echo.
echo 2. Running C test program...
test_peaks.exe
if %ERRORLEVEL% NEQ 0 (
    echo Failed to run C test program. Exiting.
    exit /b 1
)

echo.
echo 3. Running Python test program...
python test_python.py
if %ERRORLEVEL% NEQ 0 (
    echo Failed to run Python test program. Exiting.
    exit /b 1
)

echo.
echo 4. Comparing results...
python compare_results.py
if %ERRORLEVEL% NEQ 0 (
    echo Failed to compare results. Exiting.
    exit /b 1
)

echo.
echo Done! Check the comparison results and generated images.
echo. 