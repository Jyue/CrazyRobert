@echo off
echo Checking system dependencies...

echo.
echo 1. Checking for GCC compiler...
where gcc >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo GCC not found. You need to install MinGW or another GCC compiler.
    echo You can download MinGW from https://sourceforge.net/projects/mingw/
    echo After installing, make sure to add the bin directory to your PATH.
    echo.
    echo If you have Visual Studio, you can use the MSVC compiler instead:
    echo Replace 'gcc' with 'cl' in compile.bat
) else (
    echo GCC found: 
    gcc --version | findstr GCC
)

echo.
echo 2. Checking for Python...
where python >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Python not found. Please install Python 3.
    echo You can download Python from https://www.python.org/downloads/
) else (
    echo Python found:
    python --version
)

echo.
echo 3. Checking for Python packages...
python -c "import numpy" >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo NumPy not found. Run install_dependencies.bat to install.
) else (
    echo NumPy found.
)

python -c "import scipy" >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo SciPy not found. Run install_dependencies.bat to install.
) else (
    echo SciPy found.
)

python -c "import matplotlib" >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Matplotlib not found. Run install_dependencies.bat to install.
) else (
    echo Matplotlib found.
)

python -c "import pandas" >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Pandas not found. Run install_dependencies.bat to install.
) else (
    echo Pandas found.
)

echo.
echo Dependency check complete. If any dependencies are missing, please install them before running the validation. 