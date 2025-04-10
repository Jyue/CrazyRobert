@echo off
echo Compiling C program...

REM 檢查 gcc 是否可用
where gcc >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo gcc not found. Please install MinGW or another C compiler.
    exit /b 1
)

REM 編譯
gcc -Wall -c find_peaks.c -o find_peaks.o
gcc -Wall -c test_peaks.c -o test_peaks.o
gcc -Wall find_peaks.o test_peaks.o -o test_peaks.exe -lm

echo Compilation complete. 