# Makefile for detect_oscillation program (Windows)

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -std=c99

# Source and target
SRC = detect_oscillation.c
TARGET = detect_oscillation.exe

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $<

# Clean build files
clean:
	del /Q $(TARGET)

# Run the program
run: $(TARGET)
	$(TARGET)

.PHONY: all clean run 