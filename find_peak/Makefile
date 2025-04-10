CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99
LDFLAGS = -lm

all: test_peaks

find_peaks.o: find_peaks.c find_peaks.h
	$(CC) $(CFLAGS) -c find_peaks.c -o find_peaks.o

test_peaks.o: test_peaks.c find_peaks.h
	$(CC) $(CFLAGS) -c test_peaks.c -o test_peaks.o

test_peaks: test_peaks.o find_peaks.o
	$(CC) $(CFLAGS) test_peaks.o find_peaks.o -o test_peaks $(LDFLAGS)

clean:
	rm -f *.o test_peaks *.csv

.PHONY: all clean 