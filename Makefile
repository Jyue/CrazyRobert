CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99
LDFLAGS = -lm

all: find_peak_test

find_peak_test: find_peak.o find_peak_test.o
	$(CC) -o find_peak_test find_peak.o find_peak_test.o $(LDFLAGS)

find_peak.o: find_peak.c find_peak.h
	$(CC) -c find_peak.c $(CFLAGS)

find_peak_test.o: find_peak_test.c find_peak.h
	$(CC) -c find_peak_test.c $(CFLAGS)

test: find_peak_test
	./find_peak_test
	python3 compare_results.py

clean:
	rm -f *.o find_peak_test *.csv compare_results.py 