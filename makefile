
CC=clang
CFLAGS=-c -Wall -g

all: strassen

strassen: strassen.o 
	$(CC) strassen.o -o strassen

main.o: strassen.c
	$(CC) $(CFLAGS) strassen.c

clean:
	rm *.o strassen
