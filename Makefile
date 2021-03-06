CC=gcc
CFLAGS = -c -O3 -Wall -msse -msse2 
LDFLAGS = -lz -lm 
SOURCES = pairExtractor.c  
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = pairExtractor


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 
clean:
	rm -f *.o *~ \#* pairExtractor
