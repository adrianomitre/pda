CC=gcc
CFLAGS=-c -Wall -ansi -lm
LDFLAGS=-lm
SOURCES=example.c pda/pda.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=example
all: $(SOURCES) $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)
.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -f pda/pda.o example.o example
