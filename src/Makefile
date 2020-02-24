CC=gcc
CCFLAGS=-c -Wall -ansi -lm
LDFLAGS=-lm

TARGET=example
MAIN=$(addsuffix .o, $(TARGET))
SOURCES=pda/pda.c
OBJ=$(SOURCES:.c=.o) $(MAIN)
DEPS=pda/pda.h pda/_pda.h

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

$(OBJ): %.o : %.c $(DEPS)
	$(CC) -c -o $@ $< $(CCFLAGS)

clean:
	rm -f $(TARGET) $(OBJ)
