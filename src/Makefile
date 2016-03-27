default: example
example: example.o pda.o
	gcc example.o pda.o -lm -ansi -o example
example.o: example.c
	gcc -c example.c -lm -ansi
pda.o: pda/pda.c
	gcc -c pda/pda.c -lm -ansi
clean:
	rm -f pda.o example.o example
