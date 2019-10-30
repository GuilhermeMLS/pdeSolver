FLAGS=  -O3 -mavx -march=native
DEBUG= -Wall -g -DDEBUG 
CXX=gcc
RM=rm -rf

all: main.c pdelib.o  Makefile
	$(CXX) -o pdeSolver main.c pdelib.o $(FLAGS)

pdelib.o: pdelib.c  pdelib.h Makefile
	$(CXX) -c pdelib.c pdelib.h $(FLAGS)

clean:
	$(RM) pdeSolver *.o

doc: $(OBJ)
	doxygen ./Doxyfile
#doc: *.c trabalho-1.doxy *.h
#	@doxygen trabalho-1.doxy