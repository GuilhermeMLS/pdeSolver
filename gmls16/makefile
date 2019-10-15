FLAGS=  -Wall 
DEBUG= -Wall -g -DDEBUG 
CXX=gcc
RM=rm -rf

all: main.c pdelib.o  Makefile
	$(CXX) -o pdeSolver main.c pdelib.o $(FLAGS)

pdelib.o: pdelib.c  pdelib.h Makefile
	$(CXX) -c pdelib.c pdelib.h $(FLAGS)

clean:
	$(RM) pdeSolver *.o