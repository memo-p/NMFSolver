CXX=g++

LIBDIR  =  -L~/usr/lib  -L/usr/local/lib 

DEBUG=
OPTIM=-O2
# BLASOARMA=
BLASOARMA=-DARMA_DONT_USE_WRAPPER -framework Accelerate -lblas -llapack 
CFLAGS= $(DEBUG) $(OPTIM) -std=c++11 -fpermissive -Wall -w $(BLASOARMA)

LIBS= $(BLASOARMA) -lm -larmadillo -w -lhdf5

NAME=NMFSolver


libUnix: src/NMFSolver.o
	/bin/rm -f ./lib/lib$(NAME).a
	ar -cvq ./lib/lib$(NAME).a src/NMFSolver.o

src/%.o: src/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

main: src/NMFSolver.o src/main.o
	$(CXX) -o $(NAME) src/NMFSolver.o src/main.o $(CFLAGS) $(LIBS)
# main: $(OBJ) $(OBJDIR)/main.o
# 	$(CXX) -o $(NAME) $(OBJ) $(OBJDIR)/main.o $(CFLAGS) $(LIBS)

run: main
	./NMFSolver

.PHONY: clean

clean:
	rm -f src/*.o *~ ./lib/lib$(NAME).a ./NMFSolver






