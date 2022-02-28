CC=+gg
NVCC=nvcc
CUDAFLAGS= -std=c++11 -c




all: gravitation ;

gravitation: main.o engine.o body.o tree_node.o tree.o ConfigFile.o
	g++ -std=c++11 -o gravitation main.o engine.o body.o tree_node.o tree.o ConfigFile.o

main.o: main.cpp tree.h body.h engine.h
	g++ -std=c++11 -c main.cpp

body.o: body.cpp body.h
	g++ -std=c++11 -c body.cpp

tree_node.o: tree_node.cpp tree_node.h physical_constants.h body.h
	g++ -std=c++11 -c tree_node.cpp

tree.o: tree.cpp tree.h tree_node.h body.h
	g++ -std=c++11 -c tree.cpp

engine.o: engine.cpp engine.h body.h tree.h ConfigFile.h
	g++ -std=c++11 -c engine.cpp

ConfigFile.o: ConfigFile.h ConfigFile.cpp
	g++ -std=c++11 -c ConfigFile.cpp

clean:
	rm -f *.o

distclean: clean
	rm -f gravitation
