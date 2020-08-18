
# GCC 4.8+ assumed

SRC_DIR	= .
BIN_DIR = bin

CXX		= g++
XFLAGS	= -std=c++11 -march=native -Wall -fopenmp
LFLAGS	= -Lgzstream -lgzstream -lz
#OPTIM	= -O3 
OPTIM	= -O1

OBJ_FILES = $(wildcard *.o)

default:	map2train

all:	map2train

map2train:	main.cpp Grid.hpp
	cd gzstream; make
	mkdir -p $(BIN_DIR)
	$(CXX)  	main.cpp -o $(BIN_DIR)/map2train $(OPTIM) $(XFLAGS) $(LFLAGS) 
	$(CXX) -g -pg  main.cpp -o $(BIN_DIR)/map2train_debug $(XFLAGS) $(LFLAGS) 

clean:
	rm -f $(OBJ_FILES)
	rm -rf $(BIN_DIR)
	cd gzstream; make cleanall

package:	map2zernike
	mkdir -p map2zernike
	cp ./bin/map2zernike ./map2zernike
	cp ./README ./map2zernike
	cp ./LICENSE ./map2zernike
	tar czf map2zernike-distro-bin.tar.gz ./map2zernike
	rm ./map2zernike/map2zernike
	cp $(wildcard *.h) $(wildcard *.hpp) Makefile main.cpp ./map2zernike
	cp -r ./gzstream ./map2zernike/gzstream
	cp -r ./bin ./map2zernike/bin
	tar czf map2zernike-distro-full.tar.gz ./map2zernike
	rm -rf ./map2zernike
	@echo
	@echo "Created packages 'map2zernike-distro-bin.tar.gz' and 'map2zernike-distro-full.tar.gz'."
