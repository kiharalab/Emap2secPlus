# Emap2sec+ is a computational tool using deep learning that can accurately identify structures, alpha helices, beta sheets, other(coils/turns) and DNA/RNA, in cryo-Electron Microscopy (EM) maps of medium to low resolution.
# Copyright (C) 2020 Xiao Wang, Eman Alnabati, Tunde W Aderinwale, Sai Raghavendra Maddhuri, Genki Terashi, Daisuke Kihara, and Purdue University.
# License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)
# Contact: Daisuke Kihara (dkihara@purdue.edu)


#

# This program is free software: you can redistribute it and/or modify

# it under the terms of the GNU General Public License as published by

# the Free Software Foundation, version 3.

#

# This program is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of

# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

# GNU General Public License V3 for more details.

#

# You should have received a copy of the GNU v3.0 General Public License

# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.en.html.
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
	cd gzstream; make -s
	mkdir -p $(BIN_DIR)
	$(CXX)  	main.cpp -o $(BIN_DIR)/map2train $(OPTIM) $(XFLAGS) $(LFLAGS) 
	$(CXX) -g -pg  main.cpp -o $(BIN_DIR)/map2train_debug $(XFLAGS) $(LFLAGS) 

clean:
	rm -f $(OBJ_FILES)
	rm -rf $(BIN_DIR)
	cd gzstream; make -s cleanall

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
