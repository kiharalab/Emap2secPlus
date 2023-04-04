/*
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
*/

#include <iostream>
#include <istream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include "pdb.h"

#ifndef Grid_h

using std::istream;
using std::ifstream;
using std::vector;

extern bool VERBOSE;

enum {TYPE_UNKNOWN, TYPE_SITUS, TYPE_MRC};


template <typename Value=double>
class Grid{

	protected:
		Value defaultValue = 0;
		Value indicatorValue = 1000;
		vector<Value> vec;
		
		int xdim, ydim, zdim; // for some reason the file format allows signed dims
		size_t numVoxels;

		void readMRC(istream &in);
		void readSitus(istream &in);
		void interpSkewToCubic(double, double, double, double, double, double);
		void interpRectToCubic(double, double, double);


	public:
		Value stepx,stepy,stepz,origx,origy,origz;
		Value origx2,origy2,origz2;//Consider ncstart,nrstart...
		int ActV;
		//vector<atom> a;

		enum{ CONTOUR_INPLACE, CONTOUR_COPY };

		Grid(istream &in, int type);
		Grid(istream &in, int type, Value def); // create with a nonzero default default value
		Grid(Grid &g, Value contour); // get the map at a certain contour
		//Grid(Grid &g, Value contour,Vec3 pos); // get the map at a certain contour & position
		Grid(Grid &g, Value contour,double *,double); // get the map at a certain contour & position
		Grid(Grid &g, Value contour,bool);//Density Mode
		Grid(Grid &g, Value contour,vector<atom> &,double,int,int,int,bool,bool);//Genki
		Grid(Grid &g, Value contour,double,int,int,int,bool);//

		void readFrom(istream &in, int type);

		void set(size_t x, size_t y, size_t z, Value val);
		Value get(size_t x, size_t y, size_t z); // exception on out-of-bounds
		Value getOrDefault(size_t x, size_t y, size_t z); // default to default default
		Value getOrDefault(size_t x, size_t y, size_t z, Value def); // default to supplied default

		double pos_dens(Grid &g,double *cd);

		size_t sizeX() const;
		size_t sizeY() const;
		size_t sizeZ() const;
		size_t size() const;
		size_t sizeMax() const;

		void contour(Value contour); // in-place

		typename std::vector<Value>::iterator begin();
		typename std::vector<Value>::iterator end();
		typename std::vector<Value>::iterator begin(size_t k, size_t n);
		typename std::vector<Value>::iterator end(size_t k, size_t n);
};

#include "Grid.hpp"

#define Grid_h
#endif
