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
#ifndef ZernikeDescriptor_h
#define ZernikeDescriptor_h

#include <vector>
#include <complex>

#include "GeoMoments.h"
#include "ZernikeMoments.h"
#include "Grid.h"
#include "util.h"

using std::vector;
using std::complex;

template<typename T>
class ZernikeDescriptor
{
public:
	typedef vector<T> T1D;
	typedef typename T1D::iterator T1DIter;


	ZernikeDescriptor ();
	ZernikeDescriptor(Grid<T> &grid, int order);

	typename vector<T>::iterator begin();
	typename vector<T>::iterator end();
	size_t size() const;
	int getOrder() const;
 
private:
	int maxOrder;
	T1D invariants;
};



#include "ZernikeDescriptor.hpp"

#endif


