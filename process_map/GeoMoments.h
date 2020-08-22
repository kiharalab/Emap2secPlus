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
#ifndef GeoMoments_h
#define GeoMoments_h

#include <vector>
#include <stdexcept>
#include "Grid.h"
#include "Sum.h"
#include "util.h"

using std::vector;

template<typename T>
class GeoMoments
{
public:
	typedef vector<T> T1D;
	typedef NDVector<T, 3> T3D;

	typedef typename T1D::iterator T1DIter;

	GeoMoments();
	GeoMoments(Grid<T> &grid, const Vec3<T> &center = {0,0,0}, T scale = 1.0, int maxOrder = 1);

	T getMoment(int i, int j, int k);

private:	
	Vec3<T1D> samples;
	int maxOrder;
	NDVector<T, 3> moments;

	void computeMoments(Grid<T> &grid);
	void computeFirst(Grid<T> &grid);
	void computeSamples(const Grid<T> &grid, const Vec3<T> &center, T scale);
	void computeDiffs(T1DIter iter, T1DIter diffIter, int dim);
	void computeDiffsInitial(T1DIter iter, T1DIter diffIter, int dim, T &extra);
	
	T multiply(T1DIter diffIter, T1DIter sampleIter, int dim);
	T multiplyInitial(T1DIter diffIter, T1DIter sampleIter, int dim, T &extra);
};

#include "GeoMoments.hpp"

#endif

