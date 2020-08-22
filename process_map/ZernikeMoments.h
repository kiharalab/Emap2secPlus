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
#ifndef ZernikeMoments_h
#define ZernikeMoments_h

#include <complex>
#include <iostream>

#include "GeoMoments.h"
#include "util.h"

using std::vector;
using std::complex;


template<typename T>
class ZernikeMoments
{
public:
	typedef vector<T> T1D;
	typedef vector<T1D> T2D;
	typedef vector<T2D> T3D;

	typedef std::complex<T> ComplexT;
	typedef NDVector<ComplexT, 3> ComplexT3D;
	
public:
	ZernikeMoments();
	ZernikeMoments(int _order, GeoMoments<T>& _gm);

	ComplexT getMoment(int n, int l, int m);

private:
	void computeCs(GeoMoments<T>& gm);
	void computeQs(GeoMoments<T>& gm);
	void computeMoments(GeoMoments<T>& gm);
	ComplexT computeMoment(GeoMoments<T>& gm, int n, int li, int m);
	
	ComplexT3D moments;
	T3D qs;
	T2D cs;

	int order;
};

#include "ZernikeMoments.hpp"

#endif


