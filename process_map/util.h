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
#ifndef util_h
#define util_h

#include <vector>

#include "NDVector.h"

#define PI 3.141592653589793

template<typename T>
struct Vec3 {
	T x,y,z;
};

template<typename T>
bool isNaN(T x);

template<typename T>
bool isNaN(vector<T> vec);

template<typename T>
bool isNaN(vector<vector<T>> vec);

template<typename T>
bool isNaN(vector<vector<vector<T>>> vec);

template<class T>
class Binomial
{
public:
	typedef vector<T>           VectorT;
	typedef vector<VectorT>     VVectorT;

	static T get(int i, int j);
	static void computePascalsTriangle(int max);

private:
	static VVectorT triangle;
};      

template<class T>
typename Binomial<T>::VVectorT Binomial<T>::triangle;


#include "util.hpp"

#endif
