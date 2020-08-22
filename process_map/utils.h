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

#ifndef _UTILS_H_
#define _UTILS_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <set>
#include <map>
using namespace std;

#include "cp.h"
#include "pdb.h"
#include <ANN/ANN.h>

typedef double T3Matrix[3][3];
typedef double Transform[4][4];


const char kBlankChars[] = " \t\n\r";

void normalize(double*, int);
void point_transform (double*, Transform);
void normal_transform (double*, Transform);
double get_torsion_angle(double*, double*, double*, double*);

double compute_angle(double*, double*);
double compute_dotp(double*, double*);
double compute_angle(vector<double>&, vector<double>&);
double get_norm(double*, int);
void invert_normals(vector<cp>&);
bool check_cp(vector<cp>&);

double get_distance(double*, double*);
double get_distance2(double*, double*);

void filter_cp(vector<cp>&, vector<atom>&, vector<pair<string, string> >&);
string trimmed(string const&, char const*);
string to_string(int);
void apply_random_rotation(vector<atom>&, vector<cp>&);




#endif
