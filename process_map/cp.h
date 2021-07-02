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
#include <cmath>
#include <vector>
#include <string>
using namespace std;

#ifndef _CP_H_
#define _CP_H_

class cp
{
    public:
        int type;
        double xyz[3];
        vector<double> ZINV;
        double nrm[3];
        double cval;// curvature
        double mep; // mep value
        string res;
        cp(int t, double* d1, double* d2, string s, double f1, double f2)
        {
            type = t;
            for(int i=0;i<3;i++)
            {
                xyz[i] = d1[i];
                nrm[i] = d2[i];
            }
            res = s;
            cval = f1;
            mep = f2;
        }
        cp(){}
        ~cp(){}
};

#endif

