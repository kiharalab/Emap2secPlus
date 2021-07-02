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
template<typename T>
bool isNaN(T x){
	return x != x; //std library isnan only works for primitives
}

template<typename T>
bool isNaN(vector<T> vec){
	for(auto &x : vec){
		if(isNaN(x)){
			return true;
		}
	}
	return false;
}

template<typename T>
bool isNaN(vector<vector<T>> vec){
	for(auto &x : vec){
		if(isNaN(x)){
			return true;
		}
	}
	return false;
}

template<typename T>
bool isNaN(vector<vector<vector<T>>> vec){
	for(auto &x : vec){
		if(isNaN(x)){
			return true;
		}
	}
	return false;
}

template<class T>
void Binomial<T>::computePascalsTriangle(int max)
{
	triangle.resize(max+1); 

	for(int i = 0; i <= max; i++){
		triangle[i].resize(max+1-i);
		for(int j = 0; j <= max-i; j++){
			if(i==0 || j==0){
				triangle[i][j] = 1;
			}else{
				triangle[i][j] = triangle[i][j-1] + triangle[i-1][j];
			}
		}
	}
}

template<class T>
T Binomial<T>::get(int i, int j)
{
	return triangle[j][i-j];
}
