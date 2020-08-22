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
template<typename Value>
Sum::SumStream<Value>::SumStream(const SumStream<Value> &orig) : sum(orig.sum), correction(orig.correction), y(orig.y), t(orig.t) {

}

template<typename Value>
Sum::SumStream<Value>::SumStream() : sum(0), correction(0), y(0), t(0) {

}


template<typename Value>
void Sum::SumStream<Value>::operator<<(const Value &x){
	y = x - correction;
	t = sum+y;
	correction = (t-sum) - y;
	sum = t;
}

template<typename Value>
Value Sum::SumStream<Value>::getSum(){
	return sum;
}


template<typename Value>
Value Sum::kahanSum(std::vector<Value> vec){ // error only depends on Value's precision
	Value sum = 0;
	Value correction = 0;
	Value y, t;

	for(auto &v : vec){
		y = v - correction;
		t = sum + y;
		correction = (t - sum) - y;
		sum = t;
	}

	return sum;
}
