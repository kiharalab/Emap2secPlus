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
#include <list>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include "pdb.h"

using std::cerr;
using std::endl;
using std::to_string;
using std::exception;
using std::runtime_error;

#define PI 3.14159265358979323846

static void reverse4(void *p){
	unsigned char *cptr, tmp;
	cptr = (unsigned char *)p;
	tmp = cptr[0];
	cptr[0]=cptr[3];
	cptr[3]=tmp;
	tmp = cptr[1];
	cptr[1]=cptr[2];
	cptr[2]=tmp;
}

static size_t permuted_index(int ordermode, size_t count, unsigned nc, unsigned nr, unsigned ns) {
	unsigned ic, ir, is;
	unsigned long ncr, q;

	ncr = nc*nr;
	is = count / ncr;
	q = count - is*ncr;
	ir = q / nc;
	ic = q - ir*nc;

	switch(ordermode) {
		case 1:
			return ic+ir*nc+is*nc*nr;
		case 2:
			return ic+is*nc+ir*nc*ns;
		case 3:
			return ir+ic*nr+is*nr*nc;
		case 4:
			return is+ic*ns+ir*ns*nc;
		case 5:
			return ir+is*nr+ic*nr*ns;
		case 6:
			return is+ir*ns+ic*ns*nr;
		default:
			throw runtime_error("Input file gives malformed dimension ordering.");
	}
}

static float read_float_istream(istream &fin, bool swap) {
	float currfloat;

	fin.read((char*)&currfloat,4);

	if(swap){
		reverse4(&currfloat);
	}
	return currfloat;
}

static float read_short_float_istream(istream &fin, bool swap) {
	unsigned char *cptr, tmp;
	short currshort;

	fin.read((char*)&currshort,2);
	if(swap){
		cptr = (unsigned char *)&currshort;
		tmp = cptr[0];
		cptr[0]=cptr[1];
		cptr[1]=tmp;
	}
	return (float)currshort;
}

static float read_char_float_istream(istream &fin) {
	char currchar;
	fin.read((char*)&currchar, 1);
	return (float)currchar;
}

static int read_int_istream(istream &fin, bool swap) {
	int currint;
	fin.read((char*)&currint, 4);
	if(swap){
		reverse4(&currint);
	}
	return currint;
}

template<typename Value>
void Grid<Value>::set(size_t x, size_t y, size_t z, Value val){
	if(x>=xdim || y>=ydim || z>=zdim){
		throw std::out_of_range("Attempted to access invalid grid position ("+to_string(x)+","+to_string(y)+","+to_string(z)+"), dimensions are ("+to_string(xdim)+","+to_string(ydim)+","+to_string(zdim)+").");
	}
	vec[(z*ydim+ y)*xdim + x] = val;
}
template<typename Value>
Value Grid<Value>::get(size_t x, size_t y, size_t z){
	if(x>=xdim || y>=ydim || z>=zdim){
		throw std::out_of_range("Attempted to access invalid grid position ("+to_string(x)+","+to_string(y)+","+to_string(z)+"), dimensions are ("+to_string(xdim)+","+to_string(ydim)+","+to_string(zdim)+").");
	}
	return vec[(z*ydim+ y)*xdim + x];
}
template<typename Value>
Value Grid<Value>::getOrDefault(size_t x, size_t y, size_t z, Value def){
	if(x>=xdim || y>=ydim || z>=zdim){
		return def;
	}
	return vec[(z*ydim+ y)*xdim + x];
}
template<typename Value>
Value Grid<Value>::getOrDefault(size_t x, size_t y, size_t z){
	if(x>=xdim || y>=ydim || z>=zdim){
		return defaultValue;
	}
	return vec[(z*ydim+ y)*xdim + x];
}

template<typename Value>
void Grid<Value>::readFrom(istream &fin, int type){
	try{
		if(type == TYPE_SITUS){
			readSitus(fin);
		}else{
			readMRC(fin);
		}
	}catch(runtime_error &e){
		cerr << "[error] Read from input file failed (" << e.what() << ")." << endl;
		exit(1);
	}

}

template<typename Value>
void Grid<Value>::readSitus(std::istream &fin){
	//double origx, origy, origz, width;
	double width;

	try{
		fin >> width >> origx >> origy >> origz >> xdim >> ydim >> zdim;
	}catch(exception &e){
		throw runtime_error("Situs file's header appears malformed.");
	}

	numVoxels = xdim * ydim * zdim;

	vec.reserve(numVoxels);
	double tmp;
	try{
		for(size_t i=0; i < numVoxels; i++){
			fin >> tmp;
			vec.push_back(tmp);
		}
	}catch(exception &e){
		throw runtime_error("Situs file's data block appears malformed.");
	}

	// if there is stuff after numVoxels reads, likely user error
	int didCatch = 0;
	try{
		double temp;
		fin >> temp;
	}catch(istream::failure &e){
		didCatch = 1;
	}

	if(!didCatch){
		throw runtime_error("Situs file's data block has more values than it should.");
	}

	stepx=width;
	stepy=width;
	stepz=width;

}

template<typename Value>
void Grid<Value>::interpSkewToCubic(double widthx, double widthy, double widthz, double alpha, double beta, double gamma) {
	Value *pphiout;
	unsigned long numVoxels;
	int pextx, pexty, pextz;
	double pwidth;
	double bx, by, cx, cy, cz, bix, biy, cix, ciy, ciz;
	double t1x, t1y, t2x, t2y, t3x, t3y, t4x, t4y, cdet, scz;
	double ux[8], uxmin, uxmax, uymin, uymax, uzmin, uzmax;
	double xpos, ypos, zpos, gx, gy, gz, a, b, c;
	int x0, y0, z0, x1, y1, z1;
	double endox, endoy, endoz;
	int i, indx, indy, indz;
	double porigx, porigy, porigz;

	bx = std::cos(PI*gamma/180.0); by = std::sin(PI*gamma/180.0);
	t1x = std::cos(PI*(gamma+alpha)/180.0); t1y = std::sin(PI*(gamma+alpha)/180.0);
	t2x = std::cos(PI*(gamma-alpha)/180.0); t2y = std::sin(PI*(gamma-alpha)/180.0);
	t3x = std::cos(PI*beta/180.0); t3y = std::sin(PI*beta/180.0);
	t4x = std::cos(PI*beta/180.0); t4y = -1.0 * std::sin(PI*beta/180.0);
	cdet = (t4y-t3y)*(t2x-t1x)-(t2y-t1y)*(t4x-t3x);

	cx = ((t4x-t3x)*(t1y*t2x-t2y*t1x)-(t2x-t1x)*(t3y*t4x-t4y*t3x))/cdet;
	cy = ((t4y-t3y)*(t1y*t2x-t2y*t1x)-(t2y-t1y)*(t3y*t4x-t4y*t3x))/cdet;
	scz = 1.0 - (cx*cx+cy*cy);

	cz = std::sqrt(scz);

	bix = -bx / by; biy = 1.0 / by;
	cix = (bx * cy - cx * by) / (by * cz); ciy = -cy / (by * cz); ciz = 1.0 / cz;

	endox = (xdim-1) * widthx;
	endoy = (ydim-1) * widthy;
	endoz = (zdim-1) * widthz;

	ux[1] = endox;
	ux[2] = bx * endoy;
	ux[3] = cx * endoz;
	ux[4] = endox + bx * endoy;
	ux[5] = bx * endoy + cx * endoz;
	ux[6] = endox + cx * endoz;
	ux[7] = endox + bx * endoy + cx * endoz;

	uxmin = 0; for(i=1;i<8;i++) if(ux[i] < uxmin) uxmin = ux[i];
	uxmax = 0; for(i=1;i<8;i++) if(ux[i] > uxmax) uxmax = ux[i];

	uymin = std::min(0.0, cy * endoz);
	uymin = std::min(uymin, by * endoy + cy * endoz);
	uymax = std::max(by * endoy, cy * endoz);
	uymax = std::max(uymin, by * endoy + cy * endoz);

	uzmin = std::min(0.0, cz * endoz);
	uzmax = std::max(0.0, cz * endoz);

	pwidth = widthx;
	if((widthy*by) < pwidth) pwidth = widthy*by;
	if((widthz*cz) < pwidth) pwidth = widthz*cz;
	pwidth = floor(pwidth*10.0+0.5)/10.0; // keeping with original behavior

	porigx = pwidth * floor(uxmin / pwidth);
	porigy = pwidth * floor(uymin / pwidth);
	porigz = pwidth * floor(uzmin / pwidth);

	pextx = (int) (ceil(uxmax/ pwidth) - floor(uxmin/ pwidth) + 1.5);
	pexty = (int) (ceil(uymax/ pwidth) - floor(uymin/ pwidth) + 1.5);
	pextz = (int) (ceil(uzmax/ pwidth) - floor(uzmin/ pwidth) + 1.5);
	numVoxels = pextx * pexty * pextz;

	pphiout = new Value[numVoxels];
	Value *p = pphiout;
	for(indz=0;indz<(int)pextz;indz++)
		for(indy=0;indy<(int)pexty;indy++)
			for(indx=0;indx<(int)pextx;indx++){

				xpos = porigx + indx * pwidth;
				ypos = porigy + indy * pwidth;
				zpos = porigz + indz * pwidth;

				gx = (xpos + bix * ypos + cix * zpos) / widthx;
				gy = (biy * ypos + ciy * zpos) / widthy;
				gz = (ypos + ciz * zpos) / widthz;
				x0 = (int)gx;
				y0 = (int)gy;
				z0 = (int)gz;
				x1 = x0+1;
				y1 = y0+1;
				z1 = z0+1;

				if(x0>=0 && x1<(int)xdim && y0>=0 && y1<(int)ydim && z0>=0 && z1<(int)zdim) {
					a = gx-x0;
					b = gy-y0;
					c = gz-z0;
					*p =
					a * b * c * get(x1, y1, z1) +
					(1-a) * b * c * get(x0, y1, z1) +
					a * (1-b) * c * get(x1, y0, z1) +
					a * b * (1-c) * get(x1, y1, z0) +
					a * (1-b) * (1-c) * get(x1, y0, z0) +
					(1-a) * b * (1-c) * get(x0, y1, z0) +
					(1-a) * (1-b) * c * get(x0, y0, z1) +
					(1-a) * (1-b) * (1-c) * get(x0, y0, z0);
				}
				p++;
			}

	xdim = pextx;
	ydim = pexty;
	zdim = pextz;
	vec.resize(numVoxels);
	Value *cur = pphiout;
	for(auto &v : vec){
		v = *(cur++);
	}

	delete pphiout;
}

template<typename Value>
void Grid<Value>::interpRectToCubic(double in_widthx, double in_widthy, double in_widthz) {
	double *outmap;
	double out_width = (int)std::max(in_widthx, std::max(in_widthy, in_widthz));
	out_width = floor(out_width*10.0+0.5)/10.0; // keeping with original behavior
	double xpos, ypos, zpos, gx, gy, gz, a, b, c;
	int x0, y0, z0, x1, y1, z1;

	xdim = ((int)((in_widthx*(xdim-1)))/out_width)+1;
	ydim = ((int)((in_widthy*(ydim-1)))/out_width)+1;
	zdim = ((int)((in_widthz*(zdim-1)))/out_width)+1;

	numVoxels = xdim * ydim * zdim;
	outmap = new Value[numVoxels];
	Value *p = outmap;

	for(int z=0;z<zdim;z++)
		for(int y=0;y<ydim;y++)
			for(int x=0;x<xdim;x++) {

				xpos = x * out_width;
				ypos = y * out_width;
				zpos = z * out_width;

				gx = (xpos/in_widthx);
				gy = (ypos/in_widthy);
				gz = (zpos/in_widthz);

				x0 = ((int)gx);
				y0 = ((int)gy);
				z0 = ((int)gz);
				x1 = ((int)gx) + 1;
				y1 = ((int)gy) + 1;
				z1 = ((int)gz) + 1;
				a = gx-x0;
				b = gy-y0;
				c = gz-z0;

				*p =
					a * b * c * get(x1, y1, z1) +
					(1-a) * b * c * get(x0, y1, z1) +
					a * (1-b) * c * get(x1, y0, z1) +
					a * b * (1-c) * get(x1, y1, z0) +
					a * (1-b) * (1-c) * get(x1, y0, z0) +
					(1-a) * b * (1-c) * get(x0, y1, z0) +
					(1-a) * (1-b) * c * get(x0, y0, z1) +
					(1-a) * (1-b) * (1-c) * get(x0, y0, z0);
				p++;
			}

	vec.resize(numVoxels);
	Value *cur = outmap;
	for(auto &v : vec){
		v = *(cur++);
	}

	delete outmap;
}

template<typename Value>
void Grid<Value>::readMRC(std::istream &fin){
	bool swap = false;

	int mx, my, mz;
	int mode;
	float xlen, ylen, zlen;
	int mapc, mapr, maps;
	float cfunsigned, prevfloat=0.0f, pfunsigned = 0.0f;
	double totdiffsigned=0.0, totdiffunsigned=0.0;
	float alpha, beta, gamma;
	int nsymbt;
	double widthx, widthy, widthz;
	int ncstart,nrstart,nsstart;

	xdim = read_int_istream(fin,swap);
	ydim = read_int_istream(fin,swap);
	zdim = read_int_istream(fin,swap);

	mode = read_int_istream(fin,swap);

	//fin.ignore(4*3);
	ncstart = read_int_istream(fin,swap);
	nrstart = read_int_istream(fin,swap);
	nsstart = read_int_istream(fin,swap);

	mx = read_int_istream(fin,swap);
	my = read_int_istream(fin,swap);
	mz = read_int_istream(fin,swap);

	xlen = read_float_istream(fin,swap);
	ylen = read_float_istream(fin,swap);
	zlen = read_float_istream(fin,swap);

	alpha = read_float_istream(fin,swap);
	beta = read_float_istream(fin,swap);
	gamma = read_float_istream(fin,swap);

	mapc = read_int_istream(fin,swap);
	mapr = read_int_istream(fin,swap);
	maps = read_int_istream(fin,swap);

	fin.ignore(4*4);//DMIN,DMAX,DMEAN,ISPG

	nsymbt = read_int_istream(fin,swap);

	fin.ignore(4*25);//EXTRA
	origx = read_float_istream(fin,swap);
	origy = read_float_istream(fin,swap);
	origz = read_float_istream(fin,swap);
	fin.ignore(4*1);//MAP
	//fin.ignore(4*29);
	char machst[4];
	fin.read(machst,4);
	if(machst[0] == 0x44){
		swap = false;
	}else if(machst[0] == 0x11){
		swap = true;
	}
	if(swap){
		if(VERBOSE) cerr << "[info] Map file is big-endian! Converting values to little-endian." << endl;

		reverse4(&xdim);
		reverse4(&ydim);
		reverse4(&zdim);

		reverse4(&mode);

		reverse4(&mx);
		reverse4(&my);
		reverse4(&mz);

		reverse4(&xlen);
		reverse4(&ylen);
		reverse4(&zlen);

		reverse4(&alpha);
		reverse4(&beta);
		reverse4(&gamma);

		reverse4(&mapc);
		reverse4(&mapr);
		reverse4(&maps);

		reverse4(&nsymbt);
	}


	int ordermode = 0;
	//int (*perm)(int);
	if(mapc==1 && mapr==2 && maps==3) {
		ordermode = 1;
		//perm = [&ic, &ir, &is, &x, &y, &z, &ncr, &q] (size_t c) { is = c/ncr; q = c-is*ncr; ir = q/x; ic = q - ir*x; return ic+ir*x+is*x*y; };
	}
	else if(mapc==1 && mapr==3 && maps==2) {
		ordermode = 2;
		//perm = [&ic, &ir, &is, &x, &y, &z, &ncr, &q] (size_t c) { is = c/ncr; q = c-is*ncr; ir = q/x; ic = q - ir*x; return ic+is*x+ir*x*z; };
	}
	else if(mapc==2 && mapr==1 && maps==3) {
		ordermode = 3;
		//perm = [&ic, &ir, &is, &x, &y, &z, &ncr, &q] (size_t c) { is = c/ncr; q = c-is*ncr; ir = q/x; ic = q - ir*x; return ir+ic*y+is*y*x; };
	}
	else if(mapc==2 && mapr==3 && maps==1) {
		ordermode = 4;
		//perm = [&ic, &ir, &is, &x, &y, &z, &ncr, &q] (size_t c) { is = c/ncr; q = c-is*ncr; ir = q/x; ic = q - ir*x; return is+ic*z+ir*z*x; };
	}
	else if(mapc==3 && mapr==1 && maps==2) {
		ordermode = 5;
		//perm = [&ic, &ir, &is, &x, &y, &z, &ncr, &q] (size_t c) { is = c/ncr; q = c-is*ncr; ir = q/x; ic = q - ir*x; return ir+is*y+ic*y*z; };
	}
	else if(mapc==3 && mapr==2 && maps==1) {
		ordermode = 6;
		//perm = [&ic, &ir, &is, &x, &y, &z, &ncr, &q] (size_t c) { is = c/ncr; q = c-is*ncr; ir = q/x; ic = q - ir*x; return is+ir*z+ic*z*y; };
	}
	else if(ordermode == 0) {
		throw runtime_error("Input file gives malformed dimension ordering.");
	}

	numVoxels = xdim*ydim*zdim;

	vec.resize(numVoxels); // initializes everything :C

	fin.ignore(4*(256-54)+nsymbt);

	if(VERBOSE) cerr << "[info] Map float mode is " << mode << "." << endl;

	switch(mode) {
		case 0: // char - converted to float, testing for signed-ness
			for(size_t i = 0; i<numVoxels; ++i) {
				size_t qt = permuted_index(ordermode,i,xdim,ydim,zdim);
				vec[qt] = read_char_float_istream(fin);
				if(vec[qt] < 0) cfunsigned = vec[qt] + 256.0;
				else cfunsigned = vec[qt];
				totdiffsigned += std::abs(vec[qt]-prevfloat);
				totdiffunsigned += std::abs(cfunsigned-pfunsigned);
				prevfloat=vec[qt];
				pfunsigned=cfunsigned;
			}
			if(totdiffsigned > totdiffunsigned) { // guess using unconventional unsigned
				for(size_t i = 0;i<numVoxels;i++) if (vec[i]<0) vec[i]+=256.0;
			}
			break;

		case 1: // 16-bit float
			for(size_t i = 0; i<numVoxels; ++i) {
				vec[permuted_index(ordermode,i,xdim,ydim,zdim)] = read_short_float_istream(fin,swap);
			}
			break;

		case 2: // 32-bit float
			for(size_t i = 0; i<numVoxels; ++i) {
				vec[permuted_index(ordermode,i,xdim,ydim,zdim)] = read_float_istream(fin,swap);
				//cout << "read in " << vec[permuted_index(ordermode,i,xdim,ydim,zdim)] << endl;
			}
			break;
		default:
			throw runtime_error("Unknown floating-point mode specified.");
	}

	widthx = xlen / (double) mx;
	widthy = ylen / (double) my;
	widthz = zlen / (double) mz;

	//printf("#width: %f %f %f\n",widthx,widthy,widthz);

	// is the grid orthogonal? cubic?
	if(alpha != 90.0 || beta != 90.0 || gamma != 90.0){
		cerr << "[warning] Input grid is not orthogonal; it will be interpolated to a cubic grid." << endl;
		interpSkewToCubic(widthx, widthy, widthz, alpha, beta, gamma);
		if(VERBOSE) cerr << "[info] Finished interpolating to cubic grid." << endl;
	}else if(std::abs(widthx-widthy)>1e-4 || std::abs(widthy-widthz)>1e-4 || std::abs(widthx-widthz)>1e-4){
		cerr << "[warning] Input grid is rectangular but not cubic; it will be interpolated to a cubic grid." << endl;
		interpRectToCubic(widthx, widthy, widthz);
		if(VERBOSE) cerr << "[info] Finished interpolating to cubic grid." << endl;
	}

	// make sure our dimension order matches what we loaded
	int nx,ny,nz;
	switch(ordermode){
		case 1:
			nx = xdim; ny = ydim; nz = zdim;
			origx2=ncstart*widthx+origx;
			origy2=nrstart*widthy+origy;
			origz2=nsstart*widthz+origz;
			break;
		case 2:
			nx = xdim; ny = zdim; nz = ydim;
			origx2=ncstart*widthx+origx;
			origy2=nsstart*widthy+origy;
			origz2=nrstart*widthz+origz;
			break;
		case 3:
			nx = ydim; ny = xdim; nz = zdim;
			origx2=nrstart*widthx+origx;
			origy2=ncstart*widthy+origy;
			origz2=nsstart*widthz+origz;
			break;
		case 4:
			nx = zdim; ny = xdim; nz = ydim;
			origx2=nsstart*widthx+origx;
			origy2=ncstart*widthy+origy;
			origz2=nrstart*widthz+origz;
			break;
		case 5:
			nx = ydim; ny = zdim; nz = xdim;
			origx2=nsstart*widthx+origx;
			origy2=nrstart*widthy+origy;
			origz2=nsstart*widthz+origz;
			break;
		case 6:
			nx = zdim; ny = ydim; nz = xdim;
			origx2=nsstart*widthx+origx;
			origy2=nrstart*widthy+origy;
			origz2=ncstart*widthz+origz;
			break;
		default:
			throw runtime_error("Input file gives malformed dimension ordering.");
	}


	xdim = nx;
	ydim = ny;
	zdim = nz;

	stepx=widthx;
	stepy=widthy;
	stepz=widthz;
}

template<typename Value>
void Grid<Value>::contour(Value c){
	if(VERBOSE) cerr << "[info] Check density value Here!." << endl;
	for(size_t i = 0; i < vec.size(); i++){
		if(vec[i] <= c){
			vec[i] = defaultValue;
		}else{
			vec[i] = indicatorValue;
		}
	}
}

template<typename Value>
Grid<Value>::Grid(istream &fin, int type){
	readFrom(fin, type);
}

template<typename Value>
Grid<Value>::Grid(istream &fin, int type, Value def){
	defaultValue = def;
	readFrom(fin, type);
}

template<typename Value>
Grid<Value>::Grid(Grid &g, Value contour){
	if(VERBOSE) cerr << "[info] Check density value Here!2." << endl;
	xdim = g.xdim;
	ydim = g.ydim;
	zdim = g.zdim;
	numVoxels = g.numVoxels;

	vec.reserve(numVoxels);
	ActV=0;
	for(auto &v : g.vec){
		if(v <= contour){
			vec.push_back(defaultValue);
		}else{
			vec.push_back(indicatorValue);
			//vec.push_back(v); //<--Interesting!!
			ActV++; //New Genki
		}
	}
}

//NEW by Genki
template<typename Value>
Grid<Value>::Grid(Grid &g, Value contour,bool mode){
	if(VERBOSE) cerr << "[info] Check density value Here Density Mode!!." << endl;
	xdim = g.xdim;
	ydim = g.ydim;
	zdim = g.zdim;
	numVoxels = g.numVoxels;

	vec.reserve(numVoxels);
	ActV=0;
/*
	//check max-min
	double dmax=contour;
	double dmin=10000000;
	for(auto &v : g.vec){
	 	if(v > contour){
		 if(dmax<v) dmax=v;
		 if(dmin>v) dmin=v;
		}
	}
	if(VERBOSE) cerr << "[info] dmin" << dmin<< "dmax." << dmax << endl;
	double base=indicatorValue/(dmax-dmin);
*/
	for(auto &v : g.vec){
		if(v > contour){
			vec.push_back(defaultValue);
		}else{
			vec.push_back(indicatorValue);
			ActV++; //New Genki
		}
	}


/*
	for(auto &v : g.vec){
		if(v <= contour){
			vec.push_back(defaultValue);
		}else{
			//vec.push_back(indicatorValue);
			vec.push_back((v-dmin)*base); //<--Interesting!!
			ActV++; //New Genki
		}
	}
*/
}


//use coordinates and radius restraint
template<typename Value>
//Grid<Value>::Grid(Grid &g, Value contour,double x, double y, double z,double r){
Grid<Value>::Grid(Grid &g, Value contour,double *cd,double r){
	if(VERBOSE) cerr << "[info] Check density value Here!3." << endl;
	xdim = g.xdim;
	ydim = g.ydim;
	zdim = g.zdim;
	int cnt=0;
	double r2=r*r;
	double rvx2=1/(g.stepx*g.stepx);
	double rvy2=1/(g.stepy*g.stepy);
	double rvz2=1/(g.stepz*g.stepz);
	numVoxels = g.numVoxels;

	vec.reserve(numVoxels);

	cnt=0;
	auto v = g.begin();
	for(int z = 0; z < zdim; z++){
	 double dz=(z-cd[2])*(z-cd[2])*rvz2;

	 for(int y = 0; y < ydim; y++){
	  double dyz=dz+(y-cd[1])*(y-cd[1])*rvy2;
	  for(int x = 0; x < xdim; x++){
	   double dxyz=dyz+(x-cd[0])*(x-cd[0])*rvx2;
	   //if(*v  <= contour ||dxyz >=r2 ){
	   if(*v  < contour ||dxyz >=r2 ){
		vec.push_back(defaultValue);
	   }else{
		vec.push_back(indicatorValue);
	    //printf("%d %d %d %f\n",x,y,z,*v);
	    cnt++;
	   }
           v++;
	  }
	 }
	}

	vector<Value> vec2(numVoxels,defaultValue);

	vector<int> chk(numVoxels,0);
	std::list<int> tbl;

	tbl.push_back(round(cd[0]));
	tbl.push_back(round(cd[1]));
	tbl.push_back(round(cd[2]));

	ActV=0;

	int tag=round(cd[0])+round(cd[1])*xdim+round(cd[2])*xdim*ydim;
	if(vec[tag]==indicatorValue){
	chk[tag]=1;
	vec2[tag]=indicatorValue;
	ActV++;
	while(tbl.size()>0){
	 int x=tbl.front();
	 tbl.pop_front();
	 int y=tbl.front();
	 tbl.pop_front();
	 int z=tbl.front();
	 tbl.pop_front();
	 //printf("Check %d %d %d %d\n",x,y,z,tbl.size());
	 for(int ix=std::max(x-1,0);ix<std::min(xdim,x+2);ix++){
	 for(int iy=std::max(y-1,0);iy<std::min(ydim,y+2);iy++){
	 for(int iz=std::max(z-1,0);iz<std::min(zdim,z+2);iz++){
	  tag=ix+iy*xdim+iz*xdim*ydim;
	  if(chk[tag]!=0) continue;
	  chk[tag]=1;
	  if(vec[tag]==indicatorValue){
	  tbl.push_back(ix);
	  tbl.push_back(iy);
	  tbl.push_back(iz);
	  vec2[tag]=indicatorValue;
	  ActV++;
	  }
	 }}}
	}

	}

	copy(vec2.begin(),vec2.end(),vec.begin());
	//printf("Step2=%f %f %f %d\n",g.stepx,g.stepy,g.stepz,cnt);

	//printf("Base=%f %f %f\n",g.origx,g.origy,g.origz);
}

//use coordinates and radius restraint
template<typename Value>
Grid<Value>::Grid(Grid &g, Value contour,vector<atom> &P,double r,int Vw,int vstep,int sstep,bool gnorm,bool igstart){
	if(VERBOSE) cerr << "[info] Check density value Here! PDB." << endl;
	xdim = g.xdim;
	ydim = g.ydim;
	zdim = g.zdim;
	int cnt=0;
	//double r2=r*r/((g.stepx*g.stepx)+(g.stepy*g.stepy)+(g.stepz*g.stepz));
	double r2=r*r;
	double rvx2=1/(g.stepx*g.stepx);
	double rvy2=1/(g.stepy*g.stepy);
	double rvz2=1/(g.stepz*g.stepz);
	numVoxels = g.numVoxels;

	int Nclose,ListClose[100];

	//igstart?
	if(igstart==false){
	 g.origx=g.origx2;
	 g.origy=g.origy2;
	 g.origz=g.origz2;
	}


	//printf("#Steps= %f %f %f\n",g.stepx,g.stepy,g.stepz);

	vec.reserve(numVoxels);

	cnt=0;
	//Filtering
	auto v = g.begin();
	for(int z = 0; z < zdim; z++){
	 for(int y = 0; y < ydim; y++){
	  for(int x = 0; x < xdim; x++){
	   if(*v  < contour){
		//vec.push_back(defaultValue);
		vec.push_back(0);
	   }else{
		vec.push_back(*v);
	    //printf("%d %d %d %f\n",x,y,z,*v);
	    cnt++;
	   }
           v++;
	  }
	 }
	}

	double dmax,dmin;
	dmax=0;dmin=1000000;
	for(int i=0;i<numVoxels;i++){
	 if(vec[i]>dmax)
	  dmax=vec[i];
	 if(vec[i]<dmin && vec[i] != 0)
	  dmin=vec[i];
	}
	int useful_voxel=0;
	int zero_voxel=0;
	for(int i=0;i<numVoxels;i++){
	 if(vec[i]>dmin)
	  useful_voxel+=1;
	  if(vec[i]==0)
	  zero_voxel+=1;
	}

	cerr<<"we should have voxles as output:"<<useful_voxel<<endl;
	cerr<<"we have voxles removed because of contour level:"<<zero_voxel<<endl;
	cerr<<"map size: "<<zdim<<" "<<ydim<<" "<<xdim<<endl;
	printf("#dmax= %f dmin= %f\n",dmax,dmin);
	printf("#Voxel Size %d %d %d\n",(xdim-1)/sstep+1,(ydim-1)/sstep+1,(zdim-1)/sstep+1);
	printf("#Base= %f %f %f\n",g.origx,g.origy,g.origz);
	useful_voxel=0;
	//center
	for(int z = 0; z < zdim; z+=sstep){
	 for(int y = 0; y < ydim; y+=sstep){
	  for(int x = 0; x < xdim; x+=sstep){
	   int tag=xdim*ydim*z+xdim*y+x;

	   if(vec[tag]==0&&igstart)
	      { continue;}
       else
           {useful_voxel+=1;}


	   //search closest CA atom

		double cxyz[3];
		int cacnt=0;
		double dismin=r2;
		int closestid=-2;
		int tmpres_id=0;
		double d2;
		Nclose=0;
		int mark_rna_use=0;
		for(int pos=0;pos<P.size();pos++){
		 if(P[pos].atype == "CA"||P[pos].atype == "N"||P[pos].atype == "C"||P[pos].atype == "O"||P[pos].atype=="RNA"){
		  //real space based coordinates

		  cxyz[0]=P[pos].axyz[0]-g.origx;
		  cxyz[1]=P[pos].axyz[1]-g.origy;
		  cxyz[2]=P[pos].axyz[2]-g.origz;
/*
		  cxyz[0]=(P[pos].axyz[0]-g.origx)/g.stepx;
		  cxyz[1]=(P[pos].axyz[1]-g.origy)/g.stepy;
		  cxyz[2]=(P[pos].axyz[2]-g.origz)/g.stepz;
*/
/*
		  d2=(cxyz[0]-(double)x)*(cxyz[0]-(double)x)*rvx2+
		     (cxyz[1]-(double)y)*(cxyz[1]-(double)y)*rvy2+
		     (cxyz[2]-(double)z)*(cxyz[2]-(double)z)*rvz2;
*/

		  d2=(cxyz[0]-(double)x*g.stepx)*(cxyz[0]-(double)x*g.stepx)+
		     (cxyz[1]-(double)y*g.stepy)*(cxyz[1]-(double)y*g.stepy)+
		     (cxyz[2]-(double)z*g.stepz)*(cxyz[2]-(double)z*g.stepz);
		  if(d2<dismin){

		   dismin=d2;
		   //if residue id is -999, directly mark it as -999
           tmpres_id=P[pos].rnum;

           if(tmpres_id==-999){
           closestid=-1000;
           if (mark_rna_use==0){
           cerr<<"We have one hit for this voxel,resid:"<<tmpres_id<<endl;
           ListClose[Nclose]=-999;
		   Nclose++;
		   mark_rna_use=1;
		   }


           }
           else{
           cerr<<"We have one hit for this voxel,resid:"<<tmpres_id<<" cacnt:"<<cacnt<<endl;
		   closestid=cacnt;
		   ListClose[Nclose]=cacnt;
		   Nclose++;
		   }

		   //printf("#CA= %f %f %f\n",P[pos].axyz[0],P[pos].axyz[1],P[pos].axyz[2]);
		   //printf("#OnGridCA= %f %f %f\n",cxyz[0],cxyz[1],cxyz[2]);
		  }
           if(P[pos].atype == "CA"){
           cacnt++;
           }


		 }
		}



		//Outside of contour level
		if(vec[tag]==0){

		 if(igstart){
		 if(closestid!=-2){
		 cerr<<"unlucky go here because of contour level: "<<closestid<<endl;}}
//		 else{
//		 printf("#C: Res= %d d= %.3f Window= %d Step= %f",-2,-1.0,Vw*2+1,g.stepx*vstep);
//		 printf(" Base= %f %f %f ",(double)(x-Vw*vstep)*g.stepx+g.origx,(double)(y-Vw*vstep)*g.stepy+g.origy,(double)(z-Vw*vstep)*g.stepz+g.origz);
//		 printf(" pos= %d %d %d\n",x,y,z);
//		 printf("%d\n",-2);}
		 continue;
		}
		if(closestid!=-2){cerr<<"Finally We have one hit for this voxel,resid:"<<closestid<<endl;}

		printf("#C: Res= %d d= %.3f Window= %d Step= %f",closestid+1,sqrt(dismin),Vw*2+1,g.stepx*vstep);
		printf(" Base= %f %f %f ",(double)(x-Vw*vstep)*g.stepx+g.origx,(double)(y-Vw*vstep)*g.stepy+g.origy,(double)(z-Vw*vstep)*g.stepz+g.origz);
		printf(" pos= %d %d %d\n",x,y,z);
		printf("%d",closestid+1);
		for(int i=0;i<Nclose;i++){
		if ((closestid==-1000&&ListClose[i]==-999)){//already marked in closest id
		    continue;
		}
		 if(ListClose[i]!=closestid)
		  printf(";%d",ListClose[i]);
        }
		//local max and min
		if(gnorm==false){
		 dmax=0;dmin=1000000;
		for(int vz=z-Vw*vstep;vz<=z+Vw*vstep;vz+=vstep){
		for(int vy=y-Vw*vstep;vy<=y+Vw*vstep;vy+=vstep){
		for(int vx=x-Vw*vstep;vx<=x+Vw*vstep;vx+=vstep){
		 double den=0;
		 if(vz>=0 &&vy>=0 &&vx>=0 && vz<zdim && vy < ydim &&vx < xdim){
		  den=vec[xdim*ydim*vz+xdim*vy+vx];
		  if(den>dmax)
	 	   dmax=den;
	 	  if(den<dmin)
	 	   dmin=den;
		 }
		}}}
		}
		double dwin=dmax-dmin;
	   //Show dens info
		for(int vz=z-Vw*vstep;vz<=z+Vw*vstep;vz+=vstep){
		for(int vy=y-Vw*vstep;vy<=y+Vw*vstep;vy+=vstep){
		for(int vx=x-Vw*vstep;vx<=x+Vw*vstep;vx+=vstep){
		 double den=0;
		 if(vz>=0 &&vy>=0 &&vx>=0 && vz<zdim && vy < ydim &&vx < xdim){
		  den=vec[xdim*ydim*vz+xdim*vy+vx];
		 }
		 //printf("V: %3d %3d %3d %.3f %f %f\n",vx,vy,vz,den,dmax,dmin);
		 if(den<dmin)
		  printf(",%.3f",0.000);
		 else 
		  printf(",%.3f",(den-dmin)/dwin);
		}}}
		printf("\n");
	  }
	 }
	}
	cerr<<"final usefule voxel count:"<<useful_voxel<<endl;

}

//use coordinates and radius restraint
template<typename Value>
Grid<Value>::Grid(Grid &g, Value contour,double r,int Vw,int vstep,int sstep,bool gnorm){
	if(VERBOSE) cerr << "[info] Check density value Here! without PDB." << endl;
	xdim = g.xdim;
	ydim = g.ydim;
	zdim = g.zdim;
	int cnt=0;
	//double r2=r*r/((g.stepx*g.stepx)+(g.stepy*g.stepy)+(g.stepz*g.stepz));
	double r2=r*r;
	double rvx2=1/(g.stepx*g.stepx);
	double rvy2=1/(g.stepy*g.stepy);
	double rvz2=1/(g.stepz*g.stepz);
	numVoxels = g.numVoxels;

	vec.reserve(numVoxels);

	cnt=0;
	//Filtering
	auto v = g.begin();
	for(int z = 0; z < zdim; z++){
	 for(int y = 0; y < ydim; y++){
	  for(int x = 0; x < xdim; x++){
	   if(*v  < contour){
		//vec.push_back(defaultValue);
		vec.push_back(0);
	   }else{
		vec.push_back(*v);
	    //printf("%d %d %d %f\n",x,y,z,*v);
	    cnt++;
	   }
           v++;
	  }
	 }
	}

	double dmax,dmin;
	dmax=0;dmin=1000000;
	for(int i=0;i<numVoxels;i++){
	 if(vec[i]>dmax)
	  dmax=vec[i];
	 if(vec[i]<dmin && vec[i] != 0)
	  dmin=vec[i];
	}
	printf("#dmax= %f dmin= %f\n",dmax,dmin);
	//center
	for(int z = 0; z < zdim; z+=sstep){
	 for(int y = 0; y < ydim; y+=sstep){
	  for(int x = 0; x < xdim; x+=sstep){
	   int tag=xdim*ydim*z+xdim*y+x;
	   if(vec[tag]==0)
	    continue;


	   //search closest CA atom

	   	printf("#C: Res= %d d= %.3f Window= %d Step= %f",0,sqrt(dmin),Vw*2+1,g.stepx*vstep);
		printf(" Base= %f %f %f\n",(double)(x-Vw*vstep)*g.stepx+g.origx,(double)(y-Vw*vstep)*g.stepy+g.origy,(double)(z-Vw*vstep)*g.stepz+g.origz);

		printf("%d",0);
		//local max and min
		if(gnorm==false){
		 dmax=0;dmin=1000000;
		for(int vz=z-Vw*vstep;vz<=z+Vw*vstep;vz+=vstep){
		for(int vy=y-Vw*vstep;vy<=y+Vw*vstep;vy+=vstep){
		for(int vx=x-Vw*vstep;vx<=x+Vw*vstep;vx+=vstep){
		 double den=0;
		 if(vz>=0 &&vy>=0 &&vx>=0 && vz<zdim && vy < ydim &&vx < xdim){
		  den=vec[xdim*ydim*vz+xdim*vy+vx];
		  if(den>dmax)
	 	   dmax=den;
	 	  if(den<dmin)
	 	   dmin=den;
		 }
		}}}
		}
		double dwin=dmax-dmin;
	   //Show dens info
		for(int vz=z-Vw*vstep;vz<=z+Vw*vstep;vz+=vstep){
		for(int vy=y-Vw*vstep;vy<=y+Vw*vstep;vy+=vstep){
		for(int vx=x-Vw*vstep;vx<=x+Vw*vstep;vx+=vstep){
		 double den=0;
		 if(vz>=0 &&vy>=0 &&vx>=0 && vz<zdim && vy < ydim &&vx < xdim){
		  den=vec[xdim*ydim*vz+xdim*vy+vx];
		 }
		 //printf("V: %3d %3d %3d %.3f %f %f\n",vx,vy,vz,den,dmax,dmin);
		 if(den<dmin)
		  printf(",%.3f",0.000);
		 else 
		  printf(",%.3f",(den-dmin)/dwin);
		}}}
		printf("\n");
	  }
	 }
	}

}



template<typename Value>
double Grid<Value>::pos_dens(Grid &g,double *cd) {
	int tag=round(cd[0])+round(cd[1])*xdim+round(cd[2])*xdim*ydim;
	return g[tag];
}


template<typename Value>
size_t Grid<Value>::sizeX() const{
	return xdim;
}

template<typename Value>
size_t Grid<Value>::sizeY() const{
	return ydim;
}

template<typename Value>
size_t Grid<Value>::sizeZ() const{
	return zdim;
}

template<typename Value>
size_t Grid<Value>::size() const{
	return numVoxels;
}

template<typename Value>
size_t Grid<Value>::sizeMax() const{
	return max(sizeX(),max(sizeY(), sizeZ()));
}


template<typename Value>
typename vector<Value>::iterator Grid<Value>::begin(){
	return vec.begin();
}

template<typename Value>
typename vector<Value>::iterator Grid<Value>::end(){
	return vec.end();
}


template<typename Value>
typename vector<Value>::iterator Grid<Value>::begin(size_t k, size_t n){
	return vec.begin()+vec.size()*k/n;
}

template<typename Value>
typename vector<Value>::iterator Grid<Value>::end(size_t k, size_t n){
	return begin(k+1,n);
}

#undef PI
