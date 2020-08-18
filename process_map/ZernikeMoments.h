
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


