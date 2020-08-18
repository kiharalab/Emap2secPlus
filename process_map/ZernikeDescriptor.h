
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



