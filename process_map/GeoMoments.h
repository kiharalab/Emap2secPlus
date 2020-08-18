
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

