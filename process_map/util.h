
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