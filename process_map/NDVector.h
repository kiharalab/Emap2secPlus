#include <vector>
#include <initializer_list>

#ifndef NDVector_h

// NDVector implements a nested vector type.

template<typename Type, int n>
class NDVector : public std::vector<NDVector<Type, n-1> > {
	public:
		typedef std::initializer_list<typename NDVector<Type, n-1>::NDListType> NDListType;
};

template<typename Type>
class NDVector<Type, 1> : public std::vector<Type> {
	public:
		typedef std::initializer_list<Type> NDListType;
};

#define NDVector_h
#endif
