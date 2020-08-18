#include <vector>
#include <stdexcept>

#ifndef Sum_h

/*
	Tools for summation with O(1) error growth
*/

namespace Sum {

	template <typename Value=double>
	Value kahanSum(std::vector<Value> vec);

	template < typename Value=double >
	struct SumStream {
		public:
			SumStream(const SumStream<Value> &orig);
			SumStream();
			void operator<<(const Value &x);
			Value getSum();

		private:
			Value sum = 0;
			Value correction = 0;
			Value y, t;
	};


};

#include "Sum.hpp"

#define Sum_h
#endif