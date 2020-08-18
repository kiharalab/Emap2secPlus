
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