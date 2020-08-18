
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
