
template<typename T>
GeoMoments<T>::GeoMoments(){};


template<typename T>
GeoMoments<T>::GeoMoments(Grid<T> &grid, const Vec3<T> &center, T scale, int n) : 
	samples{
		T1D(grid.sizeX() + 1),
		T1D(grid.sizeY() + 1),
		T1D(grid.sizeZ() + 1)
	},
	maxOrder(n)
{
	moments.resize(n+1);
	for(int i = 0; i <= n; i++){
		moments[i].resize(n - i + 1);
		for(int j = 0; j <= n - i; j++){
			moments[i][j].resize(n - i - j + 1);
		}
	}
	//printf("#center= %f %f %f scale= %f\n",center.x, center.y, center.z,scale);
	computeSamples(grid, center, scale);
	if(VERBOSE) cerr << "[info] Compute Samples." << endl;

	computeMoments(grid);

	if(isNaN(moments)){
		throw std::runtime_error("Computed a NaN!");
	}
}

template<typename T>
void GeoMoments<T>::computeSamples(const Grid<T> &grid, const Vec3<T> &center, T scale)
{
	for(int x = 0; x <= grid.sizeX(); x++) samples.x[x] = scale * (x - center.x);
	for(int y = 0; y <= grid.sizeY(); y++) samples.y[y] = scale * (y - center.y);
	for(int z = 0; z <= grid.sizeZ(); z++) samples.z[z] = scale * (z - center.z);
}

template<typename T>
void GeoMoments<T>::computeMoments(Grid<T> &grid)
{
	int xdim = grid.sizeX();
	int ydim = grid.sizeY();
	int zdim = grid.sizeZ();

	if(VERBOSE) cerr << "[info] Compute Moments." << endl;

	if(maxOrder == 1){
	 	if(VERBOSE) cerr << "[info] First time." << endl;
		computeFirst(grid);
		return;
	}

	int arrayDim = zdim;
	int layerDim = ydim * zdim;

	T1D diffLayer((ydim + 1) * zdim);
	T1D diffArray(zdim + 1);

	T1D layer(layerDim);
	T1D array(arrayDim);
	T   moment;

	T1D diffGridExtra(layerDim);

	auto iter = grid.begin();

	// repeated application of the F-recurrence in Novotni's paper
	#pragma omp parallel for
	for(int x=0; x<layerDim; ++x)
	{
		auto it = iter + xdim*x;
		computeDiffsInitial(it, it, xdim, diffGridExtra[x]);
		//printf("DiffX %d = %f\n",x,diffGridExtra[x]);
	}

	for(int i=0; i<=maxOrder; ++i)
	{
		auto diffIter = grid.begin();

		#pragma omp parallel for
		for(int p=0; p<layerDim; ++p)
		{
			auto diffIt = diffIter + (xdim)*p;
			auto sampleIter = samples.x.begin();
			layer[p] = multiplyInitial(diffIt, sampleIter, xdim, diffGridExtra[p]);
		}              

		iter = layer.begin();
		diffIter = diffLayer.begin();

		#pragma omp parallel for
		for(int y=0; y<arrayDim; ++y)
		{
			auto it = iter + ydim*y;
			auto diffIt = diffIter + (ydim+1)*y;
			computeDiffs(it, diffIt, ydim);
		}

		for(int j=0; j<maxOrder+1-i; ++j)
		{
			diffIter = diffLayer.begin();

			#pragma omp parallel for
			for(int p=0; p<arrayDim; ++p)
			{
				auto diffIt = diffIter + (ydim+1)*p;
				auto sampleIter = samples.y.begin();
				array[p] = multiply(diffIt, sampleIter, ydim + 1);
			}

			iter = array.begin();
			diffIter = diffArray.begin();
			computeDiffs(iter, diffIter, zdim);

			for(int k=0; k<maxOrder+1-i-j; ++k)
			{
				moment = multiply(diffIter, samples.z.begin(), zdim + 1);
				moments[i][j][k] = moment / ((1+i) * (1+j) * (1+k));
				//if(i==0 && j==0 && k==0)printf("moment %d %d %d %f\n",i,j,k,moment);
			}
		}
	}
}

template<typename T>
void GeoMoments<T>::computeFirst(Grid<T> &grid)
{
	int xdim = grid.sizeX();
	int ydim = grid.sizeY();
	int zdim = grid.sizeZ();

	double sss=0;
	int cnt=0;
	// here scale must be 1, so samples are indices
	Sum::SumStream<T> mass, cx,cy,cz;
	auto vinit = grid.begin();
	auto layerDim = ydim*zdim;
	auto v = vinit;
	for(int z = 0; z < zdim; z++){
		//auto v = vinit + z*layerDim; <--Bug??????
		for(int y = 0; y < ydim; y++){
			for(int x = 0; x < xdim; x++){
				mass << *v;
				cx << (*v * (2*x+1)/2.0);
				cy << (*v * (2*y+1)/2.0);
				cz << (*v * (2*z+1)/2.0);
				//if(*v == 1000.00){
				 //printf("%d %d %d %f\n",x,y,z,*v);
				// sss+=*v;
				// cnt++;
				//}
				v++;
			}
		}
	}
	//printf("Sum= %f %f %d\n",mass.getSum(),sss,cnt);
	moments[0][0][0] = mass.getSum();
	moments[1][0][0] = cx.getSum();
	moments[0][1][0] = cy.getSum();
	moments[0][0][1] = cz.getSum();

	if(isNaN(moments)){
		throw std::runtime_error("Computed a NaN!");
	}
}


template<typename T>
void GeoMoments<T>::computeDiffs(T1DIter iter, T1DIter diffIter, int dim)
{
	diffIter[0] = -iter[0];
	for(int i=1; i<dim; ++i)
	{
		diffIter[i] = iter[i-1] - iter[i];
	}
	diffIter[dim] = iter[dim-1];
}


template<typename T>
T GeoMoments<T>::multiply(T1DIter diffIter, T1DIter sampleIter, int dim)
{
	Sum::SumStream<T> sum;
	for(int i=0; i<dim; i++)
	{
		diffIter[i] *= sampleIter[i];
		sum << diffIter[i];
	}

	return sum.getSum();
}

template<typename T>
T GeoMoments<T>::getMoment(int i, int j, int k)
{
	return moments[i][j][k];
}

template<typename T>
void GeoMoments<T>::computeDiffsInitial(T1DIter iter, T1DIter diffIter, int dim, T &extra)
{
	T temp = (T)iter[0];
	T u = (T)iter[0];

	diffIter[0] = -temp;
	for(int i=1; i<dim; ++i)
	{
		u = (T)iter[i];
		diffIter[i] = temp - u;
	 	//printf("%d diff %f = tmp %f - iter %f\n",i,diffIter[i],temp,(T)iter[i]);
		temp = u;
	}
	extra = temp;
}

template<typename T>
T GeoMoments<T>::multiplyInitial(T1DIter diffIter, T1DIter sampleIter, int dim, T &extra)
{
	Sum::SumStream<T> sum;
	for(int i=0; i<dim; ++i)
	{
		diffIter[i] *= sampleIter[i];
		sum << diffIter[i];
	}

	extra *= sampleIter[dim];
	sum << extra;

	return sum.getSum();
}
