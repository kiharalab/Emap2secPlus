

template<typename T>
ZernikeMoments<T>::ZernikeMoments(){}

template<typename T>
ZernikeMoments<T>::ZernikeMoments(int n, GeoMoments<T>& _gm)
{
	order = n;

	computeCs(_gm);
	computeQs(_gm);
	computeMoments(_gm);
}

template<typename T>
void ZernikeMoments<T>::computeCs(GeoMoments<T>& _gm)
{
	T tmp;
	cs.resize(order + 1);

	for(int l=0; l<=order; ++l){
		cs[l].resize(l + 1);

		tmp = 2*l+1;
		cs[l][0] = sqrt(tmp);

		for(int m=1; m<=l; ++m){
			tmp = tmp*(l+m)/(l-m+1);
			cs[l][m] = sqrt(tmp);
		}
	}
}

template<typename T>
void ZernikeMoments<T>::computeQs(GeoMoments<T>& _gm)
{
	qs.resize(order + 1);

	for(int n=0; n<=order; ++n)
	{
		qs[n].resize(n / 2 + 1);

		int l0 = n % 2;
		for(int l=l0; l<=n; l+=2)
		{
			int k = (n-l)/2;

			qs[n][l/2].resize(k + 1);

			for(int mu=0; mu<=k; ++mu)
			{
				T nom = Binomial<T>::get(2*k, k) *
						Binomial<T>::get(k, mu) *
						Binomial<T>::get(2 * (k + l + mu) + 1, 2 * k) *
						( 1 - 2*((k+mu) % 2) ); // negate if quantity is odd

				T den = (1L<<(2*k)) * Binomial<T>::get(k + l + mu, k);
				
				qs[n][l/2][mu] = nom / den * sqrt((T)(2 * l + 4 * k + 3) / 3.0);
			}
		}
	}
}

template<typename T>
void ZernikeMoments<T>::computeMoments(GeoMoments<T> &gm)
{
	moments.resize(order + 1);
	for(int n=0; n<=order; ++n){
		moments[n].resize(n/2 + 1);
		int li = 0, l0 = n%2;
		for(int l=l0; l<=n; ++li, l+=2){
			moments[n][li].resize(l + 1);
			for(int m=0; m<=l; ++m){
				moments[n][li][m] = computeMoment(gm, n, li, m);
			}
		}
	}
}

template<typename T>
typename ZernikeMoments<T>::ComplexT ZernikeMoments<T>::computeMoment(GeoMoments<T>& gm, int n, int li, int m){
	ComplexT zm((T)0, (T)0);

	int l = 2*li + n%2;

	T w = cs[l][m] / (1L<<m);

	int k = (n-l)/2;

	for(int nu=0; nu<=k; nu++){
		T w_Nu = w * qs[n][li][nu];
		for(int alpha=0; alpha<=nu; alpha++){
			T w_NuA = w_Nu * Binomial<T>::get(nu, alpha);
			for(int beta=0; beta<=nu-alpha; beta++){
				T w_NuAB = w_NuA * Binomial<T>::get(nu-alpha, beta);
				for(int p=0; p<=m; p++){
					T w_NuABP = w_NuAB * Binomial<T>::get(m, p);
					for(int mu = 0; mu <= (l-m)/2; mu++){
						T w_NuABPMu = w_NuABP * Binomial<T>::get(l, mu) * Binomial<T>::get(l-mu, m+mu) / (1L<<(2*mu));
						for(int q = 0; q <= mu; q++){
							T w_NuABPMuQ = w_NuABPMu * Binomial<T>::get(mu, q) * (1-2*((m-p+mu)%2));

							ComplexT c;
							switch(p % 4)
							{
								case 0: c = ComplexT(w_NuABPMuQ, (T)0); break;
								case 1: c = ComplexT((T)0, w_NuABPMuQ); break;
								case 2: c = ComplexT((T)(-1) * w_NuABPMuQ, (T)0); break;
								case 3: c = ComplexT((T)0, (T)(-1) * w_NuABPMuQ); break;
							}

							int z_i = l-m+2*(nu-alpha-beta-mu);
							int y_i = 2*(mu-q+beta)+m-p;
							int x_i = 2*q+p+2*alpha;

							zm = zm + (std::conj(c) * gm.getMoment(x_i, y_i, z_i));
						} // q
					} // mu
				} // p
			} // beta
		} // alpha
	} // nu

	return zm*((T)(3.0 / (4.0 * PI)));
}

template<typename T>
typename ZernikeMoments<T>::ComplexT ZernikeMoments<T>::getMoment(int _n, int _l, int _m)
{
	if(_m >= 0)
	{
		return moments[_n][_l/2][_m];
	}
	else
	{
		T sign;
		if(_m%2)
		{
			sign = (T)(-1);
		}
		else
		{
			sign = (T)1;
		}
		return sign * std::conj(moments[_n][_l/2][(int)abs((float)_m)]);
	}
}
