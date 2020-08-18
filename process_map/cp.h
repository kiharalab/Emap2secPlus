#include <cmath>
#include <vector>
#include <string>
using namespace std;

#ifndef _CP_H_
#define _CP_H_

class cp
{
    public:
        int type;
        double xyz[3];
        vector<double> ZINV;
        double nrm[3];
        double cval;// curvature
        double mep; // mep value
        string res;
        cp(int t, double* d1, double* d2, string s, double f1, double f2)
        {
            type = t;
            for(int i=0;i<3;i++)
            {
                xyz[i] = d1[i];
                nrm[i] = d2[i];
            }
            res = s;
            cval = f1;
            mep = f2;
        }
        cp(){}
        ~cp(){}
};

#endif

