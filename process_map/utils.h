#ifndef _UTILS_H_
#define _UTILS_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <set>
#include <map>
using namespace std;

#include "cp.h"
#include "pdb.h"
#include <ANN/ANN.h>

typedef double T3Matrix[3][3];
typedef double Transform[4][4];


const char kBlankChars[] = " \t\n\r";

void normalize(double*, int);
void point_transform (double*, Transform);
void normal_transform (double*, Transform);
double get_torsion_angle(double*, double*, double*, double*);

double compute_angle(double*, double*);
double compute_dotp(double*, double*);
double compute_angle(vector<double>&, vector<double>&);
double get_norm(double*, int);
void invert_normals(vector<cp>&);
bool check_cp(vector<cp>&);

double get_distance(double*, double*);
double get_distance2(double*, double*);

void filter_cp(vector<cp>&, vector<atom>&, vector<pair<string, string> >&);
string trimmed(string const&, char const*);
string to_string(int);
void apply_random_rotation(vector<atom>&, vector<cp>&);




#endif
