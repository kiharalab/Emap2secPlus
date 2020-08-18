/* prototypes for read.cc */
#ifndef _READ_H_
#define _READ_H_

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

#include "cp.h"
#include "pdb.h"


void read_protein(string, vector<atom>&);
void read_point_info(string, vector<cp>&);
void read_zd_info(string, vector<cp>&);
void read_interface_residues(string, vector<pair<string, string> >&);
void read_names(string, vector<string>&);

#endif


