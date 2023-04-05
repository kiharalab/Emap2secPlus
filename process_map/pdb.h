/*
# Emap2sec+ is a computational tool using deep learning that can accurately identify structures, alpha helices, beta sheets, other(coils/turns) and DNA/RNA, in cryo-Electron Microscopy (EM) maps of medium to low resolution.
# Copyright (C) 2020 Xiao Wang, Eman Alnabati, Tunde W Aderinwale, Sai Raghavendra Maddhuri, Genki Terashi, Daisuke Kihara, and Purdue University.
# License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)
# Contact: Daisuke Kihara (dkihara@purdue.edu)


#

# This program is free software: you can redistribute it and/or modify

# it under the terms of the GNU General Public License as published by

# the Free Software Foundation, version 3.

#

# This program is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of

# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

# GNU General Public License V3 for more details.

#

# You should have received a copy of the GNU v3.0 General Public License

# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.en.html.
 * can be used interchangeably with ".pqr" files
 */

#include <string>
#include <vector>

using namespace std;

typedef double Transform[4][4];
typedef double T3Matrix[3][3];



#ifndef _ATOM_H_
#define _ATOM_H_

class atom
{
    public:
        int         anum;   // id of the atom
        string      atype;  // type of the atom
        int         rnum;   // residue id
        double      axyz[3]; // xyz coordinates
        string      residue; // name of residue
        string      chain;// chain id
	    double      arad;

        atom(int m_num, string m_type, string m_residue, string m_chain, int m_rid, double* m_xyz,
             double m_rad)
        {
            anum = m_num;
            rnum = m_rid;
            atype = m_type;
            for(int i=0;i<3;i++)
            {
                axyz[i] = m_xyz[i];
            }
            residue = m_residue;
            chain = m_chain;
			arad = m_rad;
        }

        atom(){}
        ~atom(){}
};

#endif



/**
 * PDB header file
 * stores atom coordinates,types and bond information
 *
 * @author vishwesh venkatraman
 * @date   19/08/2007
 */

#ifndef _PDB_H_
#define _PDB_H_

class pdb
{
public:
    string pname;
    vector<atom> atoms;


    pdb(string m_name, vector<atom>& m_atoms)
    {
        pname = m_name;
        atoms = m_atoms;
    }
    pdb(){}
    ~pdb(){};
};

#endif


#ifndef _NEIGHBOUR_H_
#define _NEIGHBOUR_H_

class Neighbour 
{
    public:
        double x;
        double y;
        double z;
        double dot;
};

#endif


