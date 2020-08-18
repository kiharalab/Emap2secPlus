/*
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


