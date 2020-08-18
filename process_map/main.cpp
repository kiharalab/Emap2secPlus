#include <vector>
#include <sstream>
#include <iostream>
#include <istream>
#include <fstream>
#include <stdexcept>
#include <cassert>
#include <string>

#include "gzstream/gzstream.h"
#include "Grid.h"
#include "ZernikeDescriptor.h"
#include "util.h"
#include "pdb.h"

using namespace std;

bool VERBOSE = false;

static string formatFilenamePrefix(string s, string c){
	if(s.find("%s") == std::string::npos){
		return s + c;
	}else{
		char buf[1024];
		sprintf(buf, s.c_str(), c.c_str());
		return string(buf);
	}
}

static bool endsWith(const string &a, const string &b){
	return a.size() >= b.size() && a.compare(a.size()-b.size(), b.size(), b) == 0;
}

template<class T>
static void saveInvariants(ZernikeDescriptor<T> &zd, const string &path)
{
	std::ofstream outfile;
	outfile.exceptions(ios::badbit|ios::failbit);
	try{
		outfile.open(path.c_str());
	}catch(std::exception &e){
		throw std::runtime_error("Output file " + path + " could not be opened for writing.");
	}

	try{
		outfile << zd.size() << std::endl;
		for(auto &inv : zd){
			outfile << (inv/10) << std::endl;
		}
	}catch(std::exception &e){
		throw std::runtime_error("Write to output file " + path + " failed.");
	}
}

template<class T>
static void printInvariants(ZernikeDescriptor<T> &zd )
{

		//cout << zd.size() << std::endl;
	//printf("#%d",zd.size());
	int od=1;
		for(auto &inv : zd){
			printf(",%f",(inv/10));
			//printf("%d %f\n",od,(inv/10));
			//std::cout << (inv/10) << std::endl;
		 od++;
		}
	printf("\n");
	//printf("\n\n");
}
//Old format
template<class T>
static void printInvariantsV(ZernikeDescriptor<T> &zd )
{

		//cout << zd.size() << std::endl;
	printf("%d\n",zd.size());
	int od=1;
		for(auto &inv : zd){
			//printf(",%f",(inv/10));
			printf("%f\n",od,(inv/10));
			//std::cout << (inv/10) << std::endl;
		 od++;
		}
	//printf("\n");
	//printf("\n\n");
}



template<class T>
static double distInvariants(ZernikeDescriptor<T> &z1 ,ZernikeDescriptor<T> &z2 )
{
 double sum=0;
 auto Z1=z1.begin();
 auto Z2=z2.begin();
 for(int i=0;i<z1.size() && i<z2.size();i++ ){
  sum+=(Z1[i]-Z2[i])*(Z1[i]-Z2[i])*0.01;
  //printf("%f %f\n",Z1[i],Z2[i]);
 }
 return sqrt(sum);
}


// http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html


const double H_RAD = 1.20; // Hydrogen
const double C_RAD = 1.70; // carbon
const double N_RAD = 1.65; // nitrogen
const double O_RAD = 1.60; // oxygen
const double P_RAD = 1.90; // phosphorous
const double S_RAD = 1.80; // sulphur
const double DEF_RAD = 1.70; // default
const char kBlankChars[] = " \t\n\r";

/// Returns a string with leading/trailing characters of a set stripped
string trimmed(string const& str, char const* sepSet=kBlankChars)
{
    string::size_type const first = str.find_first_not_of(sepSet);
    return (first==string::npos )? string():str.substr(first, str.find_last_not_of(sepSet)-first+1);
}

string rtrimmed(string const& str, char const* sepSet)
{
    string::size_type const last = str.find_last_not_of(sepSet);
    return (last==string::npos)? string():str.substr(0, last+1);
}

string ltrimmed(string const& str, char const* sepSet)
{
    string::size_type const first = str.find_first_not_of(sepSet);
    return (first==string::npos)? string():str.substr(first);
}


void read_protein(string infile, vector<atom>& P)
{
    ifstream ifs;
    ifs.open(infile.c_str());
    if(!ifs)
    {
        cerr << "Cannot open file: " << infile  << endl;
        exit(EXIT_FAILURE);
    }

    string line;

    int atom_id, res_id;
    string res_name, atom_type;
    double d_xyz[3];
    string chain_id;
        double rad = 0.;
    int count_pdb,count_rna;
    count_pdb=0;
    count_rna=0;
    int count_model=0;
    while(getline(ifs, line))
    {
        if((line.substr(0, 5)).compare("MODEL") == 0){
            count_model+=1;
            if(count_model>1){
               cerr << "#PDB Model End dtected"<< endl;
                break;

            }
        }
        if((line.substr(0, 4)).compare("ATOM") == 0)
        {
            atom_id = atoi((line.substr(6, 5)).c_str());
            atom_type = trimmed(line.substr(12, 4), kBlankChars); // hydrogen atoms have 4 letter identification

            if(atom_type == "OT1" || atom_type == "OT2")
            {
                atom_type = "OXT";
            }

            if(atom_type.at(0) == 'H') rad = H_RAD;
            else if(atom_type.at(0) == 'C') rad = C_RAD;
            else if(atom_type.at(0) == 'S') rad = S_RAD;
                else if(atom_type.at(0) == 'P') rad = P_RAD;
                else if(atom_type.at(0) == 'O') rad = O_RAD;
                else if(atom_type.at(0) == 'N') rad = N_RAD;
                else rad = DEF_RAD;

            res_name = trimmed(line.substr(17, 3), kBlankChars);


            chain_id = trimmed(line.substr(21, 1), kBlankChars);
            res_id = atoi((line.substr(22, 4)).c_str());
            //printf("The value is %d\n",res_id);
            if((line.substr(17, 3)).compare("  A") == 0||(line.substr(17, 3)).compare("  U") == 0||
            (line.substr(17, 3)).compare("  T") == 0||(line.substr(17, 3)).compare("  C") == 0
            ||(line.substr(17, 3)).compare("  G") == 0||(line.substr(17, 3)).compare(" DG") == 0
            ||(line.substr(17, 3)).compare(" DC") == 0||(line.substr(17, 3)).compare(" DA") == 0
            ||(line.substr(17, 3)).compare(" DT") == 0||(line.substr(17, 3)).compare(" DU") == 0
            ||(line.substr(17, 3)).compare("A  ") == 0||(line.substr(17, 3)).compare("U  ") == 0||
            (line.substr(17, 3)).compare("T  ") == 0||(line.substr(17, 3)).compare("C  ") == 0
            ||(line.substr(17, 3)).compare("G  ") == 0||(line.substr(17, 3)).compare("DG ") == 0
            ||(line.substr(17, 3)).compare("DC ") == 0||(line.substr(17, 3)).compare("DA ") == 0
            ||(line.substr(17, 3)).compare("DT ") == 0||(line.substr(17, 3)).compare("DU ") == 0){
            atom_type="RNA";
            res_id=-999;//Usually the residue name is marked as
            //printf("Successful read the RNA -999");
            count_rna+=1;
            }else{
            count_pdb+=1;
            }
            //cout << chain_id << " " << res_id << endl;
            d_xyz[0] = atof((line.substr(30, 8)).c_str());
            d_xyz[1] = atof((line.substr(38, 8)).c_str());
            d_xyz[2] = atof((line.substr(46, 8)).c_str());

            // create atom object
            atom n_atom(atom_id, atom_type, res_name, chain_id, res_id, d_xyz, rad);
            P.push_back(n_atom);
        }
        if((line.substr(0, 5)).compare("ENDML") == 0){
            cerr << "#PDB Model End dtected"<< endl;
            break;
        }
       // printf ("Extracting pdb finished with  %d rna/dna and %d atoms in protein\n", count_rna, count_pdb);

//        if((line.substr(0, 6)).compare("HETATM") == 0){
//            atom_id=-999;//denotes it's a rna structure
//            atom_type="RNA";
//            res_name="HOH";
//            chain_id="XW";
//            res_id=-999;//mark as rna
//            d_xyz[0] = atof((line.substr(30, 8)).c_str());
//            d_xyz[1] = atof((line.substr(38, 8)).c_str());
//            d_xyz[2] = atof((line.substr(46, 8)).c_str());
//            rad=1.0;
//            atom n_atom(atom_id, atom_type, res_name, chain_id, res_id, d_xyz, rad);
//            P.push_back(n_atom);
//
//
//        }
    }// end while

 ifs.close();
 //cerr<<<<endl
  cerr << "protein atom count " << count_pdb<<" rna/dna count "<<count_rna<< endl;
}





int main(int argc, char** argv)
{
	assert(sizeof(char) == 1);
	assert(sizeof(int) == 4);
	assert(sizeof(float) == 4);
	assert(sizeof(double) == 8);
	assert(sizeof(long) == 8);

	Binomial<double>::computePascalsTriangle(60);

	string programName(argv[0]);
	int maxOrder = 20; // default to max stable order
	double contour = 0.0; //default to assuming preprocessed map
	string contourString = "1.0";
	string inFileName;
	string outFilePrefix;
	string pdbFile;
	bool isGzipped = false;
	bool dmode = false;
	bool imode = false;
	//Genki
	bool usePdb = false;
	double rad = 5.00;// default distance
	double MinDist = 1.0;
	int Vwindow=5;
	int vstep=1;
	int sstep=1;
	bool gnorm=true;
	bool igstart=false;
	vector<double> contours;
	vector<string> contourStrings;

	for(int i=1; i<argc; i++){
		string cur(argv[i]);

		if(cur == "-n" || cur == "--order"){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> maxOrder)){
					cerr << "[error] Expected integer in [1..20] after '" << cur << "', not '" << argv[i] << "'." << endl;
					exit(1);
				}
				if(maxOrder<1 || maxOrder>20){
					cerr << "[error] Expected integer in [1..20] after '" << cur << "', not '" << argv[i] << "'." << endl;
					exit(1);
				}
			}else{
					cerr << "[error] Expected integer in [1..20] after '" << cur << "'." << endl;
					exit(1);
			}
		}else if(cur == "-c" || cur == "--contour"){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> contour)){
					cerr << "[error] Expected real contour value after '" << cur << "'." << endl;
					exit(1);
				}
				contourString = argv[i];
			}else{
				cerr << "[error] Expected real contour value after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "--contours"){
			if((i+1) < argc){
				i++;
				stringstream ss,tt;
				ss << argv[i];
				tt << argv[i];
				double tmp;
				string stmp;
				while(ss >> tmp){
					contours.push_back(tmp);
					tt >> stmp;
					contourStrings.push_back(stmp);
				}
				if(ss.get() != EOF){
					cerr << "[error] Expected space-delimited list of contours after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected space-delimited list of contours after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-p" || cur == "--prefix"){
			if((i+1) < argc){
				i++;
				outFilePrefix = argv[i];
				if(outFilePrefix == ""){
					cerr << "[error] Expected a valid filename prefix after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected a valid filename prefix after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-r" || cur == "--radius"){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> rad)){
					cerr << "[error] Expected real radius value after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected real radius value after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-MinDist" ){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> MinDist)){
					cerr << "[error] Expected float after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected real radius value after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-vw" ){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> Vwindow )){
					cerr << "[error] Expected int after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected voxel window size value after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-vstep" ){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> vstep )){
					cerr << "[error] Expected int after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected voxel window size value after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-sstep" ){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> sstep )){
					cerr << "[error] Expected int after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected voxel window size value after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-h" || cur == "--help" || cur == "-?" || cur == "/?"){
			argc = -1; // print usage
			break;
		}else if(cur == "-g" || cur == "--gzip"){
			isGzipped = true;
		}else if(cur == "-v" || cur == "--verbose"){
			VERBOSE = true;
		}else if(cur == "-gnorm"){
			gnorm = true;
		}else if(cur == "-lnorm"){
			gnorm = false;
		}else if(cur == "-ignorestart"){
			igstart = true;
		}else if(cur == "-P" || cur == "--PDB"){
			if((i+1) < argc){
                                i++;
                                pdbFile = argv[i];
				usePdb = true;
                                if(pdbFile == ""){
                                        cerr << "[error] Expected a valid pdb file after '" << cur << "'." << endl;
                                        exit(1);
                                }
                        }else{
                                cerr << "[error] Expected a valid filename prefix after '" << cur << "'." << endl;
                                exit(1);
                        }
		}else{
			if(inFileName != ""){
				cerr << programName << ": ambiguous input file, is it '" << inFileName << "' or '" << cur << "'?" << endl;
				exit(1);
			}
			inFileName = cur;
		}
	}
	

	// if no args were given (or user asked for help)
	if(argc < 2){ 
		cerr << "EMmap -> training data ver 0.3" << endl << endl;

		cerr << "Usage: " << programName << " [OPTION]... FILE" << endl << endl;

		cerr << "FILE:" << endl;
		/*cerr << "  If FILE is -, the input map will be read from stdin." << endl;*/
		//cerr << "  Otherwise " << programName << " expects a valid filename. Supported file formats are" << endl;
		cerr << "  " << programName << " expects FILE to be a valid filename. Supported file formats are" << endl;
		cerr << "  Situs, CCP4, and MRC2000. Input may be gzipped. Format is deduced from FILE's" << endl;
		cerr << "  extension. Use the Situs tool map2map to inspect maps that fail to read." << endl << endl;

		cerr << "OPTION:" << endl;
		cerr << "  -c, --contour C      The level of isosurface to generate invariants for.def=0.0" << endl;
		cerr << "  -g, --gzip           Set this flag to force reading input as gzipped." << endl;
		cerr << "  -P PDBFILE           Use PDB file to use CA atom position." << endl;
		cerr << "  -r [float]           Max Distance from CA position.  def=5.0" << endl;
		cerr << "  -vw [int]            Size of voxel window def=5 (->11x11x11)" << endl;
		cerr << "  -sstep [int]         Step number [int] of the sliding voxel def=1" << endl;
		cerr << "  -vstep [int]         Step number [int] in the voxel def=1" << endl;
		cerr << "  -gnorm               Normalizef by Global Maximum density value def=true" << endl;
		cerr << "  -lnorm               Normalizef by Local Maximum density value def=false" << endl;
		cerr << "  -ignorestart         Ignore NSCTART... for simulated maps def=false" << endl;
		cerr << "  -h, --help, -?, /?   Displays this." << endl << endl;

		cerr << "EXAMPLES:" << endl;
		cerr << "  " << programName << " protein.situs -c 2.75" << endl;
		cerr << "  " << programName << " protein.map -P PDB.pdb -d 15.0" << endl << endl;

		return argc==-1 ? 0 : 1;
	}

	// default prefix
	if(outFilePrefix == ""){
		outFilePrefix = inFileName;
	}

	int type = TYPE_UNKNOWN;
	if(endsWith(inFileName, ".gz")){
		isGzipped = true;
		if(endsWith(inFileName, ".map.gz")||endsWith(inFileName, ".mrc.gz")){
			type = TYPE_MRC;
		}else if(endsWith(inFileName, ".situs.gz")){
			type = TYPE_SITUS;
		}
	}else{
		if(endsWith(inFileName, ".map")||endsWith(inFileName, ".mrc")){
			type = TYPE_MRC;
		}else if(endsWith(inFileName, ".situs")){
			type = TYPE_SITUS;
		}
	}

	if(type == TYPE_UNKNOWN){
		cerr << "[warning] Failed to deduce map type. Assuming default (Situs)." << endl;
		type = TYPE_SITUS;
	}


	istream *f = NULL;

	try{
		if(isGzipped){
			igzstream *ff = new igzstream();

			ff->exceptions(ios::badbit|ios::failbit);
			ff->open(inFileName.c_str(), ios::in|ios::binary);
			
			f = ff;
		}else{
			ifstream *ff = new ifstream();

			ff->exceptions(ios::badbit|ios::failbit);
			ff->open(inFileName.c_str(), ios::in|ios::binary);

			f = ff;
		}
	}catch(std::exception &e){
		cerr << "[error] Input file " << inFileName << " could not be opened (" << e.what() << ")." <<  endl;
		exit(1);
	}

	vector<atom> P;
	if(usePdb == true && pdbFile != ""){
	 if(VERBOSE) cerr << "[info] Reading pdb file." << endl;
    	 // read atoms
    	 read_protein(pdbFile, P);
	}

	try{
		Grid<double> g(*f, type);
		int Nact,Nact_pre;
		if(VERBOSE) cerr << "[info] Successfully read the input file." << endl;

		if(usePdb == true){
	 	 int cacnt=0;
		 Grid<double> pass(g,contour,P,rad,Vwindow,vstep,sstep,gnorm,igstart);
		}else{
		 Grid<double> pass(g,contour,rad,Vwindow,vstep,sstep,gnorm);
		}
	}catch(bad_alloc &e){
		cerr << "[error] Ran out of memory (" << e.what() << ")." << endl;
		exit(1);
	}catch(runtime_error &e){
		cerr << "[error] An error occurred (" << e.what() << ")." << endl;
		exit(1);
	}catch(exception &e){
		cerr << "[error] An unknown fatal error occurred (" << e.what() << ")." << endl;
		exit(1);
	}

	delete f;

	return 0;
}
