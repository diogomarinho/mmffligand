/*
 * main.cpp
 *
 *  Created on: 25/09/2011
 *      Author: diogo
 */


#include "MMFF94Parameters.h"
#include<openbabel/obconversion.h>
#include<openbabel/generic.h>
#include<openbabel/alias.h>
#include<openbabel/data.h>

#include<dirent.h>
#include<stdlib.h>
#include<iostream>

using namespace std;
using namespace OpenBabel;

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters );
void errorMessage();


void setupOptions(string lineOptions);
string vs_dir;
string ligand;
bool vs_option = false;
vector <string> vs;
bool hydrogens = false;
bool suppressOutfile = false;

int main(int argc, char **argv){

	string lineOptions;

	string aux;

	for(int i = 1; i < argc; i++){

		aux = argv[i];

		lineOptions += aux + ";";

	}

	if(argc > 1)
		setupOptions(lineOptions);

	else
		errorMessage();


	if (vs_option){

		DIR *dir;

		struct dirent *dirp;

		if((dir  = opendir(vs_dir.c_str())) == NULL) {

			cout << "Error in "<< vs_dir.c_str() << endl;

			exit(-1);

		}

		while ((dirp = readdir(dir)) != NULL){

				string filename = string(dirp->d_name);

				int lastIdx = int(filename.length()) - 1;

				string ext;

				ext.push_back(filename[lastIdx - 3]);
				ext.push_back(filename[lastIdx - 2]);
				ext.push_back(filename[lastIdx - 1]);
				ext.push_back(filename[lastIdx]);


				if (ext == ".mol" || ext == "mol2" || ext == ".pdb" || ext == ".sdf"){

					string dirName = vs_dir.c_str();

					int last = int(dirName.size()) - 1;

					if (dirName[last] !='/')
						dirName.push_back('/');

					vs.push_back(dirName+string(dirp->d_name));

				}

				ext.clear();
		    }

	}else{

		vs.push_back(ligand);

	}

	char id[] = "MMFF94s";
	MMFF94Parameters *p = new MMFF94Parameters(id, true);

	for (int i = 0; i < (int) vs.size(); i++){

		cout << "lingad: " << vs[i] << endl;

	//....................................................... Carregando molecula na memoria ...........................................................
		OBConversion conv;

		OBFormat *format_in = conv.FormatFromExt(vs[i].c_str());

		if (!format_in || !conv.SetInFormat(format_in) || !conv.SetOutFormat(format_in) ){

			cout << "Error in input file" << endl;

			exit(-1);
		}


		ifstream ifs(vs[i].c_str());

		OBMol mol;

		conv.ReadFile(&mol, vs[i]);

		if (mol.Empty()){

			cout << "Some error occur on load molecule, sorry!" << endl;

			exit(-1);

		}

		if (hydrogens){

			cout << "Add Hydrogens...\n";

			mol.AddHydrogens(false, true, 7.4);


			if (!suppressOutfile){

					string vi = "new_";

					string vi_file_name = vi + conv.GetTitle();

					conv.SetOutStream(new ofstream(vi_file_name.c_str()));

					conv.Write(&mol);

					conv.CloseOutFile();

				}
		}


	//...................................... Campo de força MMFF94s ..............................................................

//		cout << p->Description() << endl;

		p->Setup(mol);

		string ext = ".top";

		string out = mol.GetTitle(true) + ext;

		p->setOutFile(out);

		p->setSymbolsVector();

		p->setMMFFChargesList();

		p->setMMFFTorsionList();

		p->setMMFFElectrostaticList();

		p->setMMFFVdwList();

		p->setMMFFFBondsList();


	//.................................... gera o arquivo de topologia do ligante ................................................

		p->logBasicInformation();

		p->logAtomConection();

		p->logTorsionalParameters();

		p->logSelectedTorsions();

		p->setNonbondedInteractions();

	//....................................... FIM Do algoritmo e fechando o stream .................................................

		cout <<"The " + out + " file is generated.\nExiting." << endl;

		p->closeOutFile();

	}

//	delete p;

	return 0;
}


void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ") {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

void setupOptions(string lineOptions){

	vector<string>opt;

	Tokenize(lineOptions,opt,";");


	int s = (int) opt.size();

	bool condP = false;

	for (int i = 0; i < s; i++){

		if(opt[i]=="-l"){

			condP = true;

			ligand = opt[i+1];


		}if(opt[i]=="-vs"){

			condP = true;

			vs_option = true;

			vs_dir = opt[i+1];

		}else if (opt[i] == "-h"){

			hydrogens = true;

		}else if (opt[i] == "-s"){

			suppressOutfile = true;

		}
	}

	if (!(condP)){
		errorMessage();
	}

}

void errorMessage(){

	cerr << "Error in  command line.\n";

	cerr <<"You made a mistake, use: -l <ligand_file> or -vs <dir_with_ligands (pdb or mol2 or mol or sdf)>\n";

	cerr <<"Type [babel -L formats] to see the supported extensions\n";

	printf( "Use:\n"
			"mmffligand -vs dir      ..............  to generate multiple *.top files or\n"
			"mmffligand -vs dir -h   ..............  to generate multiple *.top files with hydrogens or\n"
			"mmffligand -l file.*    ..............  to generate the *.top file or\n"
			"mmffligand -l file.* -h ..............  to generate the *.top file adding hydrogens\n"
			"mmffligand -l file.* -h -s ...........  to supress the new molecule file with hydrogens\n"
			);

	exit(-1);
}

