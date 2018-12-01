/*
 * MMFF94Parameters.h
 *
 *  Created on: 25/09/2011
 *      Author: diogo
 */

#ifndef MMFF94PARAMETERS_H_
#define MMFF94PARAMETERS_H_

#include"forcefieldmmff94.h"

using namespace std;


namespace OpenBabel {

class MMFF94Parameters : public  OBForceFieldMMFF94 {

private:
	vector <string> atom_symbol_list;
	vector <double> atom_charge_list;
	vector <OBFFTorsionCalculationMMFF94> tor;
	vector <OBFFVDWCalculationMMFF94> vdw;
	vector <OBFFElectrostaticCalculationMMFF94> ele;
	vector <OBFFBondCalculationMMFF94> bds;
	ofstream topFile;

public:

	MMFF94Parameters(const char* ID, bool IsDefault=true):OBForceFieldMMFF94(ID, IsDefault){

	}

	void setOutFile(string top_file_name);

	void closeOutFile();

	double getCharge(int i);

	void setSymbolsVector();

	void setMMFFChargesList();

	void setMMFFTorsionList();

	void setMMFFVdwList();

	void setMMFFElectrostaticList();

	void setMMFFFBondsList();

	void logBasicInformation();

	void logAtomConection();

	void logTorsionalParameters();

	void logSelectedTorsions();

	void setNonbondedInteractions();
};

} /* namespace OpenBabel */
#endif /* MMFF94PARAMETERS_H_ */

