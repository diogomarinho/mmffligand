/*
 * MMFF94Parameters.cpp
 *
 *  Created on: 25/09/2011
 *      Author: diogo
 */

#include "MMFF94Parameters.h"
#include<iostream>
#include<stdio.h>

using namespace std;

namespace OpenBabel {

void MMFF94Parameters::setOutFile(string top_file_name){

	topFile.open(top_file_name.c_str());
}

double MMFF94Parameters::getCharge(int i){

	return _mol.GetAtom(i)->GetPartialCharge();

}


void MMFF94Parameters::setSymbolsVector(){

	OBElementTable table;

	FOR_ATOMS_OF_MOL(atom, _mol){

		string symbol_aux = table.GetSymbol(atom->GetAtomicNum());

		atom_symbol_list.push_back(symbol_aux);

	}

}

void MMFF94Parameters::setMMFFChargesList(){

	FOR_ATOMS_OF_MOL(atom, _mol){

		atom_charge_list.push_back(atom->GetPartialCharge());

	}

}

void MMFF94Parameters::closeOutFile(){

	topFile.close();

}

void MMFF94Parameters::setMMFFTorsionList(){

	tor = _torsioncalculations;

}


void MMFF94Parameters::setMMFFVdwList(){

	vdw = _vdwcalculations;

}

void MMFF94Parameters::setMMFFElectrostaticList(){

	ele = _electrostaticcalculations;

}

void MMFF94Parameters::setMMFFFBondsList(){

	bds = _bondcalculations;

}

void MMFF94Parameters::logBasicInformation(){

	char buffer[100];

	sprintf(buffer, "$TOTAL_ATOMS %3d\n", _mol.NumAtoms());

	topFile << buffer ;

	int i = 0;

	FOR_ATOMS_OF_MOL(atom, _mol){

		//             label  id  typ  chg     x       y        z     at_n
		sprintf(buffer, "%4s %5d %5s %12.5lf %12.5lf %12.5lf %12.5lf %5d\n",
		atom_symbol_list[i].c_str(), atom->GetIdx(), atom->GetType(), atom_charge_list[i] ,atom->GetX(), atom->GetY(), atom->GetZ(), atom->GetAtomicNum());

		topFile << buffer;

		i++;

	}
}


void MMFF94Parameters::logAtomConection(){

	char buffer[100];

	sprintf(buffer, "\n$CON %3d\n", _mol.NumBonds());

	topFile << buffer;


	FOR_BONDS_OF_MOL(bond, _mol){

		sprintf(buffer, "%5d %5d\n", bond->GetBeginAtomIdx(), bond->GetEndAtomIdx());

		topFile << buffer;

	}
}


void MMFF94Parameters::logTorsionalParameters(){

	char buffer[100];

	sprintf(buffer, "\n $TOR %5d\n", (int) tor.size());

	topFile << buffer;

	for (int idx = 0; idx < (int) tor.size(); idx++){

		OBAtom *i = tor[idx].a; OBAtom *j = tor[idx].b; OBAtom *k = tor[idx].c; OBAtom *l = tor[idx].d;

		sprintf(buffer, "%5d %5d %5d %5d %10.5lf %10.5lf %10.5lf\n",
				i->GetIdx(), j->GetIdx(), k->GetIdx(), l->GetIdx(), tor[idx].v1, tor[idx].v2, tor[idx].v3);

		topFile << buffer;

	}

}

void MMFF94Parameters::logSelectedTorsions(){

	vector<OBBond *> bondRors;

	int nrotors = 0;

	for(int i = 0; i <(int)_mol.NumBonds(); i++){

		OBBond *bond = _mol.GetBond(i);

		if(bond->IsRotor() && !bond->IsAmide()){


			bondRors.push_back(bond);

			nrotors++;

		}else{

			if (bond->IsSingle() && !bond->IsAmide()){

				OBAtom *a = bond->GetBeginAtom();

				OBAtom *b = bond->GetEndAtom();

				if ((a->IsCarbon() && b->IsHbondDonor() && b->IsOxygen()) || (b->IsCarbon() && a->IsHbondDonor() && a->IsOxygen()) ){

					bondRors.push_back(bond);

					nrotors++;

				}

			}


		}

	}


	char buffer[100];

	sprintf(buffer, "\n$SELECTED_TORSIONS %3d\n", nrotors);

	topFile << buffer;

	for (int i = 0 ; i < nrotors; i++){

		OBBond *bond= bondRors[i];

		sprintf(buffer, "%2d %5d %5d\n", i+1, bond->GetBeginAtomIdx(), bond->GetEndAtomIdx());

		topFile << buffer;

	}

}

void MMFF94Parameters::setNonbondedInteractions(){

	char buffer[100];

	topFile << "\n$NONBONDED_INTERACTIONS\n";

	for (int i = 0; i < (int) vdw.size(); i++){

		int a = vdw[i].a->GetIdx();

		int b = vdw[i].b->GetIdx();


		string scale = "1.00";

		if (vdw[i].a->IsOneFour(vdw[i].b))	scale = "0.75";

		sprintf(buffer, "%d\t %d\t  %s\n", a, b, scale.c_str());

		topFile << buffer;

	}

	sprintf(buffer, "\n$TOTAL_INTRAMOLECULAR_INTERACTIONS: %3d\n\n", (int)vdw.size());

	topFile << buffer;
}


} /* namespace OpenBabel */
