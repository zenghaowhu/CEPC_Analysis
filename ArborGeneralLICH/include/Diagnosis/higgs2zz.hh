#ifndef _higgs2zz_hh_
#define _higgs2zz_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

class TTree;

class higgs2zz  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new higgs2zz ; }

		higgs2zz();

		~higgs2zz() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;
//		TFile *tree_file;

		int _overwrite;
		TTree *_outputTree;

		unsigned int eventNr;
		int Num;
    int Z1Dau_PDG, Z2Dau_PDG;
    //for nnhllqq
    float mass_Hllqq;
    
    //for mmhnnqq
    int quark, lepton, neutrino;
    float angle_muon_pfo, angle_muon_mcp, mmMiss_mass, mmRecoil_mass;
    
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


