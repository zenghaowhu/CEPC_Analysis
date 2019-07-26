#ifndef _ZYFastJet_hh_
#define _ZYFastJet_hh_

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

class ZYFastJet  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new ZYFastJet ; }

		ZYFastJet();

		~ZYFastJet() {};

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
    float B1_En, B2_En, B1_theta, B2_theta, B1_phi, B2_phi, B1_costheta, B2_costheta, B1_mass, B2_mass;

    
    std::vector<Float_t> vB_Px;

    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


