#ifndef _softLink2quark_hh_
#define _softLink2quark_hh_

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

class softLink2quark  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new softLink2quark ; }

		softLink2quark();

		~softLink2quark() {};

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

    std::vector<Float_t> vfinalPx;
    std::vector<Float_t> vfinalPy;
    std::vector<Float_t> vfinalPz;
    std::vector<Float_t> vfinalE;
    
    std::vector<Float_t> vquarkPx;
    std::vector<Float_t> vquarkPy;
    std::vector<Float_t> vquarkPz;
    std::vector<Float_t> vquarkE;

    int PDG_quark1;
    int PDG_quark2;

    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


