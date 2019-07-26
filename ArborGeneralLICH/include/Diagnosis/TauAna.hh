#ifndef _TauAna_hh_
#define _TauAna_hh_

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

class TauAna  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new TauAna ; }

		TauAna();

		~TauAna() {};

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
        float angle_ZDau, angle_HDau, angle_ZpH1, angle_ZpH2, angle_ZmH1, angle_ZmH2;

		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


