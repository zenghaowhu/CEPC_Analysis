#ifndef _softLink_hh_
#define _softLink_hh_

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

class softLink  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new softLink ; }

		softLink();

		~softLink() {};

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

    std::vector<Float_t> v941Px;
    std::vector<Float_t> v941Py;
    std::vector<Float_t> v941Pz;
    std::vector<Float_t> v941E;
    std::vector<Float_t> v942Px;
    std::vector<Float_t> v942Py;
    std::vector<Float_t> v942Pz;
    std::vector<Float_t> v942E;
    
    std::vector<Float_t> vquarknear941Px;
    std::vector<Float_t> vquarknear941Py;
    std::vector<Float_t> vquarknear941Pz;
    std::vector<Float_t> vquarknear941E;

    std::vector<Float_t> vquarknear942Px;
    std::vector<Float_t> vquarknear942Py;
    std::vector<Float_t> vquarknear942Pz;
    std::vector<Float_t> vquarknear942E;
    
    std::vector<Float_t> vquark941Px;
    std::vector<Float_t> vquark941Py;
    std::vector<Float_t> vquark941Pz;
    std::vector<Float_t> vquark941E;
    
    std::vector<Float_t> vquark942Px;
    std::vector<Float_t> vquark942Py;
    std::vector<Float_t> vquark942Pz;
    std::vector<Float_t> vquark942E;

    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


