#ifndef _Ztoqq_hh_
#define _Ztoqq_hh_

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

class Ztoqq  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new Ztoqq ; }

		Ztoqq();

		~Ztoqq() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;

		int _overwrite;
		TTree *_outputTree;

		unsigned int eventNr;
		int Num;
    
        //about event shape
    double CParameter;
    double hemiMass1, hemiMass2;
    double hemiBroadening1, hemiBroadening2;
    double T;
    double angleTwoQ;
    double En_twoQ;
    double En_ISR;
    double Q1BackThrust;
    double Q1Thrust;
    int Q1PDG, Q2PDG;
    
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


