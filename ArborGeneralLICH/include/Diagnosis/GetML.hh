#ifndef _GetML_hh_
#define _GetML_hh_

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

class GetML  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new GetML ; }

		GetML();

		~GetML() {};

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
    float Px_RecoJet1, Py_RecoJet1, Pz_RecoJet1, En_RecoJet1;
    float Px_RecoJet2, Py_RecoJet2, Pz_RecoJet2, En_RecoJet2;
    float Px_RecoJet3, Py_RecoJet3, Pz_RecoJet3, En_RecoJet3;
    float Px_RecoJet4, Py_RecoJet4, Pz_RecoJet4, En_RecoJet4;

    float Px_GenJet1, Py_GenJet1, Pz_GenJet1, En_GenJet1;
    float Px_GenJet2, Py_GenJet2, Pz_GenJet2, En_GenJet2;
    float Px_GenJet3, Py_GenJet3, Pz_GenJet3, En_GenJet3;
    float Px_GenJet4, Py_GenJet4, Pz_GenJet4, En_GenJet4;
    
    float Px_q1, Py_q1, Pz_q1, En_q1;
    float Px_q2, Py_q2, Pz_q2, En_q2;
    float Px_q3, Py_q3, Pz_q3, En_q3;
    float Px_q4, Py_q4, Pz_q4, En_q4;
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


