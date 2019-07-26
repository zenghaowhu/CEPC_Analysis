#ifndef _nngamma_hh_
#define _nngamma_hh_

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

class nngamma  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new nngamma ; }

		nngamma();

		~nngamma() {};

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
        float select_j1j2, select_j3j4;
        int eventType;
        float costheta_1, costheta_2, costheta_3, costheta_4;
    
    float Gselect_j1j2, Gselect_j3j4;
    int GeventType;
    float Gcostheta_1, Gcostheta_2, Gcostheta_3, Gcostheta_4;
    
    float REn_barrel, REn_endcap, REn_overlap;
    float En_barrel, En_endcap, En_overlap;
    float REn_charge, REn_neutron;
    float En_charge, En_neutron;
    float En_ISR, Perp_ISR;
    float En_neutrino, Perp_Neutrino;
    float En_RecoJet1, En_RecoJet2, En_RecoJet3, En_RecoJet4;
    float En_GenJet1, En_GenJet2, En_GenJet3, En_GenJet4;
    float MCP_energy[100], MCP_costheta[100];
    float Reco_energy[100], Reco_costheta[100];
    int Rgamma_size;
    int gamma_size;
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


