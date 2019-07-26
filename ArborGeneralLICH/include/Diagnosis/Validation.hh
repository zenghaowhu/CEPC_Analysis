#ifndef _Validation_hh_
#define _Validation_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class Validation  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new Validation ; }

		Validation();

		~Validation() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;

		int _overwrite, _leptonID;
		float _cmsE; 
		TTree *_outputTree, *_outputPFO;

//   for MCParticle
    int type;
    float N3En, N3Pt, ISRPt, MaxJetCosTheta;
    int num_higgsDau;
        int num_mcparticle;
        int num_gamma_mcp, num_neutral_hadron_mcp, num_charge_mcp;
        int lepton_plus_count_mcp, lepton_minus_count_mcp;
        float total_energy_mcp, mcp_left_mass;
        float mcp_lepton_minus_energy, mcp_lepton_plus_energy;
        float mcp_dilepton_mass;
    int lepton_plus_decayedfromHiggs_count_mcp, lepton_minus_decayedfromHiggs_count_mcp;
    float mcp_lepton_minus_decayedfromHiggs_energy, mcp_lepton_plus_decayedfromHiggs_energy;

//   for ArborPFO
        int num_pfo;
        int num_gamma_pfo, num_neutron_pfo, num_charge_pfo;
        int lepton_plus_count_pfo, lepton_minus_count_pfo;
		float total_energy_pfo, total_invmass_pfo, mass_dilepton_pfo, mass_left_pfo, momentum_pfo, mass_recoil_pfo;
        float reco_lepton_plus_energy, reco_lepton_minus_energy;
    int count_neutral_particle;
    float energy_neutron_first, energy_neutron_second, Mass_neutron_value;
    
		int _Num;

		unsigned int _eventNr;




		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


