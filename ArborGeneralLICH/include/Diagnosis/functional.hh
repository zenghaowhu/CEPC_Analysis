#ifndef _functional_hh_
#define _functional_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class functional  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new functional ; }

		functional();

		~functional() {};

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

        int _Num;
        unsigned int _eventNr;

    //MCParticle
    
        float ISREn, ISRPt;
		float NEn, NPt, NEn_mcp, NPt_mcp;
        int type;
        int PID1_decayedfromH, PID2_decayedfromH;
        int PID2_decayedfromZ2, PID1_decayedfromZ1, PID2_decayedfromZ1, PID1_decayedfromZ2;
        float mcp_Z1D1_energy, mcp_Z1D2_energy, mcp_Z2D1_energy, mcp_Z2D2_energy;
		int OriQuarkID;
		float HDir;
		float MaxJetCosTheta;
        float mcp_muminus_energy, mcp_muplus_energy, mcp_HD1_energy, mcp_HD2_energy, mcp_HDplus_energy, mcp_HDminus_energy;
        float mass_dimuon1_mcp;
  
    
    //ArborPFO
        int NMuplusCount, NMuminusCount, NNeutronCount;
        float energy_neutron_first, energy_neutron_second, energy_muminus, energy_muplus,energy_muminus_H, energy_muplus_H;
        float angle_minus, angle_plus, H_angle_minus, H_angle_plus;
        float Mass_neutron_value;
        float Mass_a;
        float mass_dilepton1_pfo, mass_dilepton2_pfo, mass_recoil1, mass_recoil2;
        float missing_mass, visible_mass;
    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


