#ifndef _BushConnect_hh_
#define _BushConnect_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <IMPL/LCFlagImpl.h>
#include <TNtuple.h>
#include <TObject.h>

#include <TTree.h>
#include <TFile.h>
#include <TH3.h>
#include <TVector3.h>

class TTree;

class BushConnect  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new BushConnect ; }

		BushConnect();

		~BushConnect() {};

		void init();

		void Clean();

		void MCPTrackSort(LCEvent* evtPP);

		void TrackSort(EVENT::LCEvent* evtPP);

		void BushSelfMerge(EVENT::LCEvent* evtPP);

		void TagCore(EVENT::LCEvent* evtPP);

		void ParticleReco(LCEvent * evtPP);

		void processEvent( LCEvent * evtP );   

		void end();


	protected:

		std::vector<std::string> _BushCollections; 
		std::vector<std::string> _HitCollections;
		LCFlagImpl Cluflag;
		std::ostream *_output;
		std::vector<Cluster*> SortedSMBushes;
		std::vector<Track*> SortedTracks;
		std::map<Track*, float> Track_Energy;
		std::map<Track*, TVector3> Track_P3;
		std::map<Track*, int> Track_Type;
		std::map<Track*, float> Track_Theta;
		std::map<Track*, float> Track_Phi;	

		std::map<Cluster*, int> ClusterType_1stID;
		std::map<ReconstructedParticle*, int> ChCoreID; 

		std::vector<Cluster*> ecalchcore_tight;         //TightCores
		std::vector<Cluster*> ecalchcore_medium;
		std::vector<Cluster*> ecalchcore_loose;         //LooseCores    Let's also get
		std::vector<Cluster*> ecalchcore; 		    //Above three
		std::vector<Cluster*> ecalnecore;
		std::vector<Cluster*> ecalnecore_EM;
		std::vector<Cluster*> ecalnecore_NonEM;
		std::vector<Cluster*> ecalfrag;
		std::vector<Cluster*> ecalundef;
		std::vector<Cluster*> ecalfrag_TBM_CH;
		std::vector<Cluster*> ecalfrag_TBM_NE_EM;
		std::vector<Cluster*> ecalfrag_TBM_NE_NonEM;
		std::vector<Cluster*> ecalundef_iso;
		std::vector<Cluster*> ecalpotentialbackscattering;

		std::vector<Cluster*> chargedclustercore;
		std::vector<Cluster*> chargedclustercore_abs;

		std::vector<Cluster*> selfmergedcluster; 
		std::vector<Cluster*> non_chargedclustercore;
		std::vector<Cluster*> onlyNeutralCore;

		std::vector<Cluster*> non_charged_pem_neutral_core;
		std::vector<Cluster*> pem_neutral_core;

		std::map<Track*, int>MCPTrack_Type;
		std::map<Track*, TVector3> Track_EndPoint;       //Last hit
		std::map<Track*, TVector3> TrackStartPoint;
		std::map<Cluster*, float> CluFD; 
		std::map<Cluster*, float> CluT0;
		std::map<Cluster*, float> Clu_Depth; 
		std::map<Cluster*, TVector3> CluCoG;
		float _HitAbsCut, _TimeCut;
		int _FLAG_DIAGNOSIS, _FLAG_MCPMIMIC; 

		int _eventNr;		

};

#endif


