////////////////////////////////////////////////////////////////////////
// Class:      WireAna 
// Plugin Type: analyzer (art v3_00_00)
// File:        WireAna_module.cc
// Written by Tejin Cai
// Reach out for questions/issues/bugs
////////////////////////////////////////////////////////////////////////
#ifndef WIREANA_H
#define WIREANA_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardataobj/RawData/RDTimeStamp.h"


#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"

#include "art_root_io/TFileService.h"

#include "c2numpy.h"

// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "WireAna_Utils.h"

//Others
#define DEFAULT_VALUE -99999

namespace wireana {

  //using recob::SpacePoint;



  struct DataBlock_Truth
  {
    //Truth information
    //Incident Particle
    int truth_intType; //ES,CC,Rad: 0,1,2
    double truth_nu_momentum_x;
    double truth_nu_momentum_y;
    double truth_nu_momentum_z;
    double truth_nu_momentum_e;
    double truth_nu_vtx_x;
    double truth_nu_vtx_y;
    double truth_nu_vtx_z;
    double truth_nu_vtx_t;
    double truth_nu_PDG;
    // // Electron Info
    double truth_lep_momentum_x;
    double truth_lep_momentum_y;
    double truth_lep_momentum_z;
    double truth_lep_momentum_e;
    double truth_lep_vtx_x;
    double truth_lep_vtx_y;
    double truth_lep_vtx_z;
    double truth_lep_vtx_t;
    double truth_lep_PDG;
    // // Hadron Info
    double truth_had_momentum_x;
    double truth_had_momentum_y;
    double truth_had_momentum_z;
    double truth_had_momentum_e;
    double truth_had_vtx_x;
    double truth_had_vtx_y;
    double truth_had_vtx_z;
    double truth_had_vtx_t;
    double truth_had_PDG;
  };


  class WireAna;

}

class wireana::WireAna : public art::EDAnalyzer {
public:
  explicit WireAna(fhicl::ParameterSet const& pset);
  WireAna(WireAna const&) = delete;
  WireAna(WireAna&&) = delete;
  WireAna& operator=(WireAna const&) = delete;
  WireAna& operator=(WireAna&&) = delete;

  /////////////////////////////////////////////
  // Required functions.
  void analyze(art::Event const& evt) override;

  /////////////////////////////////////////////
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  /////////////////////////////////////////////
  void reset();

private:
  /////////////////////////////////////////////
  // Geometry Options && Tool options
  int fNPlanes;
  int fNChanPerApa;
  int fNTicksPerWire;
  int fChannelDistance;
  int fTickDistance;
  int fMinClusterSize;

  unsigned int fMinPts;
  float fEps;
  float fDrift;
  float fPitch;

  float fDeltaMetric;

  int image_channel_width;
  int image_tick_width;
  int image_rebin_tick;
  int image_size;


  /////////////////////////////////////////////
  // Backtracker services

  /////////////////////////////////////////////
  // config
  int fLogLevel;
  bool fDoAssns;

  /////////////////////////////////////////////
  // Wire Filtering/Clustering Functions
  template<class T>
  void SortWirePtrByChannel( std::vector<art::Ptr<T>> &vec, bool increasing );
  std::vector<wireana::wirecluster> BuildInitialClusters( std::vector<art::Ptr<recob::Wire>> &vec, int dW, int dTick );


  void BuildPlaneViewROIMap(  std::vector<art::Ptr<recob::Wire>> &wires );
  void BuildInitialROIClusters();
  std::vector<art::Ptr<recob::Wire>> FilterWires(std::vector<art::Ptr<recob::Wire>> &vec, int dC1, int dT1, int dCn, int dTn );
  bool HasHit( const art::Ptr<recob::Wire> &wire, int minTick );


  void PrintClusters( std::vector<wirecluster> &clusters );
  void PrintROIs( const std::vector<roi> &ROIs);

  void WriteNumPy( matchedroicluster& cluster, std::vector<art::Ptr<recob::Wire>>& wirelist );

  /////////////////////////////////////////////
  // Cluster Matching

  /////////////////////////////////////////////
  // Declare output data
  TTree *fTree;





  /////////////////////////////////////////////
  //Module Labels and Settting
  const art::InputTag fWireProducerLabel; 
  const art::InputTag fSimChannelLabel;
  const art::InputTag fSimulationProducerLabel;

  std::string fDumpFileName;
  int fDumpMaxRow;


  /////////////////////////////////////////////
  // Event level information
  int run;
  int subrun;
  int event;
  int MC;
  int signal_Apa;
 
  ////////////////////////////////////////////
  // Truth Operation
  // declare truth branch
  wireana::DataBlock_Truth truth_data;

  void DeclareTruthBranches(TTree*t, DataBlock_Truth &blk);
  void FillTruthBranches(art::Event const& evt, TTree*t, DataBlock_Truth &blk);
  void ResetTruthBranches(DataBlock_Truth &blk);

  // Truth tagging
  void FillTrackIDtoLabelMap( art::Event const& evt );
  void TagAllROITruth( const detinfo::DetectorClocksData &clock );
  void TagROITruth( wireana::roicluster &roi, const detinfo::DetectorClocksData &clock );
  art::ServiceHandle<cheat::ParticleInventoryService> PIS;

  ////////////////////////////////////////////
  // Internal Data Structure
  std::map<raw::ChannelID_t, std::pair<art::Ptr<recob::Wire>, art::Ptr<sim::SimChannel>>> 
    ch_w_sc;
  std::map<int, std::string> 
    trkid_to_label_map;
  PlaneViewROIMap plane_view_roi_map; //type is std::map<int, std::map< geo::View_t, std::vector<roi> > >
  PlaneViewROIClusterMap plane_view_roicluster_map; //type is std::map<int, std::map< geo::View_t, std::vector<roi> > >

  ////////////////////////////////////////////
  // Image Forming:
  std::vector<double> GetArrayFromWire( std::vector<art::Ptr<recob::Wire>> &wirelist, wireana::roicluster &cluster, int channel_width, int new_tickwidth );
  std::vector<double> CombineTicks( const std::vector<double> &input, int channel_width, int nticks);
  std::vector<double> ScaleArray( const std::vector<double> &input, double min, double max );
  int CalculateIndex( int c, int t, int c_width, int t_width );
  std::pair<int,int> CalculateCT( int index, int c_width );

  std::map<int, std::string> fViewMap;

  c2numpy_writer npywriter;
};





#endif
