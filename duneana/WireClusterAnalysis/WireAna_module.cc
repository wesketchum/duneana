#include "WireAna_module.h"
//Constructor for cnn struct


wireana::WireAna::WireAna(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}  ,
  fWireProducerLabel(pset.get< art::InputTag >("InputWireProducerLabel")),
  fSimChannelLabel(pset.get< art::InputTag >("InputSimChannelLabel", "elecDrift"))
{
  fLogLevel           = pset.get<int>("LogLevel", 10);
  fDoAssns            = pset.get<bool>("DoAssns", false);
  fNChanPerApa        = pset.get<int>("ChannelPerApa", 2560);
  fNTicksPerWire      = pset.get<int>("TickesPerWire", 6000);
  fIsMC               = pset.get<bool>("IsMC", true);
  auto const* geo = lar::providerFrom<geo::Geometry>();
  fNPlanes = geo->Nplanes();


}

void wireana::WireAna::analyze(art::Event const & evt) {
  //reset containers
  reset();


  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  std::cout<<"########## EvtNo."<<event<<std::endl;

  if( !evt.isRealData() ) MC = 1;
  else MC = 0;

  art::Handle<std::vector<recob::Wire>> wireListHandle;
  std::vector<art::Ptr<recob::Wire>> wirelist;
  if (evt.getByLabel(fWireProducerLabel, wireListHandle)) art::fill_ptr_vector(wirelist, wireListHandle);
  art::Handle<std::vector<sim::SimChannel>> simChannelListHandle;
  std::vector<art::Ptr<sim::SimChannel>> channellist;
  std::vector<int> matchedChannelID;
  if (evt.getByLabel(fSimChannelLabel, simChannelListHandle)) art::fill_ptr_vector(channellist, simChannelListHandle);

  //First sort wires by channel:
  SortWirePtrByChannel( wirelist, true );
  SortWirePtrByChannel( channellist, true );

  //Get channel-> wire,simchannel map
  std::map<raw::ChannelID_t, std::pair<art::Ptr<recob::Wire>, art::Ptr<sim::SimChannel>>> ch_w_sc;
  for( auto w: wirelist ) ch_w_sc[ w->Channel() ].first = w;
  for( auto w: channellist )
  {
    raw::ChannelID_t c = w->Channel();
    if ( ch_w_sc.find( c ) != ch_w_sc.end() ) ch_w_sc[c].second = w;
  }

  //We can now print truth particles


  std::vector<wireana::wirecluster> clusters = BuildInitialClusters( wirelist, 0, 0 );
  //PrintClusters( clusters );

  std::vector<wireana::roicluster> roi_clusters = BuildInitialROIClusters( wirelist, 0, 0 );

  if (fIsMC)
  {
    //Fill Truth info to roicluster

  }

  truth_data.truth_intType = 1;

  fTree->Fill();
  std::unique_ptr<std::vector<recob::Wire>> outwires(new std::vector<recob::Wire>);
}

void wireana::WireAna::beginJob() {

  gROOT->SetBatch(1);

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("wireana","wireana tree");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("MC", &MC);
  DeclareTruthBranches(fTree, truth_data);
}

void wireana::WireAna::endJob()
{
}



void wireana::WireAna::reset()
{
  ResetTruthBranches(this->truth_data);
}



//======================================================================================
//======================================================================================
void wireana::WireAna::DeclareTruthBranches(TTree*t, DataBlock_Truth &blk)
{
  t->Branch("truth_intType",        &blk.truth_intType);         
  t->Branch("truth_nu_momentum_x",  &blk.truth_nu_momentum_x); 
  t->Branch("truth_nu_momentum_y",  &blk.truth_nu_momentum_y); 
  t->Branch("truth_nu_momentum_z",  &blk.truth_nu_momentum_z); 
  t->Branch("truth_nu_momentum_e",  &blk.truth_nu_momentum_e); 
  t->Branch("truth_nu_vtx_x",       &blk.truth_nu_vtx_x); 
  t->Branch("truth_nu_vtx_y",       &blk.truth_nu_vtx_y); 
  t->Branch("truth_nu_vtx_z",       &blk.truth_nu_vtx_z); 
  t->Branch("truth_nu_vtx_t",       &blk.truth_nu_vtx_t); 
  t->Branch("truth_nu_PDG",         &blk.truth_nu_PDG); 
  t->Branch("truth_lep_momentum_x", &blk.truth_lep_momentum_x); 
  t->Branch("truth_lep_momentum_y", &blk.truth_lep_momentum_y); 
  t->Branch("truth_lep_momentum_z", &blk.truth_lep_momentum_z); 
  t->Branch("truth_lep_momentum_e", &blk.truth_lep_momentum_e); 
  t->Branch("truth_lep_vtx_x",      &blk.truth_lep_vtx_x); 
  t->Branch("truth_lep_vtx_y",      &blk.truth_lep_vtx_y); 
  t->Branch("truth_lep_vtx_z",      &blk.truth_lep_vtx_z); 
  t->Branch("truth_lep_vtx_t",      &blk.truth_lep_vtx_t); 
  t->Branch("truth_lep_PDG",        &blk.truth_lep_PDG); 
  t->Branch("truth_had_momentum_x", &blk.truth_had_momentum_x); 
  t->Branch("truth_had_momentum_y", &blk.truth_had_momentum_y); 
  t->Branch("truth_had_momentum_z", &blk.truth_had_momentum_z); 
  t->Branch("truth_had_momentum_e", &blk.truth_had_momentum_e); 
  t->Branch("truth_had_vtx_x",      &blk.truth_had_vtx_x); 
  t->Branch("truth_had_vtx_y",      &blk.truth_had_vtx_y); 
  t->Branch("truth_had_vtx_z",      &blk.truth_had_vtx_z); 
  t->Branch("truth_had_vtx_t",      &blk.truth_had_vtx_t); 
  t->Branch("truth_had_PDG",        &blk.truth_had_PDG); 
}

void wireana::WireAna::ResetTruthBranches(DataBlock_Truth &blk)
{
  blk.truth_intType = DEFAULT_VALUE;         
  blk.truth_nu_momentum_x = DEFAULT_VALUE; 
  blk.truth_nu_momentum_y = DEFAULT_VALUE; 
  blk.truth_nu_momentum_z = DEFAULT_VALUE; 
  blk.truth_nu_momentum_e = DEFAULT_VALUE; 
  blk.truth_nu_vtx_x = DEFAULT_VALUE; 
  blk.truth_nu_vtx_y = DEFAULT_VALUE; 
  blk.truth_nu_vtx_z = DEFAULT_VALUE; 
  blk.truth_nu_vtx_t = DEFAULT_VALUE; 
  blk.truth_nu_PDG = DEFAULT_VALUE; 
  blk.truth_lep_momentum_x = DEFAULT_VALUE; 
  blk.truth_lep_momentum_y = DEFAULT_VALUE; 
  blk.truth_lep_momentum_z = DEFAULT_VALUE; 
  blk.truth_lep_momentum_e = DEFAULT_VALUE; 
  blk.truth_lep_vtx_x = DEFAULT_VALUE; 
  blk.truth_lep_vtx_y = DEFAULT_VALUE; 
  blk.truth_lep_vtx_z = DEFAULT_VALUE; 
  blk.truth_lep_vtx_t = DEFAULT_VALUE; 
  blk.truth_lep_PDG = DEFAULT_VALUE; 
  blk.truth_had_momentum_x = DEFAULT_VALUE; 
  blk.truth_had_momentum_y = DEFAULT_VALUE; 
  blk.truth_had_momentum_z = DEFAULT_VALUE; 
  blk.truth_had_momentum_e = DEFAULT_VALUE; 
  blk.truth_had_vtx_x = DEFAULT_VALUE; 
  blk.truth_had_vtx_y = DEFAULT_VALUE; 
  blk.truth_had_vtx_z = DEFAULT_VALUE; 
  blk.truth_had_vtx_t = DEFAULT_VALUE; 
  blk.truth_had_PDG = DEFAULT_VALUE; 
}


void wireana::WireAna::FillTruthBranches(art::Event const& evt, TTree*t, DataBlock_Truth &blk)
{


}

template<class T>
void wireana::WireAna::SortWirePtrByChannel( std::vector<art::Ptr<T>> &vec, bool increasing )
{
  if (increasing)
  {
    std::sort(vec.begin(), vec.end(), [](art::Ptr<T> &a, art::Ptr<T> &b) { return a->Channel() < b->Channel(); });
  }
  else
  {
    std::sort(vec.begin(), vec.end(), [](art::Ptr<T> &a, art::Ptr<T> &b) { return a->Channel() > b->Channel(); });
  }
}

template void wireana::WireAna::SortWirePtrByChannel<>( std::vector<art::Ptr<recob::Wire>> &vec, bool increasing );
template void wireana::WireAna::SortWirePtrByChannel<>( std::vector<art::Ptr<sim::SimChannel>> &vec, bool increasing );

std::vector<wireana::wirecluster> 
wireana::WireAna::BuildInitialClusters( std::vector<art::Ptr<recob::Wire>> &vec, int dC, int dTick )
{
  std::vector<wireana::wirecluster> ret;

  bool newcluster = true;
  for( auto it = vec.begin(); it!=vec.end(); )
  {
    if (newcluster) 
    {
      wirecluster clus;
      ret.push_back( clus );
      ret.back().push_back( *it );
      newcluster = false;
      ++it;
      continue;
    }

    if( abs( int( (*it)->Channel() - ret.back().back()->Channel()) ) == 1 )
    {
      ret.back().push_back( *it );
      ++it;
    } else 
    {
      newcluster = true;
    }
  }

  //now ret is a vector of wirecluster object
  //each wirecluster contains set of wires adjacent to each other
  return ret;
}



bool 
wireana::WireAna::HasHit( const art::Ptr<recob::Wire> &wire, int minTick )
{
  if( wire->SignalROI().n_ranges() == 0 ) return false;
  for( auto itr = wire->SignalROI().begin_range(); itr!=wire->SignalROI().end_range(); ++itr )
  {
    int dTick = itr->end_index() - itr->begin_index();
    if (dTick > minTick) return true;
  }
  return false;

}

void wireana::WireAna::PrintClusters( std::vector<wirecluster> &clusters )
{
  for( auto it = clusters.begin(); it != clusters.end(); it++ )
  {
    std::cout<<"Printing Cluster "<<it-clusters.begin()<<std::endl;
    std::cout<<Form("Nwires: %d, Min Channel: %d, Max Channel: %d", it->nWires, it->wire_pointers.front()->Channel(), it->wire_pointers.back()->Channel() )<<std::endl;
  }

}


std::vector<wireana::roicluster> 
wireana::WireAna::BuildInitialROIClusters( std::vector<art::Ptr<recob::Wire>> &wires, int dW, int dTick )
{
  //First get all rois
  //I assume the wires are already grouped by APA and view, because we don't want to
  //search over the entire DUNE volume all at once.
  std::vector<wireana::roicluster> ret;
  std::vector<wireana::roi> allROIs;
  wireana::PlaneViewROIMap plane_view_roi_map;
  wireana::WireAnaDBSCAN scanner;
  for( auto wire: wires )
  {
    raw::ChannelID_t chan = wire->Channel();
    int planeID = chan/fNChanPerApa;
    recob::Wire::RegionsOfInterest_t signalROI = wire->SignalROI();
    for ( auto datarange: signalROI.get_ranges() ) {
      plane_view_roi_map[planeID][wire->View()].push_back(wireana::roi( wire, datarange ) );
    }
      //allROIs.push_back( wireana::roi( wire, datarange ) );

  }

  for( auto p: plane_view_roi_map )
  {
    std::cout<<"==============================="<<std::endl;
    std::cout<<"Print Plane: "<< p.first <<std::endl;
    for( auto v: p.second )
    {
      scanner.SetParameters(v.second,2,4.5,1.6,3.0);
      scanner.run();
      std::cout<<"    View: "<<v.first<<std::endl;
      //PrintROIs( scanner.GetROIs() );
      PrintROIs( v.second );
    }

  }
  return ret;

}


void 
wireana::WireAna::PrintROIs( std::vector<wireana::roi> &ROIs)
{
  std::cout<<"Printing ROIs"<<std::endl;
  for( auto roi: ROIs )
  {
    if (roi.clusterID == NOISE ) continue;
    std::cout<<Form(
        "\t\tClusID: %d \n\t\t\t\tChannel: %d, Index: (%d,%d,%d), |Sum|: %f", 
        roi.clusterID, 
        roi.channel, roi.begin_index, roi.end_index, roi.end_index-roi.begin_index,roi.abs_sum
        )
      <<std::endl;
  }

}



DEFINE_ART_MODULE(wireana::WireAna)
