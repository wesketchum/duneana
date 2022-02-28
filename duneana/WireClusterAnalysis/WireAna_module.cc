#include "WireAna_module.h"
//Constructor for cnn struct


wireana::WireAna::WireAna(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}  ,
  fWireProducerLabel(pset.get< art::InputTag >("InputWireProducerLabel", "caldata")),
  fSimChannelLabel(pset.get< art::InputTag >("SimChannelLabel", "elecDrift")),
  fSimulationProducerLabel(pset.get< art::InputTag >("SimulationProducerLabel", "largeant"))
{
  fLogLevel           = pset.get<int>("LogLevel", 10);
  fDoAssns            = pset.get<bool>("DoAssns", false);
  fNChanPerApa        = pset.get<int>("ChannelPerApa", 2560);
  fNTicksPerWire      = pset.get<int>("TickesPerWire", 6000);
  auto const* geo = lar::providerFrom<geo::Geometry>();
  fNPlanes = geo->Nplanes();

  fMinPts  =  pset.get<unsigned int>("DBSCAN_MinPts", 2);
  fEps     =  pset.get<float>("DBSCAN_Eps", 4.5);
  fDrift   =  pset.get<float>("DBSCAN_Drift", 1.6);
  fPitch   =  pset.get<float>("DBSCAN_Pitch", 3.0);

}

void wireana::WireAna::analyze(art::Event const & evt) {
  //reset containers
  reset();
  //get detector property
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  //get event data
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  std::cout<<"########## EvtNo."<<event<<std::endl;

  MC = !evt.isRealData();

  ////////////////////////////////////////////////////////////
  //Build Wire SimChanel List
  //parse Wire/SimChannel Info
  //Sort wires/channels
  art::Handle<std::vector<recob::Wire>> wireListHandle;
  std::vector<art::Ptr<recob::Wire>> wirelist;
  if (evt.getByLabel(fWireProducerLabel, wireListHandle)) 
    art::fill_ptr_vector(wirelist, wireListHandle);
  SortWirePtrByChannel( wirelist, true );

  art::Handle<std::vector<sim::SimChannel>> simChannelListHandle;
  std::vector<art::Ptr<sim::SimChannel>> channellist;
  if (evt.getByLabel(fSimChannelLabel, simChannelListHandle)) 
    art::fill_ptr_vector(channellist, simChannelListHandle);
  SortWirePtrByChannel( channellist, true );

  // //Get channel-> wire,simchannel map
  for( auto w: wirelist ) 
    ch_w_sc[ w->Channel() ].first = w;
  for( auto w: channellist )
  {
    if ( ch_w_sc.find( w->Channel() ) != ch_w_sc.end() ) 
    {
      ch_w_sc[ w->Channel() ].second = w;
      if( fLogLevel>=10 ) 
      {
        std::cout<<"Filled Channel "<<w->Channel()<<std::endl;
      }
    }
  }

  /// FillTruthInfo to internal data objects
  /// i.e. trkid_to_label_map
  FillTrackIDtoLabelMap( evt );
  //We can now print truth particles

  //! Since ROI Tag Truth will supply the leading trkid, now need a method to convert trkid to MCPart and the truthlabel.
  //! The method to use seems to be PIS->TrackIdToMotherParticle(trkid
  //const simb::MCParticle *  TrackIdToParticle_P (int id) const
  //simb::MCParticle   TrackIdToParticle (int const id) const
  //const simb::MCParticle *  TrackIdToMotherParticle_P (int id) const
  //simb::MCParticle   TrackIdToMotherParticle (int const id) const
  //const art::Ptr< simb::MCTruth > &   TrackIdToMCTruth_P (int id) const
  //simb::MCTruth  TrackIdToMCTruth )
  

  /////////////////////////////////////////////////////////////
  // Build ROICluster
  BuildPlaneViewROIMap( wirelist );
  BuildInitialROIClusters();
  if ( MC )
  {
    TagAllROITruth( clockData );
  }

  /////////////////////////////////////////////////////////////
  // Build ROICluster
  // Set Truth Info on roicluster

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
  ch_w_sc.clear();
  trkid_to_label_map.clear();
  plane_view_roi_map.clear();
  plane_view_roicluster_map.clear();
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
  if( fLogLevel >= 3 ) 
  {
    std::cout<<"Entering SortWirePtrByChannel, sorting "<<vec.size()<<"channels."<<std::endl;
  }
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


void
wireana::WireAna::BuildPlaneViewROIMap(  std::vector<art::Ptr<recob::Wire>> &wires )
{
  //First get all rois
  //I assume the wires are already grouped by APA and view, because we don't want to
  //search over the entire DUNE volume all at once.
  for( const auto &wire: wires )
  {
    raw::ChannelID_t chan = wire->Channel();
    int planeID = chan/fNChanPerApa;
    recob::Wire::RegionsOfInterest_t signalROI = wire->SignalROI();
    for ( auto datarange: signalROI.get_ranges() ) {
      wireana::roi thisroi( wire, datarange );
      plane_view_roi_map[planeID][wire->View()].push_back( thisroi );
    }
  }
}

void wireana::WireAna::BuildInitialROIClusters()
{
  wireana::WireAnaDBSCAN scanner;
  for( auto &p: plane_view_roi_map ) //loop over planeid
  {
    //Internal clustering of ROIs
    for( auto &v: p.second ) //loop over view
    {
      scanner.SetParameters(v.second, fMinPts, fEps, fDrift, fPitch );
      scanner.run();
      std::map<int, wireana::roicluster> idclusmap;
      for( unsigned int i = 0; i<v.second.size(); i++ )
      {
        idclusmap[ v.second[i].clusterID ].AddROI( v.second[i], p.first );
      }
      for( auto idclus: idclusmap )
      {
        plane_view_roicluster_map[p.first][v.first].push_back(idclus.second);
      }
      std::sort(plane_view_roicluster_map[p.first][v.first].begin(), plane_view_roicluster_map[p.first][v.first].end(),[]( auto &a, auto &b ){ return a.nWires > b.nWires; } );
  //std::sort( trkID_sum.begin(), trkID_sum.end(), [](auto &a, auto &b){ return a.second.first > b.second.first;} );
    }
  }

  //if( fLogLevel >= 3 )
  //{
  //  for( const auto &p: plane_view_roi_map )
  //  {
  //    std::cout<<"==============================="<<std::endl;
  //    std::cout<<"Print Plane: "<< p.first <<std::endl;
  //    for( const auto &v: p.second )
  //    {
  //      std::cout<<"    View: "<<v.first<<std::endl;
  //      PrintROIs( v.second );
  //    }
  //  }
  //}

  if( fLogLevel >= 3 )
  {
    for( const auto &p: plane_view_roicluster_map)
    {
      std::cout<<"==============================="<<std::endl;
      std::cout<<"Print Plane: "<< p.first <<std::endl;
      for( const auto &v: p.second )
      {
        std::cout<<"    View: "<<v.first<<std::endl;
        for ( const auto &cluster: v.second ) PrintROIs( cluster.ROIs );
      }
    }
  }


  return;
}


void 
wireana::WireAna::PrintROIs( const std::vector<wireana::roi> &ROIs)
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



void
wireana::WireAna::FillTrackIDtoLabelMap( art::Event const& evt )
{
  ////////////////////////////////////////////////////////////
  //Build trackid -> truth particle map
  //parse MCParticles
  art::Handle<std::vector<simb::MCParticle>> particleHandle;
  if (!evt.getByLabel(fSimulationProducerLabel, particleHandle)) {
    throw cet::exception("AnalysisExample")
      << " No simb::MCParticle objects in this event - "
      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  } 

  /////////////////////////////////////////////////////////////
  // Get MCTruth handle
  // Get Trackid->generator info
  auto mcHandles = evt.getMany<std::vector<simb::MCTruth>>();

  //! Get a map of Art::Ptr<MCParticle> to the truth label
  for (auto const& mcHandle : mcHandles) {
    const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
    art::FindManyP<simb::MCParticle> findMCParts(mcHandle, evt, fSimulationProducerLabel);
    std::vector<art::Ptr<simb::MCParticle>> mcParts = findMCParts.at(0);
    for (const art::Ptr<simb::MCParticle> ptr : mcParts) {
      trkid_to_label_map[ptr->TrackId()] = sModuleLabel;
    }
  }
}





void 
wireana::WireAna::TagAllROITruth(const detinfo::DetectorClocksData &clock )
{
  for( auto &pvv : plane_view_roicluster_map )
  {
    for ( auto & vv : pvv.second )
    {
      for (auto & cluster : vv.second )
      {
        TagROITruth( cluster, clock );
      }
    }
  }
  if (fLogLevel >= 3 )
  {
    std::cout<<"========================================="<<std::endl;
    std::cout<<"Printing Truth:"<<std::endl;
    for( auto &pvv : plane_view_roicluster_map )
    {
      for ( auto & vv : pvv.second )
      {
        std::cout<<Form("Plane %d, View %d: ", pvv.first,vv.first)<<std::endl;
        for (auto & cluster : vv.second )
        {
          std::cout<<
            Form("  Cluster ID: %d, (wMin, wMax):(%d, %d), (tMin, tMax):(%d, %d)",
              cluster.clusID,
              cluster.channel_min, cluster.channel_max, 
              cluster.begin_index, cluster.end_index )
            <<std::endl;
          std::cout<<
            Form("   Truth Label: %s, Leading PDG: %d, Leading #electron: %f",
               cluster.label.c_str(), cluster.pdg_energy_list.front().first,
                cluster.pdg_energy_list.front().second.first)
            <<std::endl;
        }
      }
    }
  }

}

//! Tagging ROI Truth --
//! Look at each tick and search for true particle deposit and the
//! generator tag
void
wireana::WireAna::TagROITruth( wireana::roicluster &cluster, const detinfo::DetectorClocksData &clock )
{
  if (fLogLevel>=3) std::cout<<"Entering TagROITruth"<<std::endl;

  std::vector< std::pair< int, std::pair<float,float> > > trkID_sum; 
  for ( auto roi: cluster.ROIs ) // loop 1
  {
    int roi_channel = roi.wire->Channel();
    if (fLogLevel>=10) 
    {
      std::cout<<"look at channel "<<roi_channel<<std::endl;
      std::cout<<"Map Size: "<<ch_w_sc.size()<<std::endl;
      std::cout<<"Channel ID Exists? "<< (ch_w_sc.find(roi_channel) != ch_w_sc.end())<<std::endl;
      std::cout<<"Wire Exist? "<<ch_w_sc[ roi_channel ].first->Channel()<<std::endl;
      std::cout<<"SimChannel Exist? "<<ch_w_sc[ roi_channel ].second->Channel()<<std::endl;
    }
    art::Ptr<sim::SimChannel> &simchannel = ch_w_sc[ roi_channel ].second;

    //! Get the TPC time from the ticks in the ROI datarange_t
    double startTDC = clock.TPCTick2TDC( roi.begin_index );
    double endTDC   = clock.TPCTick2TDC( roi.end_index );
    if (fLogLevel>=10) 
    {
      std::cout<<"SimChannel exist? "<<simchannel->Channel()<<std::endl;
      std::cout<<"startTDC, endTDC: "<<startTDC<<", "<<endTDC<<std::endl;
    }

    //! Get all the IDE (Ionization at a point of the TPC sensitive volume. ) that deposited energy in the datarange_t
    std::vector< sim::IDE >  ide_vec = simchannel->TrackIDsAndEnergies(startTDC, endTDC);

    if( ide_vec.size() == 0 ) continue;

    //! This ROI has true signal
    roi.hasTrueSignal = true;

    //!pair<trkid, pair<numElectrons, energy>> , sum up n electrons and total energy within the ROI for each trkid

    //! fill vector of ide
    for( auto ide : ide_vec ) //loop 5
    {
      //! check if track has been saved in the trkID_sum vector.
      auto it = trkID_sum.begin();
      for( ; it!=trkID_sum.end(); ++it ) {
        if( it->first == ide.trackID ) break; 
      }

      if (it == trkID_sum.end()) {
        //! if not, add the track id and initialize with the IDE data
        std::pair< int, std::pair<float,float> > data;
        data.first = ide.trackID;
        data.second.first =ide.numElectrons;
        data.second.second = ide.energy;
        trkID_sum.push_back(data);
      } else {
        //! if yes, sum nElectron and energy to the track id.
        it->second.first+=ide.numElectrons;
        it->second.second+=ide.energy;
      }
    }// end loop 5

    //! sort vector of (trkid,(#electron, energy)) by # electrons
    //! could also sort by energy but #e should be a more direct observable?
  }//end loop 1


  // Once all ROI has been analyzed, sort trkID by #electrons
  std::sort( trkID_sum.begin(), trkID_sum.end(), [](auto &a, auto &b){ return a.second.first > b.second.first;} );

  cluster.trkID_sum = trkID_sum;
  cluster.label = this->trkid_to_label_map[ trkID_sum[0].first ];

  //Grep the most deposited particle first
  int trkid = trkID_sum.front().first;
  const simb::MCParticle * part = PIS->TrackIdToMotherParticle_P(trkid);
  art::Ptr<simb::MCTruth> truth = PIS->TrackIdToMCTruth_P(trkid);
  std::string label = this->trkid_to_label_map[trkid];
  if ( truth->NeutrinoSet() )
  {
    cluster.truthFromNeutrino = true;
    const simb::MCNeutrino neutrino = truth->GetNeutrino();
    cluster.momentum_neutrino = neutrino.Nu().Momentum();
  }
  cluster.momentum_part = part->Momentum();

  for( auto &trkid : trkID_sum )
  {
    auto part1 = PIS->TrackIdToMotherParticle_P(trkid.first);
    int pdg = part1->PdgCode();
    cluster.pdg_energy_list.push_back( std::make_pair( pdg, trkid.second ));
  }

}


DEFINE_ART_MODULE(wireana::WireAna)
