
cet_build_plugin(dunezsanalysis art::module
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  nusimdata::SimulationBase
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  BASENAME_ONLY
  )

cet_build_plugin(dEdx art::module
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  lardataobj_RecoBase
  lardata_Utilities
  larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
  nusimdata::SimulationBase
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  BASENAME_ONLY
  )

cet_build_plugin(EMEnergyCalib art::module
  lardataalg_DetectorInfo
  lardataobj_RawData
  lardataobj_RecoBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
  lardata_Utilities
  nusimdata::SimulationBase
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  BASENAME_ONLY
  )

cet_build_plugin(EMPi0Energy art::module
  lardataobj_RawData
  lardataobj_RecoBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
  lardata_Utilities
  nusimdata::SimulationBase
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  BASENAME_ONLY
  )

install_headers()
install_fhicl()
install_source()
install_scripts()
