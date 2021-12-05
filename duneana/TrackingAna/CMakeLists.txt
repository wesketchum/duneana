art_make(BASENAME_ONLY
	 MODULE_LIBRARIES
         lardataalg_DetectorInfo
	 lardataobj_RawData
	 lardataobj_RecoBase
	 larcorealg_Geometry
	 larcore_Geometry_Geometry_service
	 lardata_Utilities
	 larsim_MCCheater_BackTrackerService_service
	 larsim_MCCheater_ParticleInventoryService_service
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
	 ROOT_GEOM
)

install_headers()
install_fhicl()
install_scripts()
