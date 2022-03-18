#ifndef ABSRUNNINGSUMTPFINDERPASS1_H
#define ABSRUNNINGSUMTPFINDERPASS1_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "duneana/DAQSimAna/absRunningSumHitFinder/absRunningSumTPFinderTool.h"

class absRunningSumTPFinderPass1 : public absRunningSumTPFinderTool {
 public:
  explicit absRunningSumTPFinderPass1(fhicl::ParameterSet const & p);

  virtual std::vector<absRunningSumTPFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& adc_samples);
    

 protected:
  std::vector<short> downSample  (const std::vector<short>& orig);
  std::vector<short> findPedestal(const std::vector<short>& orig);
  std::vector<short> findPedestalAbsRS(const std::vector<short>& orig);
  std::vector<short> filter      (const std::vector<short>& orig);
  std::vector<short> absRS       (const std::vector<short>& orig);
  void hitFinding(const std::vector<short>& waveform,
		  std::vector<absRunningSumTPFinderTool::Hit>& hits, 
		  int channel);

  unsigned int       m_threshold;
  bool               m_useSignalKill;
  short              m_signalKillLookahead;
  short              m_signalKillThreshold;
  short              m_signalKillNContig;
  short              m_frugalNContig;
  short              m_frugalNContigAbsRS; 
  bool               m_doFiltering;
  unsigned int       m_downsampleFactor;
  std::vector<short> m_filterTaps;
  int                m_multiplier;
  float              m_R;
};

#endif
