////////////////////////////////////////////////////////////////////////
// Class:       absRunningSumTPFinderPass2
// 
// Hit finder for induction/collection channels passed through the absRS algorithm 
// Hit finding done when signals are x*sigma above noise
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "duneana/DAQSimAna/absRunningSumHitFinder/absRunningSumTPFinderPass1.h"
#include "duneana/DAQSimAna/absRunningSumHitFinder/absRSAlgParts.h"
#include <algorithm> // for std::transform


class absRunningSumTPFinderPass2 : public absRunningSumTPFinderPass1 {
public:
  explicit absRunningSumTPFinderPass2(fhicl::ParameterSet const & p);

  virtual std::vector<absRunningSumTPFinderTool::Hit>
  findHits(const std::vector<unsigned int>& channel_numbers,
	   const std::vector<std::vector<short>>& adc_samples);


protected:
  void hitFinding(const std::vector<short>& waveform,
                  const std::vector<short>& iqr_abs,
                  std::vector<absRunningSumTPFinderTool::Hit>& hits,
                  int channel);

private:
  float m_sigmaThreshold;

};


absRunningSumTPFinderPass2::absRunningSumTPFinderPass2(fhicl::ParameterSet const & p)
  : absRunningSumTPFinderPass1(p),
    m_sigmaThreshold(p.get<float>  ("ThresholdInSigma", 5))
    
{
  std::cout << "Threshold in sigma is " << m_sigmaThreshold << " (ignore the static threshold on previous line)\n"; 
}


void
absRunningSumTPFinderPass2::hitFinding(const std::vector<short>& waveform,
				       const std::vector<short>& iqr,
				    std::vector<absRunningSumTPFinderTool::Hit>& hits,
				    int channel) {

  //---------------------------------------------
  // Hit finding
  //---------------------------------------------
  bool is_hit  = false;
  bool was_hit = false;
  absRunningSumTPFinderTool::Hit hit(channel, 0, 0, 0);
  for(size_t isample=0; isample<waveform.size()-1; ++isample){
    int   sample_time = isample * m_downsampleFactor;
    short adc         = waveform[isample];
    //check if we have a hit (adc above threshold)
    is_hit = (float)adc >m_sigmaThreshold*iqr[isample];
    if(is_hit && !was_hit) {
      hit.startTime         = sample_time;
      hit.charge            = adc;
      hit.timeOverThreshold = m_downsampleFactor;
    }
    if(is_hit && was_hit) {
      hit.charge            += adc*m_downsampleFactor;
      hit.timeOverThreshold += m_downsampleFactor;
    }
    if(!is_hit && was_hit) {
      hit.charge /= m_multiplier;
      hits.push_back(hit);
    }
    was_hit = is_hit;
  }
}

std::vector<absRunningSumTPFinderTool::Hit>
absRunningSumTPFinderPass2::findHits(const std::vector<unsigned int>& channel_numbers, 
				  const std::vector<std::vector<short>>& adc_samples) {

  auto hits = std::vector<absRunningSumTPFinderTool::Hit>();
  std::cout << "findHits called with "      << adc_samples.size()
	    << " channels. First chan has " << adc_samples[0].size() << " samples" << std::endl;

  for(size_t ich=0; ich<adc_samples.size(); ++ich){

    const std::vector<short>& waveformOrig = adc_samples[ich];
    std::vector<short> waveform  = downSample  (waveformOrig);
    std::vector<short> pedestal  = findPedestal(waveform    );
    std::vector<short> pedsub(waveform.size(), 0);
    for(size_t i=0; i<pedsub.size(); ++i){
      pedsub[i]=waveform[i]-pedestal[i];
    }
    std::vector<short> filtered  = filter(pedsub);
    //get the integrated waveform, subtract the pedestal again +  run the hit finding 
    std::vector<short> absRSwaveform = absRS(filtered);
    std::vector<short> absRSpedestal  = findPedestalAbsRS(absRSwaveform);
    std::vector<short> absRSwaveform_pedsub(absRSwaveform.size(), 0);
    for(size_t i=0; i<absRSwaveform_pedsub.size(); ++i){
      absRSwaveform_pedsub[i] = absRSwaveform[i]-absRSpedestal[i];
    }

    std::vector<short> iqr_absRS = frugal_iqr_absRS(absRSwaveform_pedsub,absRSpedestal, m_frugalNContigAbsRS);

    hitFinding( absRSwaveform_pedsub, iqr_absRS, hits, channel_numbers[ich]);
  }
  std::cout << "Returning " << hits.size() << " hits" << std::endl;
  std::cout << "hits/channel=" << float(hits.size())/adc_samples   .size() << std::endl;
  std::cout << "hits/tick="    << float(hits.size())/adc_samples[0].size() << std::endl;
  return hits;
}

DEFINE_ART_CLASS_TOOL(absRunningSumTPFinderPass2)
