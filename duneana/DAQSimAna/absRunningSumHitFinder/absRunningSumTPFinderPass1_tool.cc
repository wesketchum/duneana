////////////////////////////////////////////////////////////////////////
// Class:       absRunningSumTPFinderPass1
// Plugin Type: service (art v2_10_03)
// File:        absRunningSumTPFinderPass1_service.cc
//
// Hit finder for induction/collection channels passed through the absRS algorithm 
// Hit finding done using a set ADC threshold 
////////////////////////////////////////////////////////////////////////

#include "duneana/DAQSimAna/absRunningSumHitFinder/absRunningSumTPFinderPass1.h"
#include "duneana/DAQSimAna/absRunningSumHitFinder/absRSAlgParts.h"
#include <algorithm> // for std::transform
#include <numeric> // for std::accumulate


absRunningSumTPFinderPass1::absRunningSumTPFinderPass1(fhicl::ParameterSet const & p)
  : m_threshold          (p.get<unsigned int>      ("Threshold"            , 200)),
    m_useSignalKill      (p.get<bool>              ("UseSignalKill"        , true)),
    m_signalKillLookahead(p.get<short>             ("SignalKillLookahead"  , 5)),
    m_signalKillThreshold(p.get<short>             ("SignalKillThreshold"  , 15)),
    m_signalKillNContig  (p.get<short>             ("SignalKillNContig"    ,  1)),
    m_frugalNContig      (p.get<short>             ("FrugalPedestalNContig", 10)),
    m_frugalNContigAbsRS (p.get<short>             ("FrugalPedestalNContigAbsRS",  15)),
    m_doFiltering        (p.get<bool>              ("DoFiltering"          ,  true)),
    m_downsampleFactor   (p.get<unsigned int>      ("DownsampleFactor"     ,  1)),
    m_filterTaps         (p.get<std::vector<short>>("FilterCoeffs"         , {2,  9, 23, 31, 23,  9,  2})),
    m_multiplier         (std::accumulate(m_filterTaps.begin(), m_filterTaps.end(), 0)),
    m_R                  (p.get<float>             ("R"                    , 0.81))
    
{

}

std::vector<short> absRunningSumTPFinderPass1::downSample(const std::vector<short>& orig) {

  //---------------------------------------------
  // Do the downsampling
  //---------------------------------------------
  if (m_downsampleFactor==1) {
    return orig;
  }
  else {
    std::vector<short> waveform;
    for(size_t i=0; i<orig.size(); i+=m_downsampleFactor) {
      waveform.push_back(orig[i]);
    }
    return waveform;
  }
}

std::vector<short> absRunningSumTPFinderPass1::findPedestal(const std::vector<short>& waveform)
{
  //---------------------------------------------
  // Pedestal subtraction
  //---------------------------------------------
  const std::vector<short>& pedestal=m_useSignalKill ?
    frugal_pedestal_sigkill(waveform,
			    m_signalKillLookahead,
			    m_signalKillThreshold,
			    m_signalKillNContig) :
    frugal_pedestal(waveform, m_frugalNContig);
  return pedestal;
}


std::vector<short> absRunningSumTPFinderPass1::findPedestalAbsRS(const std::vector<short>& waveform)
{
  //---------------------------------------------
  //absRS Pedestal subtraction
  //---------------------------------------------
  const std::vector<short>& pedestal=m_useSignalKill ?
    frugal_pedestal_sigkill(waveform,
			    m_signalKillLookahead,
			    m_signalKillThreshold,
			    m_signalKillNContig) :
    frugal_pedestal_absRS(waveform, m_frugalNContigAbsRS);
  return pedestal;
}

std::vector<short> absRunningSumTPFinderPass1::filter(const std::vector<short>& pedsub) {

  //---------------------------------------------
  // Filtering
  //---------------------------------------------
  const size_t ntaps = m_filterTaps.size();
  const short*  taps = m_filterTaps.data();

  std::vector<short> filtered(m_doFiltering ? 
			      apply_fir_filter(pedsub, ntaps, taps) :
			      pedsub);
  if (!m_doFiltering) {
    std::transform(filtered.begin(), filtered.end(),
		   filtered.begin(), 
		   [=](short a) { return a*m_multiplier; });
  }
  return filtered;
}


std::vector<short> absRunningSumTPFinderPass1::absRS(const std::vector<short>& filtered) {

  //---------------------------------------------
  //Absolute Running Sum Algorithm
  //---------------------------------------------
  //scale all values by 10 to make sure we don't exceed numeric limits for short 
  std::vector<short> absRS_waveform(filtered.size(), 0); absRS_waveform[0] = filtered[0]/10; 

  for (size_t i = 0; i < filtered.size(); ++i){
    absRS_waveform[i] = m_R*absRS_waveform[i-1] + std::abs(filtered[i])/10; 
  }
  return absRS_waveform;
}

void
absRunningSumTPFinderPass1::hitFinding(const std::vector<short>& waveform,
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
    is_hit = adc >  (short)m_threshold;
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
absRunningSumTPFinderPass1::findHits(const std::vector<unsigned int>& channel_numbers, 
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
    hitFinding(absRSwaveform_pedsub, hits, channel_numbers[ich]);
  }
  std::cout << "Returning " << hits.size() << " hits" << std::endl;
  std::cout << "hits/channel=" << float(hits.size())/adc_samples   .size() << std::endl;
  std::cout << "hits/tick="    << float(hits.size())/adc_samples[0].size() << std::endl;
  return hits;
}

DEFINE_ART_CLASS_TOOL(absRunningSumTPFinderPass1)
