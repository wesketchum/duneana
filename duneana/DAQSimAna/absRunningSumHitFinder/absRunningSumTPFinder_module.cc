////////////////////////////////////////////////////////////////////////
// Class:       absRunningSumTPFinder
// Plugin Type: producer (art v2_10_03)
// File:        absRunningSumTPFinder_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/HitCreator.h"

#include "duneana/DAQSimAna/absRunningSumHitFinder/absRunningSumTPFinderTool.h"

#include <memory>

class absRunningSumTPFinder;


class absRunningSumTPFinder : public art::EDProducer {
public:
  explicit absRunningSumTPFinder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  absRunningSumTPFinder(absRunningSumTPFinder const &) = delete;
  absRunningSumTPFinder(absRunningSumTPFinder &&) = delete;
  absRunningSumTPFinder & operator = (absRunningSumTPFinder const &) = delete;
  absRunningSumTPFinder & operator = (absRunningSumTPFinder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
    // The module name of the raw digits we're reading in
    std::string m_inputTag;
    // The actual Service that's doing the trigger primitive finding
    std::unique_ptr<absRunningSumTPFinderTool> m_finder1;
};


absRunningSumTPFinder::absRunningSumTPFinder(fhicl::ParameterSet const & p)
  : EDProducer{p}, m_inputTag(p.get<std::string>("InputTag", "daq")),
  m_finder1{art::make_tool<absRunningSumTPFinderTool>(p.get<fhicl::ParameterSet>("finder1"))}
{
    produces<std::vector<recob::Hit>>();
    produces<art::Assns<raw::RawDigit, recob::Hit>>();
}

void absRunningSumTPFinder::produce(art::Event & e)
{
    auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);
    auto& digits_in =*digits_handle;

    art::ServiceHandle<geo::Geometry> geo;
    std::vector<std::vector<short>>  induction_samples;
    std::vector<std::vector<short>> collection_samples;
    std::vector<unsigned int>  induction_channel_numbers;
    std::vector<unsigned int> collection_channel_numbers;
    std::map<raw::ChannelID_t, const raw::RawDigit*> indChanToDigit;
    std::map<raw::ChannelID_t, const raw::RawDigit*> colChanToDigit;
    for(auto&& digit: digits_in){
        // Select just the collection channels for the primitive-finding algorithm
        const geo::SigType_t sigType = geo->SignalType(digit.Channel());
        if(sigType==geo::kInduction){
            indChanToDigit[digit.Channel()]=&digit;
            induction_channel_numbers.push_back(digit.Channel());
            induction_samples.push_back(digit.ADCs());
        }
        if(sigType==geo::kCollection){
            colChanToDigit[digit.Channel()]=&digit;
            collection_channel_numbers.push_back(digit.Channel());
            collection_samples.push_back(digit.ADCs());
        }
    }

    // Pass the full list of collection and induction channels to the hit finding algorithm
    std::vector<absRunningSumTPFinderTool::Hit> hits1=m_finder1->findHits(induction_channel_numbers,  induction_samples);
    std::vector<absRunningSumTPFinderTool::Hit> hits2=m_finder1->findHits(collection_channel_numbers, collection_samples);

    // Loop over the returned induction trigger primitives and turn them into recob::Hits
    recob::HitCollectionCreator hcol(e, false /* doWireAssns */, true /* doRawDigitAssns */);
    for(auto const& hit : hits1){
        const raw::RawDigit* digit=indChanToDigit[hit.channel];
        if(!digit){
            std::cout << "No digit with channel " << hit.channel << " found. Did you set the channel correctly?" << std::endl;
        }
        std::vector<geo::WireID> wids = geo->ChannelToWire(hit.channel);
        geo::WireID wid = wids[0];

        recob::HitCreator lar_hit(*digit,                           //RAW DIGIT REFERENCE.
                              wid,                                  //WIRE ID.
                              hit.startTime,                        //START TICK.
                              hit.startTime+hit.timeOverThreshold,  //END TICK. 
                              hit.timeOverThreshold,                //RMS.
                              hit.startTime,                        //PEAK_TIME.
                              0,                                    //SIGMA_PEAK_TIME.
                              0,                                    //PEAK_AMPLITUDE.
                              0,                                    //SIGMA_PEAK_AMPLITUDE.
                              hit.charge,                           //HIT_INTEGRAL.
                              0,                                    //HIT_SIGMA_INTEGRAL.
                              hit.charge,                           //SUMMED CHARGE. 
                              0,                                    //MULTIPLICITY.
                              0,                                    //LOCAL_INDEX.
                              0,                                    //WIRE ID.
                              0                                     //DEGREES OF FREEDOM.
            );
        hcol.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits_handle, 0});
    }
    // Loop over the returned collection trigger primitives and turn them into recob::Hits
    for(auto const& hit : hits2){
        const raw::RawDigit* digit=colChanToDigit[hit.channel];
        if(!digit){
            std::cout << "No digit with channel " << hit.channel << " found. Did you set the channel correctly?" << std::endl;
        }
        std::vector<geo::WireID> wids = geo->ChannelToWire(hit.channel);
        geo::WireID wid = wids[0];

        recob::HitCreator lar_hit(*digit,                           //RAW DIGIT REFERENCE.
                              wid,                                  //WIRE ID.
                              hit.startTime,                        //START TICK.
                              hit.startTime+hit.timeOverThreshold,  //END TICK. 
                              hit.timeOverThreshold,                //RMS.
                              hit.startTime,                        //PEAK_TIME.
                              0,                                    //SIGMA_PEAK_TIME.
                              0,                                    //PEAK_AMPLITUDE.
                              0,                                    //SIGMA_PEAK_AMPLITUDE.
                              hit.charge,                           //HIT_INTEGRAL.
                              0,                                    //HIT_SIGMA_INTEGRAL.
                              hit.charge,                           //SUMMED CHARGE. 
                              0,                                    //MULTIPLICITY.
                              0,                                    //LOCAL_INDEX.
                              0,                                    //WIRE ID.
                              0                                     //DEGREES OF FREEDOM.
            );
        hcol.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits_handle, 0});
    }
    hcol.put_into(e);
}

DEFINE_ART_MODULE(absRunningSumTPFinder)
