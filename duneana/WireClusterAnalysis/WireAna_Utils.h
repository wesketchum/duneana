#ifndef WIREANA_UTILS_H
#define WIREANA_UTILS_H

#include <vector>
#include <cmath>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

namespace wireana{
  using std::string;
  using std::vector;
  using std::pair;
  using std::map;
  using std::set;

  using art::Ptr;
  using art::Event;
  using art::ServiceHandle;

  using cheat::BackTrackerService;
  using cheat::ParticleInventoryService;

  using simb::MCParticle;
  using sim::SimChannel;

  using detinfo::DetectorClocksData;

  using recob::Wire;
  using recob::Hit;

  enum PType{ kUnknown=0, kMarl, kAPA, kCPA, kAr39, kAr42, kNeutron, kKryp, kPlon, kRdon, kNPTypes };
  enum EType{
    kES,
    kCC,
    kRad
  };

  

  struct wirecluster{
    double totalCharge=0.;
    double totalAbsCharge=0.;
    int nWires=0;
    std::vector<art::Ptr<recob::Wire> > wire_pointers;
    double avgCharge(){ return totalCharge/nWires; }
    double avgAbsCharge() {return totalAbsCharge/nWires; }
    void push_back( art::Ptr<recob::Wire> &p )
    {
      this->wire_pointers.push_back( p );
      this->nWires++;
    }
    art::Ptr<recob::Wire> back()
    {
      return this->wire_pointers.back();
    }
    
  };

  struct roi{
    roi( art::Ptr<recob::Wire> &wire, recob::Wire::RegionsOfInterest_t::datarange_t r)
    {
      updateroi(wire, r);
    }
    roi(){};

    void updateroi( art::Ptr<recob::Wire> &wire, recob::Wire::RegionsOfInterest_t::datarange_t r)
    {
      this->range = r;
      this->channel = wire->Channel();
      this->wire = wire;
      this->begin_index = r.begin_index();
      this->end_index = r.end_index();
      this->width = this->end_index - this->begin_index;
      this->sum = this->Sum();
      this->abs_sum = this->Sum(true);
      this->centroid = this->Centroid();
      this->abs_centroid = this->Centroid(true);
    }

    art::Ptr<recob::Wire> wire;
    recob::Wire::RegionsOfInterest_t::datarange_t range;
    int channel=-1;
    int begin_index=-1;
    int end_index=-1;
    int width=-1;
    int clusterID=-1;
    double sum=0;
    double abs_sum=0;
    double centroid=0;
    double abs_centroid=0;
    double Sum(bool useabs=false)
    {
      double ret=0;
      for( auto it = this->range.begin(); it!=this->range.end(); ++it ) ret+= (useabs)? abs(*it):*it;
      return ret;
    }
    double Centroid(bool useabs=false)
    {
      double denum = 0;
      double num = 0;
      for( size_t i = this->range.begin_index(); i != this->range.end_index(); i++ )
      {
        double v = this->range[i];
        if (useabs) v = abs(v);
        denum+=v;
        num+=(i*v);
      }
      if(denum!=0) return num/denum;
      else return -999;
    }
    std::string label;
    std::vector< std::pair<int, double> > pdg_energy_list;
  };



  struct roicluster{
    int nWires;
    std::string mainlabel;
    std::vector< std::pair<int, double> > pdg_energy_list;
    std::vector< roi > ROIs;

  };


  typedef std::map<int, std::map< geo::View_t, std::vector<roi> > > PlaneViewROIMap;

  struct ApaROIConstainer
  {
    int PlaneID;
    std::map<int, std::vector<roi>> ViewROIMap;
  };


  class WireAnaDBSCAN;
}


class wireana::WireAnaDBSCAN {
public:    
    WireAnaDBSCAN(){
    }
    ~WireAnaDBSCAN(){}

    void SetParameters(std::vector<wireana::roi> &ROIs, unsigned int minPts=2, float eps=4.5, float drift=1.6, float pitch=3.0 ){
        m_minPoints = minPts;
        m_epsilon = eps;
        m_ROIs = &ROIs;
        m_roiSize = ROIs.size();
        m_drift = drift;
        m_pitch = pitch;
    }

    int run();
    void promoteSingleROI( int minWidth=150, double minAbsSum = 80 );

    vector<int> calculateCluster(wireana::roi &r);
    int expandCluster(wireana::roi &r, int clusterID);

    float calculateTickDistance( const wireana::roi& r1, const wireana::roi&r2 );
    float calculateDistance( const wireana::roi &r1, const wireana::roi &r2, float driftspeed=1.6, float pitch=3.0 );

    void setDriftAndPitch( float drift, float pitch ){ m_drift=drift; m_pitch = pitch ; }
    int getTotalROISize() {return m_roiSize;}
    int getMinimumClusterSize() {return m_minPoints;}
    int getEpsilonSize() {return m_epsilon;}

    std::vector<wireana::roi> &GetROIs(){ return *m_ROIs; }

    std::vector<wireana::roi> *m_ROIs;
    
private:    
    unsigned int m_roiSize;
    unsigned int m_minPoints;
    float m_epsilon;
    float m_drift;
    float m_pitch;
};




int 
wireana::WireAnaDBSCAN::run()
{
  int clusterID = 1;
  std::vector<roi>::iterator iter;
  if (!m_ROIs) return -1;
  for(iter = m_ROIs->begin(); iter != m_ROIs->end(); ++iter)
  {
    if ( iter->clusterID == UNCLASSIFIED )
    {
      if ( expandCluster(*iter, clusterID) != FAILURE )
      {
        clusterID += 1;
      }
    }
  }
  return 0;
}


std::vector<int> 
wireana::WireAnaDBSCAN::calculateCluster(wireana::roi &r)
{
  int index = 0;
  std::vector<wireana::roi>::iterator iter;
  std::vector<int> clusterIndex;
  for( iter = m_ROIs->begin(); iter != m_ROIs->end(); ++iter)
  {
    if ( calculateDistance(r, *iter, m_drift, m_pitch) <= m_epsilon )
    {
      clusterIndex.push_back(index);
    }
    index++;
  }
  return clusterIndex;
}

int 
wireana::WireAnaDBSCAN::expandCluster(wireana::roi &r, int clusterID)
{  
  std::vector<int> clusterSeeds = calculateCluster(r);

  if ( clusterSeeds.size() < m_minPoints )
  {
    r.clusterID = NOISE;
    return FAILURE;
  }
  else
  {
    int index = 0, indexCorePoint = 0;
    std::vector<int>::iterator iterSeeds;
    for( iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds)
    {
      m_ROIs->at(*iterSeeds).clusterID = clusterID;
      if (m_ROIs->at(*iterSeeds).channel == r.channel && m_ROIs->at(*iterSeeds).begin_index == r.begin_index && m_ROIs->at(*iterSeeds).end_index == r.end_index )
      {
        indexCorePoint = index;
      }
      ++index;
    }
    clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);

    for( std::vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i )
    {
      std::vector<int> clusterNeighors = calculateCluster(m_ROIs->at(clusterSeeds[i]));

      if ( clusterNeighors.size() >= m_minPoints )
      {
        std::vector<int>::iterator iterNeighors;
        for ( iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors )
        {
          if ( m_ROIs->at(*iterNeighors).clusterID == UNCLASSIFIED || m_ROIs->at(*iterNeighors).clusterID == NOISE )
          {
            if ( m_ROIs->at(*iterNeighors).clusterID == UNCLASSIFIED )
            {
              clusterSeeds.push_back(*iterNeighors);
              n = clusterSeeds.size();
            }
            m_ROIs->at(*iterNeighors).clusterID = clusterID;
          }
        }
      }
    }
    return SUCCESS;
  }
}

float 
wireana::WireAnaDBSCAN::calculateTickDistance(const wireana::roi &r1, const wireana::roi &r2)
  {
    roi rr1 = (r1.begin_index < r2.begin_index)? r1 : r2;
    roi rr2 = (r1.begin_index < r2.begin_index)? r2 : r1;
    float ret = rr2.begin_index - rr1.end_index;
    return (ret>0)? ret: 0;
  }
float 
wireana::WireAnaDBSCAN::calculateDistance( const roi &r1, const roi &r2, float driftspeed, float pitch )
{
  float dx = this->calculateTickDistance(r1,r2)*driftspeed/2.;
  float dz = (r1.channel - r2.channel)*pitch;
  return pow(dx*dx+dz*dz,0.5);
}



#endif
