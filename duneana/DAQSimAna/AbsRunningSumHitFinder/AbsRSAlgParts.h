#ifndef ABSRSALGPARTS_H
#define ABSRSALGPARTS_H

#include <vector>

// AbsRS frugal pedestal update. Main difference --> need to increment by a value > 1. 
void do_frugal_update_absRS( short& median, int& runningDiff, const short sample,
			     const int ncontig, const int increment)
{

  if (sample>median) ++runningDiff;
  if (sample<median) --runningDiff;

  if (runningDiff > ncontig){
    median += increment; 
    runningDiff = 0;
  }

  if (runningDiff < -1*ncontig){
    median -= increment; 
    runningDiff = 0;
  }

}

std::vector<short> frugal_pedestal_absRS( const std::vector<short>& raw_in,
					  const int ncontig)
{
  short median = raw_in[0];
  std::vector<short> ped(raw_in.size(), 0);
  int increment = 50; // increment larger than ncontig, but hey it works 
  int runningDiff = 0;

  for (size_t i = 0; i < raw_in.size(); ++i){

    short s = raw_in[i]; //current sample
    do_frugal_update_absRS(median, runningDiff, s, ncontig, increment);

    ped[i] = median; 
  }

  return ped; 

}

#endif
