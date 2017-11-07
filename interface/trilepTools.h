#ifndef trilepTools_H
#define trilepTools_H
//include c++ library utilities
#include <utility>
//include ROOT utilities
#include "TLorentzVector.h"
namespace trilep{
    std::pair<unsigned, unsigned> bestZ(const TLorentzVector*, const unsigned); //return indices of leptons forming best Z-candidate
    unsigned flavorChargeComb(const std::vector<unsigned>& ind, const unsigned*, const int*, const unsigned); //Determine whether event contains OSSF or OSOF pair, or not 
}
#endif
