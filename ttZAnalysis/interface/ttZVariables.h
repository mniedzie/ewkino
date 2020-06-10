#ifndef ttZVariables_H
#define ttZVariables_H

//include c++ library classes
#include <string>
#include <map>

//include other parts of framework
#include "../../Event/interface/Event.h"

namespace ttZ{
    std::map< std::string, double > computeVariables( Event& event, const std::string& unc );
}


#endif
