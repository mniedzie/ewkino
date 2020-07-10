#include "../interface/ttZVariables.h"

//include c++ library classes
#include <stdexcept>

//include other parts of framework
#include "../interface/ttZSelection.h"


std::map< std::string, double > ttZ::computeVariables( Event& event, const std::string& unc ){
    Met variedMet = ttZ::variedMet( event, unc );
    JetCollection variedJetCollection = ttZ::variedJetCollection( event, unc );
    PhysicsObject leptonSum = event.leptonCollection().objectSum();
    double mll;//, mtW;
    try{
        mll = event.bestZBosonCandidateMass();
//        mtW = mt( event.WLepton(), variedMet ); 
    } catch( std::domain_error& ){
        mll = ( event.lepton( 0 ) + event.lepton( 1 ) ).mass();
//        mtW = mt( event.lepton( 2 ), variedMet );
    }
    std::map< std::string, double > ret = {
        { "met", variedMet.pt() },
        { "nTI", event.isData() ? 0 : event.generatorInfo().numberOfTrueInteractions() },
        { "mll", mll },
//        { "mtW", mtW },
        { "ltmet", event.LT() + variedMet.pt() },
        { "m3l", leptonSum.mass() },
        { "mt3l", mt( leptonSum, variedMet ) },
        { "ht", variedJetCollection.scalarPtSum()},
        { "numberOfJets", variedJetCollection.size() },
//        { "numberOfBJets", variedJetCollection.numberOfTightBTaggedJets() }
        { "numberOfBJets", variedJetCollection.numberOfMediumBTaggedJets() }
    };
    return ret;
}
