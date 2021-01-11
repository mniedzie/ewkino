#ifndef ttZSelection_H
#define ttZSelection_H

//include ROOT classes
#include "TH2.h"

//include other parts of framework
#include "../../Event/interface/Event.h"

namespace ttZ{
    void applyBaselineObjectSelection( Event& event, const bool allowUncertainties = false );
    bool passLowMllVeto( const Event& event, const double vetoValue = 12. );
    bool passBaselineSelection( Event& event, const bool allowUncertainties = false, const bool mllVeto = true );
    JetCollection variedJetCollection( const Event& event, const std::string& uncertainty );
    JetCollection::size_type numberOfVariedBJets( const Event& event, const std::string& uncertainty );
    JetCollection::size_type numberOfVariedJets( const Event& event, const std::string& uncertainty );
    Met variedMet( const Event& event, const std::string& uncertainty );
    bool passVariedSelection( const Event& event, const std::string& uncertainty );
    bool passVariedSelectionWZCR( Event& event, const std::string& uncertainty );
    bool passVariedSelectionTTZCR( Event& event, const std::string& uncertainty );

    bool passSelectionLNumber( Event& event );
    bool passSelectionTTZ( Event& event, const std::string& uncertainty );
    bool passSelectionTTZclean( Event& event, const std::string& uncertainty );
    bool passSelectionWZCR( Event& event, const std::string& uncertainty );
    bool passSelectionDYCR( Event& event, const std::string& uncertainty );
    bool passSelectionttbarCR( Event& event, const std::string& uncertainty );
    bool passSelectionZZCR( Event& event, const std::string& uncertainty );

    bool passVariedSelectionNPCR( Event& event, const std::string& uncertainty );
    bool passVariedSelectionXGammaCR( Event& event, const std::string& uncertainty );
    bool passTriggerSelection( const Event& event );
    bool passPtCuts( const Event& event );
    bool leptonsArePrompt( const Event& event );
    bool leptonsAreTight( const Event& event );
    double fakeRateWeight( const Event& event, const std::shared_ptr< TH2 >& muonMap, const std::shared_ptr< TH2 >& electronMap );
    bool passPhotonOverlapRemoval( const Event& event );

    unsigned SR_main( const Event& event, const int nL, const int nJ, const int nB);
    unsigned ttZFlavPlot( const Event& event );
    unsigned ttZFlavPlot3L( const Event& event );
    unsigned ttZFlavPlot4L( const Event& event );
}

#endif
