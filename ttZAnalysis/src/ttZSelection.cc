#include "../interface/ttZSelection.h"

//include c++ library classes
#include <stdexcept>

//include other parts of framework
#include "../../Tools/interface/histogramTools.h"
#include "../../Tools/interface/stringTools.h"
#include "../../constants/particleMasses.h"



void ttZ::applyBaselineObjectSelection( Event& event, const bool allowUncertainties ){

    event.selectLooseLeptons();
    event.removeTaus();
    event.cleanElectronsFromLooseMuons();
//    event.cleanTausFromLooseLightLeptons();
    event.cleanJetsFromFOLeptons();
    if( allowUncertainties ){
        event.jetCollection().selectGoodAnyVariationJets();
    } else {
        event.selectGoodJets();
    }
}


bool ttZ::passLowMllVeto( const Event& event, const double vetoValue ){

    for( const auto& leptonPtrPair : event.lightLeptonCollection().pairCollection() ){

        //veto OSSF pairs of low mass
        Lepton& lepton1 = *( leptonPtrPair.first );
        Lepton& lepton2 = *( leptonPtrPair.second );
        if( !oppositeSignSameFlavor( lepton1, lepton2 ) ){
            continue;
        }
        if( ( lepton1 + lepton2 ).mass() < vetoValue ){
            return false;
        }
    }
    return true;
}


bool ttZ::passBaselineSelection( Event& event, const bool allowUncertainties, const bool mllVeto ){
    
	// only loose leptons are selected
//            if( !event.isData() ) std::cout << "mc in baseline selection funciton" <<  std::endl;
    applyBaselineObjectSelection( event, allowUncertainties );
	// accept event with at least 3 loose leptons
    if( event.numberOfLeptons() < 3 ) return false;
//            if( !event.isData() ) std::cout << "mc with 3 loose leptons " <<  std::endl;
    if( mllVeto && !passLowMllVeto( event, 12 ) ) return false;
//            if( !event.isData() ) std::cout << "mc with no light mll" <<  std::endl;

    event.selectFOLeptons();
    if( event.numberOfLeptons() < 3 ) return false;
//            if( !event.isData() ) std::cout << "mc with 3 FO leptons" <<  std::endl;
    event.applyLeptonConeCorrection();
    event.sortLeptonsByPt();
    return true;
}


JetCollection ttZ::variedJetCollection( const Event& event, const std::string& uncertainty ){
    if( uncertainty == "nominal" ){
        return event.jetCollection().goodJetCollection();
    } else if( uncertainty == "JECDown" ){
        return event.jetCollection().JECDownCollection().goodJetCollection();
    } else if( uncertainty == "JECUp" ){
        return event.jetCollection().JECUpCollection().goodJetCollection();
    } else if( uncertainty == "JERDown" ){
        return event.jetCollection().JERDownCollection().goodJetCollection();
    } else if( uncertainty == "JERUp" ){
        return event.jetCollection().JERUpCollection().goodJetCollection();
    } else if( uncertainty == "UnclDown" ){
        return event.jetCollection().goodJetCollection();
    } else if( uncertainty == "UnclUp" ){
        return event.jetCollection().goodJetCollection();
    } else {
        throw std::invalid_argument( "Uncertainty source " + uncertainty + " is unknown." );
    }
}


JetCollection::size_type ttZ::numberOfVariedBJets( const Event& event, const std::string& uncertainty ){
    return ttZ::variedJetCollection( event, uncertainty ).numberOfTightBTaggedJets();
}

JetCollection::size_type ttZ::numberOfVariedJets( const Event& event, const std::string& uncertainty ){
    return ttZ::variedJetCollection( event, uncertainty ).numberOfGoodAnyVariationJets();
}


Met ttZ::variedMet( const Event& event, const std::string& uncertainty ){
    if( uncertainty == "nominal" ){
        return event.met();
    } else if( uncertainty == "JECDown" ){
        return event.met().MetJECDown();
    } else if( uncertainty == "JECUp" ){
        return event.met().MetJECUp();
    } else if( uncertainty == "JERDown" ){
        return event.met();
    } else if( uncertainty == "JERUp" ){
        return event.met();
    } else if( uncertainty == "UnclDown" ){
        return event.met().MetUnclusteredDown();
    } else if( uncertainty == "UnclUp" ){
        return event.met().MetUnclusteredUp();
    } else {
        throw std::invalid_argument( "Uncertainty source " + uncertainty + " is unknown." );
    }
}


bool ttZ::passVariedSelection( const Event& event, const std::string& uncertainty ){
    static constexpr double metCut = 50;
    if( numberOfVariedBJets( event, uncertainty ) > 0 ) return false;
    if( variedMet( event, uncertainty ).pt() < metCut ) return false;
    return true;
}


bool ttZ::passVariedSelectionWZCR( Event& event, const std::string& uncertainty ){
    static constexpr double minMet = 30;
    static constexpr double maxMet = 100;
    static constexpr double minMT = 50;
    static constexpr double maxMT = 100;
    if( numberOfVariedBJets( event, uncertainty ) > 0 ) return false;
    Met met = variedMet( event, uncertainty );
    if( met.pt() < minMet || met.pt() > maxMet ) return false;
    if( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 15 ) return false;
    double mTW = mt( met, event.WLepton() );
    if( mTW < minMT || mTW > maxMT ) return false;
    return true;
}


bool ttZ::passVariedSelectionTTZCR( Event& event, const std::string& uncertainty ){
    static constexpr size_t numberOfBJets = 1;
    if( numberOfVariedBJets( event, uncertainty ) < numberOfBJets ) return false;
    if( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 15 ) return false;
    if( std::abs( event.leptonCollection().objectSum().mass() - particle::mZ ) < 15 ) return false;
    return true;
}

bool ttZ::passSelectionTTZ( Event& event, const std::string& uncertainty ){
    if ( numberOfVariedJets( event, uncertainty ) < 1 ) return false;
    if ( !( event.hasOSSFLightLeptonPair() ) ) return false;
    if ( event.numberOfFOLeptons() == 3 && ( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 10 ) ) return false;
    if ( event.numberOfFOLeptons() == 4 && ( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 20 ) ) return false;
    if ( event.numberOfFOLeptons() < 3 ) return false;
    if ( event.numberOfFOLeptons() == 3 && numberOfVariedBJets( event, uncertainty ) > 0 && numberOfVariedJets( event, uncertainty ) < 2  ) return false;
    if ( event.numberOfFOLeptons() == 4 && numberOfVariedJets( event, uncertainty ) < 2 ) return false;
    return true;
}


bool ttZ::passSelectionTTZNP( Event& event, const std::string& uncertainty ){
    if ( numberOfVariedJets( event, uncertainty ) < 1 ) return false;
    if ( !( event.hasOSSFLightLeptonPair() ) ) return false;
    if ( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 10 ) return false;
    if ( event.numberOfFOLeptons() == 3 && ( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 10 ) ) return false;
    if ( event.numberOfFOLeptons() == 4 && ( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 20 ) ) return false;
    if ( event.numberOfFOLeptons() < 3 ) return false;
    if ( event.numberOfFOLeptons() == 3 && numberOfVariedBJets( event, uncertainty ) > 0 && numberOfVariedJets( event, uncertainty ) < 2  ) return false;
    if ( event.numberOfFOLeptons() == 4 && numberOfVariedJets( event, uncertainty ) < 2 ) return false;
    return true;
}

bool ttZ::passVariedSelectionNPCR( Event& event, const std::string& uncertainty ){
    static constexpr size_t numberOfBJets = 1;
    if( numberOfVariedBJets( event, uncertainty ) < numberOfBJets ) return false;
    return true;
}


bool ttZ::passVariedSelectionXGammaCR( Event& event, const std::string& uncertainty ){
    if( numberOfVariedBJets( event, uncertainty ) > 0 ) return false;
    if( variedMet( event, uncertainty ).pt() >= 50 ) return false;
    if( event.bestZBosonCandidateMass() >= 75 ) return false;
    if( std::abs( event.leptonCollection().objectSum().mass() - particle::mZ ) >= 15 ) return false;
    return true;
}


bool ttZ::passTriggerSelection( const Event& event ){
//    if( event.numberOfMuons() >= 1 ){
//        if( event.passTriggers_m() ) return true;
//    } 
//    if( event.numberOfMuons() >= 2 ){
//        if( event.passTriggers_mm() ) return true;
//    }
//    if( event.numberOfMuons() >= 3 ){
//        if( event.passTriggers_mmm() ) return true;
//    }
//    if( event.numberOfElectrons() >= 1 ){
//        if( event.passTriggers_e() ) return true;
//    }
//    if( event.numberOfElectrons() >= 2 ){
//        if( event.passTriggers_ee() ) return true;
//    }
//    if( event.numberOfElectrons() >= 3 ){
//        if( event.passTriggers_eee() ) return true;
//    }
//    if( ( event.numberOfMuons() >= 1 ) && ( event.numberOfElectrons() >= 1 ) ){
//        if( event.passTriggers_em() ) return true;
//    }
//    if( ( event.numberOfMuons() >= 2 ) && ( event.numberOfElectrons() >= 1 ) ){
//        if( event.passTriggers_emm() ) return true;
//    }
//    if( ( event.numberOfMuons() >= 1 ) && ( event.numberOfElectrons() >= 2 ) ){
//        if( event.passTriggers_eem() ) return true;
//    }

        if( event.passTriggers_m() || 
		    event.passTriggers_mm() ||
            event.passTriggers_mmm() ||
            event.passTriggers_e() ||
            event.passTriggers_ee() ||
            event.passTriggers_eee() ||
            event.passTriggers_em() ||
            event.passTriggers_emm() ||
            event.passTriggers_eem() ) return true;

    return false;
}


bool ttZ::passPtCuts( const Event& event ){
    
    //assume leptons were ordered while applying baseline selection
    
    //leading lepton
    if( event.lepton( 0 ).isMuon() && event.lepton( 0 ).pt() <= 40 ) return false;
    if( event.lepton( 0 ).isElectron() && event.lepton( 0 ).pt() <= 40 ) return false;

    //subleading lepton
    if( event.lepton( 1 ).isMuon() && event.lepton( 1 ).pt() <= 20 ) return false;
    if( event.lepton( 1 ).isElectron() && event.lepton( 1 ).pt() <= 20 ) return false;

    //trailing lepton
    if( event.lepton( 2 ).isMuon() && event.lepton( 2 ).pt() <= 10 ) return false;
    if( event.lepton( 2 ).isElectron() && event.lepton( 2 ).pt() <= 10 ) return false;
    return true;
}


bool ttZ::leptonsArePrompt( const Event& event ){
    for( const auto& leptonPtr : event.leptonCollection() ){
        if( leptonPtr->isFO() && !leptonPtr->isPrompt() ) return false;
    }
    return true;
}


bool ttZ::leptonsAreTight( const Event& event ){
    for( const auto& leptonPtr : event.leptonCollection() ){
        if( leptonPtr->isFO() && !leptonPtr->isTight() ) return false;
    }
    return true;
}


bool leptonFromMEExternalConversion( const Lepton& lepton ){
    if( !( lepton.matchPdgId() == 22 ) ) return false;
    if( !( lepton.isPrompt() && lepton.provenanceConversion() == 0 ) ) return false;
    return true;
}


bool ttZ::passPhotonOverlapRemoval( const Event& event ){
    bool isPhotonSample = false;
    bool isInclusiveSample = false;
    std::string sampleName = event.sample().fileName();
    if( stringTools::stringContains( sampleName, "DYJetsToLL" ) || stringTools::stringContains( sampleName, "TTTo" ) || stringTools::stringContains( sampleName, "TTJets" ) ){
        isInclusiveSample = true;
    } else if( stringTools::stringContains( sampleName, "TTGamma" ) || stringTools::stringContains( sampleName, "ZGToLLG" ) || stringTools::stringContains( sampleName, "WGToLNuG" ) ){
        isPhotonSample = true;
    }

    if( !( isPhotonSample || isInclusiveSample ) ){
        return true;
    }

    bool hasMEExternalConversion = false;
    for( const auto& leptonPtr : event.leptonCollection() ){
        if( leptonPtr->isFO() && leptonFromMEExternalConversion( *leptonPtr ) ){
            hasMEExternalConversion = true;
            break;
        }
    }
    if( isInclusiveSample ){
        return !hasMEExternalConversion;
    } else if( isPhotonSample ){
        return hasMEExternalConversion;
    }
    return true;
}


double ttZ::fakeRateWeight( const Event& event, const std::shared_ptr< TH2 >& muonMap, const std::shared_ptr< TH2 >& electronMap ){
    static constexpr double maxPt = 44.;

    double weight = -1.;
    for( const auto& leptonPtr : event.leptonCollection() ){
        if( !leptonPtr->isFO() ) continue;
        if( leptonPtr->isTight() ) continue;
        double fr;
        double pt = std::min( leptonPtr->pt(), maxPt );
        if( leptonPtr->isMuon() ){
            fr = histogram::contentAtValues( muonMap.get(), pt, leptonPtr->absEta() );
        } else if( leptonPtr->isElectron() ){
            fr = histogram::contentAtValues( electronMap.get(), pt, leptonPtr->absEta() );
        } else {
            throw std::invalid_argument( "we are not considering taus for now" );
        }
        weight *= - fr / ( 1. - fr );
    }
    return weight;
}



unsigned ttZ::SR_main(const int nL, const int nJ, const int nB){ //3 light leptons, OSSF
    unsigned sr = 0;
    if( nL == 3){
        if( nB == 0 ){
            if ( nJ == 2 ) sr += 1;
            else if ( nJ == 3 ) sr+= 2;
            else if ( nJ > 3 ) sr+= 3;
        }
        else if( nB == 1 ){
            sr += 4;
            if ( nJ == 3 ) sr += 1;
            else if ( nJ == 4 ) sr+= 2;
            else if ( nJ > 4 ) sr+= 3;
        }
        else{
            sr += 8;
            if ( nJ == 3 ) sr += 1;
            else if ( nJ == 4 ) sr+= 2;
            else if ( nJ > 4 ) sr+= 3;
        }
    }
    else if ( nL == 4 ){
        sr += 12;
        if ( nB > 0 ) sr += 1;
    }
    return sr;
}



unsigned ttZ::ttZFlavPlot( const Event& event ){
    if( event.numberOfLightLeptons() == 3 ){
        if( event.numberOfMuons() == 3 ) return 0;
        if( event.numberOfMuons() == 2 ) return 1;
        if( event.numberOfMuons() == 1 ) return 2;
        if( event.numberOfElectrons() == 3 ) return 3;
    }
    if( event.numberOfLightLeptons() == 4 ){
        if( event.numberOfMuons() == 4 ) return 0;
        if( event.numberOfMuons() == 3 ) return 1;
        if( event.numberOfElectrons() == 4 ) return 3;
        if( event.numberOfElectrons() > 2 ) return 2;
    }
    return 0;
}
