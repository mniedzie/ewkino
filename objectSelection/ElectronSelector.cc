#include "../objects/interface/ElectronSelector.h"

//include c++ library classes
#include <cmath>

//include other parts of framework
#include "bTagWP.h"


/*
loose electron selection
*/


double leptonMVACutElectron(){
    return 0.0;
}


bool ElectronSelector::isLooseBase() const{
    if( electronPtr->uncorrectedPt() <= 10 ) return false;
    if( electronPtr->absEta() >= 2.5 ) return false;
    if( fabs( electronPtr->dxy() ) >= 0.05 ) return false;
    if( fabs( electronPtr->dz() ) >= 0.1 ) return false;
    if( electronPtr->sip3d() >= 8 ) return false;
    if( electronPtr->numberOfMissingHits() >= 2 ) return false;
    if( electronPtr->miniIso() >= 0.4 ) return false;
//    if( !electronPtr->passElectronMVAFall17NoIsoLoose() ) return false;
    return true;
}



bool ElectronSelector::isLoose2016() const{ 
    return true;
}


bool ElectronSelector::isLoose2017() const{
    return true;
}


bool ElectronSelector::isLoose2018() const{
    return true;
}


/*
FO electron selection
*/

bool ElectronSelector::isFOBase() const{
    if( !isLoose() ) return false;
    if( electronPtr->uncorrectedPt() < 10 ) return false;
    if( electronPtr->numberOfMissingHits() > 0 ) return false;   // added

    if( electronPtr->hOverE() >= 0.1 ) return false;
    if( electronPtr->inverseEMinusInverseP() <= -0.04 ) return false;
    if( electronPtr->etaSuperCluster() <= 1.479 ){
        if( electronPtr->sigmaIEtaEta() >= 0.011 ) return false;
    } else {
        if( electronPtr->sigmaIEtaEta() >= 0.030 ) return false;
    }
    if( electronPtr->leptonMVATOP() <= leptonMVACutElectron() ){
        if( electronPtr->ptRatio() <= 0.34 ) return false;   // changed from 0.4
        if( !electronPtr->passElectronMVAFall17NoIsoLoose() ) return false; // added, other options WP80, WP90
//        if( !electronPtr->passElectronMVAFall17NoIsoWP80() ) return false; // added, other options WP80, WP90
//        if( !electronPtr->passElectronMVAFall17NoIsoWP90() ) return false; // added, other options WP80, WP90
//        double shift = 0.0;
//        if( electronPtr->pt() < 10){
//            if( electronPtr->eta() < 0.8 ){
//                if( electronPtr->electronMVAFall17NoIso() < ( 0.700642584415 + shift ) )  return false;
//            }else if( electronPtr->eta() < 1.479 ){
//                if( electronPtr->electronMVAFall17NoIso() < ( 0.739335420875 + shift ) )  return false;
//            }
//            else{
//                if( electronPtr->electronMVAFall17NoIso() < ( 1.45390456109 + shift ) )  return false;
//            }
//        }else{
//            if( electronPtr->eta() < 0.8 ){
//                if( electronPtr->electronMVAFall17NoIso() < ( -0.146270871164 + shift ) )  return false;
//            }else if( electronPtr->eta() < 1.479 ){
//                if( electronPtr->electronMVAFall17NoIso() < ( -0.0315850882679 + shift ) )  return false;
//            }
//            else{
//                if( electronPtr->electronMVAFall17NoIso() < ( -0.0321841194737 + shift ) )  return false;
//            }
//        }
    }
    if( !electronPtr->passConversionVeto() ) return false;             // applied only to 2L selection in ttZ
//    if( !electronPtr->electronChargeConst() ) return false;            // applied only to 2L selection in ttZ
    return true;
}


bool ElectronSelector::isFO2016() const{
    if( electronPtr->leptonMVATOP() <= leptonMVACutElectron() ){
        // 0.0614-0.3093-0.7221
//        if( electronPtr->closestJetDeepFlavor() >= 0.0714 ) return false;
//        if( electronPtr->closestJetDeepFlavor() >= (0.4093+0.1) ) return false;
//        if( electronPtr->closestJetDeepFlavor() >= 0.7221 ) return false;
        if( electronPtr->closestJetDeepFlavor() >= slidingDeepFlavorThreshold( 0.3493, 0.4693, electronPtr->uncorrectedPt() ) ) return false;
    }
    return true;
}


bool ElectronSelector::isFO2017() const{
    if( electronPtr->leptonMVATOP() <= leptonMVACutElectron() ){
        // 0.0521-0.3033-07489
//        if( electronPtr->closestJetDeepFlavor() >= 0.0521 ) return false;
//        if( electronPtr->closestJetDeepFlavor() >= (0.4533+0.1) ) return false;
//        if( electronPtr->closestJetDeepFlavor() >= 0.7489 ) return false;
        if( electronPtr->closestJetDeepFlavor() >= slidingDeepFlavorThreshold( 0.3833, 0.4433, electronPtr->uncorrectedPt() ) ) return false;
    }
    return true;
}


bool ElectronSelector::isFO2018() const{
    if( electronPtr->leptonMVATOP() <= leptonMVACutElectron() ){
        // 0.0494-0.2770-0.7264
//        if( electronPtr->closestJetDeepFlavor() >= 0.0494 ) return false;
//        if( electronPtr->closestJetDeepFlavor() >= (0.3770+0.1) ) return false;
//        if( electronPtr->closestJetDeepFlavor() >= 0.7264 ) return false;
        if( electronPtr->closestJetDeepFlavor() >= slidingDeepFlavorThreshold( 0.4170, 0.3170, electronPtr->uncorrectedPt() ) ) return false;
    }
    return true;
}


/*
tight electron selection
*/

bool ElectronSelector::isTightBase() const{
    if( !isFO() ) return false;
//    if( electronPtr->leptonMVAtZq() <= leptonMVACutElectron() ) return false;
    if( electronPtr->leptonMVATOP() <= leptonMVACutElectron() ) return false;
    return true;
}


bool ElectronSelector::isTight2016() const{
    return true;
}


bool ElectronSelector::isTight2017() const{
    return true;
}


bool ElectronSelector::isTight2018() const{
    return true;
}


/*
cone correction
*/

double ElectronSelector::coneCorrection() const{
    return ( 0.6495 / electronPtr->ptRatio() );
}
