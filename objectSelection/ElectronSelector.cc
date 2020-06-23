#include "../objects/interface/ElectronSelector.h"

//include c++ library classes
#include <cmath>

//include other parts of framework
#include "bTagWP.h"


/*
loose electron selection
*/


double leptonMVACutElectron(){
    return 0.4;
}


bool ElectronSelector::isLooseBase() const{
    if( electronPtr->uncorrectedPt() <= 10 ) return false;
    if( electronPtr->absEta() >= 2.5 ) return false;
    if( fabs( electronPtr->dxy() ) >= 0.05 ) return false;
    if( fabs( electronPtr->dz() ) >= 0.1 ) return false;
    if( electronPtr->sip3d() >= 8 ) return false;
    if( electronPtr->numberOfMissingHits() >= 2 ) return false;
    if( electronPtr->miniIso() >= 0.4 ) return false;
    if( !electronPtr->passElectronMVAFall17NoIsoLoose() ) return false;
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

//if( event.numberOfLeptons() < 3 ) return false;

//    if( electronPtr->hOverE() >= 0.1 ) return false;
//    if( electronPtr->inverseEMinusInverseP() <= -0.04 ) return false;
//    if( electronPtr->etaSuperCluster() <= 1.479 ){
//        if( electronPtr->sigmaIEtaEta() >= 0.011 ) return false;
//    } else {
//        if( electronPtr->sigmaIEtaEta() >= 0.030 ) return false;
//    }
    if( electronPtr->leptonMVAtZq() <= leptonMVACutElectron() ){
        if( electronPtr->ptRatio() <= 0.4 ) return false;
    }
//    if( electronPtr->leptonMVAtZq() <= leptonMVACutElectron() ){
//        if( !electronPtr->passElectronMVAFall17NoIsoWP80() ) return false;
//        if( electronPtr->ptRatio() <= 0.7 ) return false;
//    }
//    if( electronPtr->numberOfMissingHits() > 0 ) return false;         // applied only to 2L selection in ttZ
//    if( !electronPtr->passConversionVeto() ) return false;             // applied only to 2L selection in ttZ
//    if( !electronPtr->electronChargeConst() ) return false;            // applied only to 2L selection in ttZ
    return true;
}


bool ElectronSelector::isFO2016() const{
//    if( electronPtr->closestJetDeepFlavor() >= bTagWP::mediumDeepFlavor2016() ) return false;
    double electronMVAvalue = electronPtr->electronMVASummer16GP();
    if( electronPtr->closestJetDeepCSV() > bTagWP::tightDeepCSV2016() ) return false;
	if( electronMVAvalue < ( -0.1 + (electronPtr->absEta() >= 1.479)*0.8 ) ) return false;
    if( electronPtr->leptonMVAtZq() <= leptonMVACutElectron() ){
        if( electronPtr->closestJetDeepCSV() > 0.4 ) return false;
    }
    return true;
}


bool ElectronSelector::isFO2017() const{
//    if( electronPtr->closestJetDeepFlavor() >= bTagWP::mediumDeepFlavor2017() ) return false;
    double electronMVAvalue = electronPtr->electronMVAFall17NoIso();
    if( electronPtr->closestJetDeepCSV() > bTagWP::tightDeepCSV2017() ) return false;
	if( electronMVAvalue < ( -0.3 + (electronPtr->absEta() >= 1.479)*0.6 ) ) return false;
    if( electronPtr->leptonMVAtZq() <= leptonMVACutElectron() ){
        if( electronPtr->closestJetDeepCSV() > 0.4 ) return false;
    }
    return true;
}


bool ElectronSelector::isFO2018() const{
//    if( electronPtr->closestJetDeepFlavor() >= bTagWP::mediumDeepFlavor2018() ) return false;
    double electronMVAvalue = electronPtr->electronMVAFall17NoIso();
    if( electronPtr->closestJetDeepCSV() > bTagWP::tightDeepCSV2018() ) return false;
	if( electronMVAvalue < ( -0.3 + (electronPtr->absEta() >= 1.479)*0.6 ) ) return false;
    if( electronPtr->leptonMVAtZq() <= leptonMVACutElectron() ){
        if( electronPtr->closestJetDeepCSV() > 0.4 ) return false;
    }
    return true;
}


/*
tight electron selection
*/

bool ElectronSelector::isTightBase() const{
    if( !isFO() ) return false;
    if( electronPtr->leptonMVAtZq() <= leptonMVACutElectron() ) return false;
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
    return ( 0.9 / electronPtr->ptRatio() );
}
