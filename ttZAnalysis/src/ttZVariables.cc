#include "../interface/ttZVariables.h"

//include c++ library classes
#include <stdexcept>

//include other parts of framework
#include "../interface/ttZSelection.h"


std::map< std::string, double > ttZ::computeVariables( Event& event, const std::string& unc ){
    Met variedMet = ttZ::variedMet( event, unc );
    JetCollection variedJetCollection = ttZ::variedJetCollection( event, unc );
    PhysicsObject leptonSum = event.leptonCollection().objectSum();
    double mll;
    try{
        mll = event.bestZBosonCandidateMass();
    } catch( std::domain_error& ){
        mll = ( event.lepton( 0 ) + event.lepton( 1 ) ).mass();
    }

//    Lepton Zboson = event.leptonCollection()[ event.bestZBosonCandidateIndices().first ] ;
    Lepton& lep1 =  event.leptonCollection()[ event.bestZBosonCandidateIndices().first ] ;
    Lepton& lep2 =  event.leptonCollection()[ event.bestZBosonCandidateIndices().second ] ;

    double ptZ, cosTheta;

     TVector3 Z, l;
     Z.SetPtEtaPhi( ( lep1 + lep2 ).pt(), ( lep1 + lep2 ).eta(), ( lep1 + lep2 ).phi());
     if ( lep1.charge() == -1 ) l.SetPtEtaPhi( lep1.pt(), lep1.eta(), lep1.phi());
     else l.SetPtEtaPhi( lep2.pt(), lep2.eta(), lep2.phi());

    cosTheta = Z*l / (sqrt(Z*Z) * sqrt(l*l)); 

    double gamma   = TMath::Sqrt( 1 + TMath::Power(Z.Pt()/event.bestZBosonCandidateMass() ,2) * TMath::Power(TMath::CosH(Z.Eta()),2) );

    double beta    = sqrt( 1 - 1/TMath::Power(gamma,2) );
    
    double cosThetaStar = (-beta + cosTheta) / (1 - beta*cosTheta);



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
        { "numberOfBJets", variedJetCollection.numberOfMediumBTaggedJets() },
        { "ptZ", ptZ },
        { "cosThetaStar", cosThetaStar },
    };
    return ret;
}
