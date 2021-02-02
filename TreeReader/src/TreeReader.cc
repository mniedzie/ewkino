#include "../interface/TreeReader.h"

//include c++ library classes 
#include <fstream>
#include <iostream>
#include <typeinfo>

//include other parts of analysis framework
#include "../../Tools/interface/analysisTools.h"
#include "../../Tools/interface/stringTools.h"
#include "../../Tools/interface/systemTools.h"
#include "../../Tools/interface/analysisTools.h"
#include "../../Event/interface/Event.h"
#include "../../constants/luminosities.h"


TreeReader::TreeReader( const std::string& sampleListFile, const std::string& sampleDirectory ){
    readSamples( sampleListFile, sampleDirectory );
}


void TreeReader::readSamples( const std::string& list, const std::string& directory, std::vector<Sample>& sampleVector ){

    //clean current sample list 
    sampleVector.clear();

    //read list of samples from file
    sampleVector = readSampleList(list, directory);

    //print sample information
    for(auto& sample : sampleVector){
        std::cout << sample << std::endl;
    }
}


void TreeReader::readSamples( const std::string& list, const std::string& directory ){
    readSamples( list, directory, this->samples );
}


void TreeReader::readSamples2016( const std::string& list, const std::string& directory ){
    std::cout << "########################################" << std::endl;
    std::cout << "         2016 samples                   " << std::endl;
    std::cout << "########################################" << std::endl;

    readSamples( list, directory, this->samples2016 );

    //add the 2016 samples to the total sample list 
    this->samples.insert( samples.end(), samples2016.begin(), samples2016.end() );

    //check for errors
    checkSampleEraConsistency();

}


void TreeReader::readSamples2017( const std::string& list, const std::string& directory ){
    std::cout << "########################################" << std::endl;
    std::cout << "         2017 samples                   " << std::endl;
    std::cout << "########################################" << std::endl;

    readSamples( list, directory, this->samples2017 );

    //add the 2017 samples to the total sample list
    this->samples.insert( samples.end(), samples2017.begin(), samples2017.end() );

    //check for errors 
    checkSampleEraConsistency();
}


std::pair< std::map< std::string, bool >, std::map< std::string, TBranch* > > buildBranchMap( TTree* treePtr, const std::string& nameIdentifier, const std::string& antiIdentifier = "" ){   
    std::map< std::string, bool > decisionMap;
    std::map< std::string, TBranch* > branchMap;
    TObjArray* branch_list = treePtr->GetListOfBranches();
    for( const auto& branchPtr : *branch_list ){
        std::string branchName = branchPtr->GetName();
        if( stringTools::stringContains( branchName, nameIdentifier ) ){
            if( antiIdentifier != "" && stringTools::stringContains( branchName, antiIdentifier ) ) continue;
            decisionMap[ branchName ] = false;
            branchMap[ branchName ] = nullptr;
        }
    }
    return { decisionMap, branchMap };
}


void TreeReader::initializeTriggerMap( TTree* treePtr ){
    auto triggerMaps = buildBranchMap( treePtr, "HLT", "prescale" );
    _triggerMap = triggerMaps.first;
    b__triggerMap = triggerMaps.second;
}


void TreeReader::initializeMetFilterMap( TTree* treePtr ){

    //WARNING: Currently one MET filter contains 'updated' rather than 'Flag' in the name. If this changes, make sure to modify the code here!
    auto filterMaps = buildBranchMap( treePtr, "Flag" );
    auto filterMaps_part2 = buildBranchMap( treePtr , "updated" );
    filterMaps.first.insert( filterMaps_part2.first.cbegin(), filterMaps_part2.first.cend() );
    filterMaps.second.insert( filterMaps_part2.second.cbegin(), filterMaps_part2.second.cend() );

    _MetFilterMap = filterMaps.first;
    b__MetFilterMap = filterMaps.second;
}


bool treeHasBranchWithName( TTree* treePtr, const std::string& nameToFind ){
	TObjArray* branch_list = treePtr->GetListOfBranches();
    for( const auto& branchPtr : *branch_list ){
        std::string branchName = branchPtr->GetName();
		if( stringTools::stringContains( branchName, nameToFind ) ){
			return true;
		}
	}
	return false;
}


bool TreeReader::containsGeneratorInfo() const{
    return treeHasBranchWithName( _currentTreePtr, "_gen_" );
}


bool TreeReader::containsSusyMassInfo() const{
    return treeHasBranchWithName( _currentTreePtr, "_mChi" );
}


bool TreeReader::containsTriggerInfo( const std::string& triggerPath ) const{
    return treeHasBranchWithName( _currentTreePtr, triggerPath );
}


bool TreeReader::isData() const{
    if( _currentSamplePtr ){
        return _currentSamplePtr->isData();
    } else {
        return !containsGeneratorInfo();
    }
}


bool TreeReader::isMC() const{
    return !isData();
}


void TreeReader::checkCurrentSample() const{
    if( !_currentSamplePtr ){
        throw std::domain_error( "pointer to current Sample is nullptr." );
    }
}


void TreeReader::checkCurrentTree() const{
    if( !_currentTreePtr ){
        throw std::domain_error( "pointer to current TTree is nullptr." );
    }
}


void TreeReader::checkCurrentFile() const{
	if( !_currentFilePtr ){
		throw std::domain_error( "pointer to current TFile is nullptr." );
	}
}


bool TreeReader::is2016() const{
    checkCurrentSample();
    return _currentSamplePtr->is2016();
}


bool TreeReader::is2017() const{
    checkCurrentSample();
    return _currentSamplePtr->is2017();
}


bool TreeReader::is2018() const{
    checkCurrentSample();
    return _currentSamplePtr->is2018();
}


bool TreeReader::isSMSignal() const{
    checkCurrentSample();
    return _currentSamplePtr->isSMSignal();
}


bool TreeReader::isNewPhysicsSignal() const{
    checkCurrentSample();
    return _currentSamplePtr->isNewPhysicsSignal();
}


long unsigned TreeReader::numberOfEntries() const{
    checkCurrentTree();
    return _currentTreePtr->GetEntries();
}


void TreeReader::initSample( const Sample& samp ){ 

    //update current sample
    //I wonder if the extra copy can be avoided here, its however hard if we want to keep the functionality of reading the sample vector, and also having the function initSampleFromFile. It's not clear how we can make a new sample in one of them and refer to an existing one in the other. It can be done with a static Sample in 'initSampleFromFile', but this makes the entire TreeReader class unthreadsafe, so no parallel sample processing in one process can be done 
    _currentSamplePtr = std::make_shared< Sample >( samp );
    _currentFilePtr = samp.filePtr();

    //Warning: this pointer is overwritten, but it is not a memory leak. ROOT is dirty and deletes the previous tree upon closure of the TFile it belongs to.
    //The previous TFile is closed by the std::shared_ptr destructor, implicitly called above when opening a new TFile.
    _currentTreePtr = (TTree*) _currentFilePtr->Get( "blackJackAndHookers/blackJackAndHookersTree" );
    checkCurrentTree();
    initTree();
    if( !samp.isData() ){

        //read sum of simulated event weights
        TH1D* hCounter = new TH1D( "hCounter", "Events counter", 1, 0, 1 );
        _currentFilePtr->cd( "blackJackAndHookers" );
        hCounter->Read( "hCounter" ); 
        double sumSimulatedEventWeights = hCounter->GetBinContent(1);
        delete hCounter;

        //event weights set with lumi depending on sample's era 
        double dataLumi;
        if( is2016() ){
            dataLumi = lumi::lumi2016;
        } else if( is2017() ){
            dataLumi = lumi::lumi2017;
        } else {
            dataLumi = lumi::lumi2018;
        }
        scale = samp.xSec()*dataLumi*1000 / sumSimulatedEventWeights;
    }

    //check whether current sample is a SUSY sample
    _isSusy = containsSusyMassInfo();
}


//initialize the next sample in the list
void TreeReader::initSample(){
    initSample( samples[ ++currentSampleIndex ] );
}


//initialize the current Sample directly from a root file, this is used when skimming
void TreeReader::initSampleFromFile( const std::string& pathToFile, const bool is2017, const bool is2018, const bool resetTriggersAndFilters ){

    //check if file exists 
    if( !systemTools::fileExists( pathToFile ) ){
        throw std::invalid_argument( "File '" + pathToFile + "' does not exist." );
    }

    _currentFilePtr = std::shared_ptr< TFile >( new TFile( pathToFile.c_str() ) );

    //Warning: this pointer is overwritten, but it is not a memory leak. ROOT is dirty and deletes the previous tree upon closure of the TFile it belongs to.
    //The previous TFile is closed by the std::shared_ptr destructor, implicitly called above when opening a new TFile.
    _currentTreePtr = (TTree*) _currentFilePtr->Get( "blackJackAndHookers/blackJackAndHookersTree" );
    checkCurrentTree();

    //make a new sample, and make sure the pointer remains valid
    //new is no option here since this would also require a destructor for the class which does not work for the other initSample case
    _currentSamplePtr = std::make_shared< Sample >( pathToFile, is2017, is2018, isData() );

    //initialize tree
    initTree( resetTriggersAndFilters );

    //check whether current sample is a SUSY sample
    _isSusy = containsSusyMassInfo();
}


//automatically determine whether sample is 2017 or 2018 from file name 
void TreeReader::initSampleFromFile( const std::string& pathToFile, const bool resetTriggersAndFilters ){
    
    std::pair< bool, bool > is2017Or2018 = analysisTools::fileIs2017Or2018( pathToFile );
    
    initSampleFromFile( pathToFile, is2017Or2018.first, is2017Or2018.second, resetTriggersAndFilters );
}


void TreeReader::GetEntry( const Sample& samp, long unsigned entry ){
    checkCurrentTree();

    _currentTreePtr->GetEntry( entry );

    //Set up correct event weight
    if( !samp.isData() ){
        _scaledWeight = _weight*scale;
    } else{
        _scaledWeight = 1;
    }
}


//use the currently initialized sample when running in serial
void TreeReader::GetEntry( long unsigned entry ){
    GetEntry( *_currentSamplePtr, entry );
}


Event TreeReader::buildEvent( const Sample& samp, long unsigned entry, const bool readIndividualTriggers, const bool readIndividualMetFilters ){
    GetEntry( samp, entry );
    return Event( *this, readIndividualTriggers, readIndividualMetFilters );
}


Event TreeReader::buildEvent( long unsigned entry, const bool readIndividualTriggers, const bool readIndividualMetFilters ){
    GetEntry( entry );
    return Event( *this, readIndividualTriggers, readIndividualMetFilters );
}


template< typename T > void setMapBranchAddresses( TTree* treePtr, std::map< std::string, T >& variableMap, std::map< std::string, TBranch* > branchMap ){
    for( const auto& variable : variableMap ){
        treePtr->SetBranchAddress( variable.first.c_str(), &variableMap[ variable.first ], &branchMap[ variable.first ] );
    }
}


//WARNING: this piece of code has not been generalized to other types than 'bool' because the correct ROOT identifier string needs to be used
void setMapOutputBranches( TTree* treePtr, std::map< std::string, bool >& variableMap ){
    for( const auto& variable : variableMap ){
        treePtr->Branch( variable.first.c_str(), &variableMap[ variable.first ], ( variable.first + "/O" ).c_str() );
    }    
}


void TreeReader::initTree( const bool resetTriggersAndFilters ){

    // Set branch addresses and branch pointers
    checkCurrentTree();

    _currentTreePtr->SetMakeClass(1);

    _currentTreePtr->SetBranchAddress("_runNb", &_runNb, &b__runNb);
    _currentTreePtr->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
    _currentTreePtr->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
    _currentTreePtr->SetBranchAddress("_nVertex", &_nVertex, &b__nVertex);    
    _currentTreePtr->SetBranchAddress("_passTrigger_e", &_passTrigger_e, &b__passTrigger_e);
    _currentTreePtr->SetBranchAddress("_passTrigger_ee", &_passTrigger_ee, &b__passTrigger_ee);
    _currentTreePtr->SetBranchAddress("_passTrigger_eee", &_passTrigger_eee, &b__passTrigger_eee);
    _currentTreePtr->SetBranchAddress("_passTrigger_em", &_passTrigger_em, &b__passTrigger_em);
    _currentTreePtr->SetBranchAddress("_passTrigger_m", &_passTrigger_m, &b__passTrigger_m);
    _currentTreePtr->SetBranchAddress("_passTrigger_eem", &_passTrigger_eem, &b__passTrigger_eem);
    _currentTreePtr->SetBranchAddress("_passTrigger_mm", &_passTrigger_mm, &b__passTrigger_mm);
    _currentTreePtr->SetBranchAddress("_passTrigger_emm", &_passTrigger_emm, &b__passTrigger_emm);
    _currentTreePtr->SetBranchAddress("_passTrigger_mmm", &_passTrigger_mmm, &b__passTrigger_mmm);
    _currentTreePtr->SetBranchAddress("_passTrigger_et", &_passTrigger_et, &b__passTrigger_et);
    _currentTreePtr->SetBranchAddress("_passTrigger_mt", &_passTrigger_mt, &b__passTrigger_mt);
    _currentTreePtr->SetBranchAddress("_passTrigger_FR", &_passTrigger_FR, &b__passTrigger_FR);
    _currentTreePtr->SetBranchAddress("_passTrigger_FR_iso", &_passTrigger_FR_iso, &b__passTrigger_FR_iso);
    _currentTreePtr->SetBranchAddress("_passMETFilters", &_passMETFilters, &b__passMETFilters);
    _currentTreePtr->SetBranchAddress("_nL", &_nL, &b__nL);
    _currentTreePtr->SetBranchAddress("_nMu", &_nMu, &b__nMu);
    _currentTreePtr->SetBranchAddress("_nEle", &_nEle, &b__nEle);
    _currentTreePtr->SetBranchAddress("_nLight", &_nLight, &b__nLight);
    _currentTreePtr->SetBranchAddress("_nTau", &_nTau, &b__nTau);
    _currentTreePtr->SetBranchAddress("_lPt", _lPt, &b__lPt);
    _currentTreePtr->SetBranchAddress("_lPtCorr", _lPtCorr, &b__lPtCorr);
    _currentTreePtr->SetBranchAddress("_lEta", _lEta, &b__lEta);
    _currentTreePtr->SetBranchAddress("_lEtaSC", _lEtaSC, &b__lEtaSC);
    _currentTreePtr->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
    _currentTreePtr->SetBranchAddress("_lE", _lE, &b__lE);
    _currentTreePtr->SetBranchAddress("_lECorr", _lECorr, &b__lECorr);
    _currentTreePtr->SetBranchAddress("_lFlavor", _lFlavor, &b__lFlavor);
    _currentTreePtr->SetBranchAddress("_lCharge", _lCharge, &b__lCharge);
    _currentTreePtr->SetBranchAddress("_dxy", _dxy, &b__dxy);
    _currentTreePtr->SetBranchAddress("_dz", _dz, &b__dz);
    _currentTreePtr->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
    _currentTreePtr->SetBranchAddress("_3dIPSig", _3dIPSig, &b__3dIPSig);
    _currentTreePtr->SetBranchAddress("_lElectronSummer16MvaGP", _lElectronSummer16MvaGP, &b__lElectronSummer16MvaGP);
    _currentTreePtr->SetBranchAddress("_lElectronSummer16MvaHZZ", _lElectronSummer16MvaHZZ, &b__lElectronSummer16MvaHZZ);
    _currentTreePtr->SetBranchAddress("_lElectronMvaFall17Iso", _lElectronMvaFall17Iso, &b__lElectronMvaFall17Iso);
    _currentTreePtr->SetBranchAddress("_lElectronMvaFall17NoIso", _lElectronMvaFall17NoIso, &b__lElectronMvaFall17NoIso);
    _currentTreePtr->SetBranchAddress("_lElectronPassMVAFall17NoIsoWPLoose", _lElectronPassMVAFall17NoIsoWPLoose, &b__lElectronPassMVAFall17NoIsoWPLoose);
    _currentTreePtr->SetBranchAddress("_lElectronPassMVAFall17NoIsoWP90", _lElectronPassMVAFall17NoIsoWP90, &b__lElectronPassMVAFall17NoIsoWP90);
    _currentTreePtr->SetBranchAddress("_lElectronPassMVAFall17NoIsoWP80", _lElectronPassMVAFall17NoIsoWP80, &b__lElectronPassMVAFall17NoIsoWP80);
    _currentTreePtr->SetBranchAddress("_lElectronPassEmu", _lElectronPassEmu, &b__lElectronPassEmu);
    _currentTreePtr->SetBranchAddress("_lElectronPassConvVeto", _lElectronPassConvVeto, &b__lElectronPassConvVeto);
    _currentTreePtr->SetBranchAddress("_lElectronChargeConst", _lElectronChargeConst, &b__lElectronChargeConst);
    _currentTreePtr->SetBranchAddress("_lElectronMissingHits", _lElectronMissingHits, &b__lElectronMissingHits);
    _currentTreePtr->SetBranchAddress("_lElectronEInvMinusPInv", _lElectronEInvMinusPInv, &b__lElectronEInvMinusPInv);
    _currentTreePtr->SetBranchAddress("_lElectronHOverE", _lElectronHOverE, &b__lElectronHOverE);
    _currentTreePtr->SetBranchAddress("_lElectronSigmaIetaIeta", _lElectronSigmaIetaIeta, &b__lElectronSigmaIetaIeta);
    _currentTreePtr->SetBranchAddress("_leptonMvaTTH", _leptonMvaTTH, &b__leptonMvaTTH);
    _currentTreePtr->SetBranchAddress("_leptonMvatZq", _leptonMvatZq, &b__leptonMvatZq);
    _currentTreePtr->SetBranchAddress("_leptonMvaTOP", _leptonMvaTOP, &b__leptonMvaTOP);
    _currentTreePtr->SetBranchAddress("_lPOGVeto", _lPOGVeto, &b__lPOGVeto);
    _currentTreePtr->SetBranchAddress("_lPOGLoose", _lPOGLoose, &b__lPOGLoose);
    _currentTreePtr->SetBranchAddress("_lPOGMedium", _lPOGMedium, &b__lPOGMedium);
    _currentTreePtr->SetBranchAddress("_lPOGTight", _lPOGTight, &b__lPOGTight);

    _currentTreePtr->SetBranchAddress("_tauDecayMode", _tauDecayMode, &b__tauDecayMode);
    _currentTreePtr->SetBranchAddress("_decayModeFinding", _decayModeFinding, &b__decayModeFinding);
    _currentTreePtr->SetBranchAddress("_decayModeFindingNew", _decayModeFindingNew, &b__decayModeFindingNew);
    _currentTreePtr->SetBranchAddress("_tauMuonVetoLoose", _tauMuonVetoLoose, &b__tauMuonVetoLoose);
    _currentTreePtr->SetBranchAddress("_tauMuonVetoTight", _tauMuonVetoTight, &b__tauMuonVetoTight);
    _currentTreePtr->SetBranchAddress("_tauEleVetoVLoose", _tauEleVetoVLoose, &b__tauEleVetoVLoose);
    _currentTreePtr->SetBranchAddress("_tauEleVetoLoose", _tauEleVetoLoose, &b__tauEleVetoLoose);
    _currentTreePtr->SetBranchAddress("_tauEleVetoMedium", _tauEleVetoMedium, &b__tauEleVetoMedium);
    _currentTreePtr->SetBranchAddress("_tauEleVetoTight", _tauEleVetoTight, &b__tauEleVetoTight);
    _currentTreePtr->SetBranchAddress("_tauEleVetoVTight", _tauEleVetoVTight, &b__tauEleVetoVTight);
    _currentTreePtr->SetBranchAddress("_tauPOGVLoose2015", _tauPOGVLoose2015, &b__tauPOGVLoose2015);
    _currentTreePtr->SetBranchAddress("_tauPOGLoose2015", _tauPOGLoose2015, &b__tauPOGLoose2015);
    _currentTreePtr->SetBranchAddress("_tauPOGMedium2015", _tauPOGMedium2015, &b__tauPOGMedium2015);
    _currentTreePtr->SetBranchAddress("_tauPOGTight2015", _tauPOGTight2015, &b__tauPOGTight2015);
    _currentTreePtr->SetBranchAddress("_tauPOGVTight2015", _tauPOGVTight2015, &b__tauPOGVTight2015);
    _currentTreePtr->SetBranchAddress("_tauVLooseMvaNew2015", _tauVLooseMvaNew2015, &b__tauVLooseMvaNew2015);
    _currentTreePtr->SetBranchAddress("_tauLooseMvaNew2015", _tauLooseMvaNew2015, &b__tauLooseMvaNew2015);
    _currentTreePtr->SetBranchAddress("_tauMediumMvaNew2015", _tauMediumMvaNew2015, &b__tauMediumMvaNew2015);
    _currentTreePtr->SetBranchAddress("_tauTightMvaNew2015", _tauTightMvaNew2015, &b__tauTightMvaNew2015);
    _currentTreePtr->SetBranchAddress("_tauVTightMvaNew2015", _tauVTightMvaNew2015, &b__tauVTightMvaNew2015);
    _currentTreePtr->SetBranchAddress("_tauPOGVVLoose2017v2", _tauPOGVVLoose2017v2, &b__tauPOGVVLoose2017v2);
    _currentTreePtr->SetBranchAddress("_tauPOGVTight2017v2", _tauPOGVTight2017v2, &b__tauPOGVTight2017v2);
    _currentTreePtr->SetBranchAddress("_tauPOGVVTight2017v2", _tauPOGVVTight2017v2, &b__tauPOGVVTight2017v2);
    _currentTreePtr->SetBranchAddress("_tauVLooseMvaNew2017v2", _tauVLooseMvaNew2017v2, &b__tauVLooseMvaNew2017v2);
    _currentTreePtr->SetBranchAddress("_tauLooseMvaNew2017v2", _tauLooseMvaNew2017v2, &b__tauLooseMvaNew2017v2);
    _currentTreePtr->SetBranchAddress("_tauMediumMvaNew2017v2", _tauMediumMvaNew2017v2, &b__tauMediumMvaNew2017v2);
    _currentTreePtr->SetBranchAddress("_tauTightMvaNew2017v2", _tauTightMvaNew2017v2, &b__tauTightMvaNew2017v2);
    _currentTreePtr->SetBranchAddress("_tauVTightMvaNew2017v2", _tauVTightMvaNew2017v2, &b__tauVTightMvaNew2017v2);

    _currentTreePtr->SetBranchAddress("_relIso", _relIso, &b__relIso);
    _currentTreePtr->SetBranchAddress("_relIso0p4", _relIso0p4, &b__relIso0p4);
    _currentTreePtr->SetBranchAddress("_relIso0p4MuDeltaBeta", _relIso0p4MuDeltaBeta, &b__relIso0p4MuDeltaBeta);
    _currentTreePtr->SetBranchAddress("_miniIso", _miniIso, &b__miniIso);
    _currentTreePtr->SetBranchAddress("_miniIsoCharged", _miniIsoCharged, &b__miniIsoCharged);
    _currentTreePtr->SetBranchAddress("_ptRel", _ptRel, &b__ptRel);
    _currentTreePtr->SetBranchAddress("_ptRatio", _ptRatio, &b__ptRatio);
    _currentTreePtr->SetBranchAddress("_closestJetCsvV2", _closestJetCsvV2, &b__closestJetCsvV2);
    _currentTreePtr->SetBranchAddress("_closestJetDeepCsv_b", _closestJetDeepCsv_b, &b__closestJetDeepCsv_b);
    _currentTreePtr->SetBranchAddress("_closestJetDeepCsv_bb", _closestJetDeepCsv_bb, &b__closestJetDeepCsv_bb);
    _currentTreePtr->SetBranchAddress("_closestJetDeepFlavor_b", _closestJetDeepFlavor_b, &b__closestJetDeepFlavor_b);
    _currentTreePtr->SetBranchAddress("_closestJetDeepFlavor_bb", _closestJetDeepFlavor_bb, &b__closestJetDeepFlavor_bb);
    _currentTreePtr->SetBranchAddress("_closestJetDeepFlavor_lepb", _closestJetDeepFlavor_lepb, &b__closestJetDeepFlavor_lepb);
    _currentTreePtr->SetBranchAddress("_selectedTrackMult", _selectedTrackMult, &b__selectedTrackMult);
    _currentTreePtr->SetBranchAddress("_lMuonSegComp", _lMuonSegComp, &b__lMuonSegComp);
    _currentTreePtr->SetBranchAddress("_lMuonTrackPt", _lMuonTrackPt, &b__lMuonTrackPt);
    _currentTreePtr->SetBranchAddress("_lMuonTrackPtErr", _lMuonTrackPtErr, &b__lMuonTrackPtErr);
    _currentTreePtr->SetBranchAddress("_nJets", &_nJets, &b__nJets);
    _currentTreePtr->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
    _currentTreePtr->SetBranchAddress("_jetSmearedPt", _jetSmearedPt, &b__jetSmearedPt);
    _currentTreePtr->SetBranchAddress("_jetSmearedPt_JECDown", _jetSmearedPt_JECDown, &b__jetSmearedPt_JECDown);
    _currentTreePtr->SetBranchAddress("_jetSmearedPt_JECUp", _jetSmearedPt_JECUp, &b__jetSmearedPt_JECUp);
    _currentTreePtr->SetBranchAddress("_jetSmearedPt_JERDown", _jetSmearedPt_JERDown, &b__jetSmearedPt_JERDown);
    _currentTreePtr->SetBranchAddress("_jetSmearedPt_JERUp", _jetSmearedPt_JERUp, &b__jetSmearedPt_JERUp);
    _currentTreePtr->SetBranchAddress("_jetPt_JECUp", _jetPt_JECUp, &b__jetPt_JECUp);
    _currentTreePtr->SetBranchAddress("_jetPt_JECDown", _jetPt_JECDown, &b__jetPt_JECDown);
    _currentTreePtr->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
    _currentTreePtr->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
    _currentTreePtr->SetBranchAddress("_jetE", _jetE, &b__jetE);
    _currentTreePtr->SetBranchAddress("_jetPt_Uncorrected",_jetPt_Uncorrected, &b__jetPt_Uncorrected);
    _currentTreePtr->SetBranchAddress("_jetPt_L1", _jetPt_L1, &b__jetPt_L1);
    _currentTreePtr->SetBranchAddress("_jetPt_L2", _jetPt_L2, &b__jetPt_L2);
    _currentTreePtr->SetBranchAddress("_jetPt_L3", _jetPt_L3, &b__jetPt_L3);
    _currentTreePtr->SetBranchAddress("_jetCsvV2", _jetCsvV2, &b__jetCsvV2);
    _currentTreePtr->SetBranchAddress("_jetDeepCsv_udsg", _jetDeepCsv_udsg, &b__jetDeepCsv_udsg);
    _currentTreePtr->SetBranchAddress("_jetDeepCsv_b", _jetDeepCsv_b, &b__jetDeepCsv_b);
    _currentTreePtr->SetBranchAddress("_jetDeepCsv_c", _jetDeepCsv_c, &b__jetDeepCsv_c);
    _currentTreePtr->SetBranchAddress("_jetDeepCsv_bb", _jetDeepCsv_bb, &b__jetDeepCsv_bb);
    _currentTreePtr->SetBranchAddress("_jetDeepFlavor_b", _jetDeepFlavor_b, &b__jetDeepFlavor_b);
    _currentTreePtr->SetBranchAddress("_jetDeepFlavor_bb", _jetDeepFlavor_bb, &b__jetDeepFlavor_bb);
    _currentTreePtr->SetBranchAddress("_jetDeepFlavor_lepb", _jetDeepFlavor_lepb, &b__jetDeepFlavor_lepb);
    _currentTreePtr->SetBranchAddress("_jetHadronFlavor", _jetHadronFlavor, &b__jetHadronFlavor);
    _currentTreePtr->SetBranchAddress("_jetIsTight", _jetIsTight, &b__jetIsTight);
    _currentTreePtr->SetBranchAddress("_jetIsTightLepVeto", _jetIsTightLepVeto, &b__jetIsTightLepVeto);
    _currentTreePtr->SetBranchAddress("_jetNeutralHadronFraction", _jetNeutralHadronFraction, &b__jetNeutralHadronFraction);
    _currentTreePtr->SetBranchAddress("_jetChargedHadronFraction", _jetChargedHadronFraction, &b__jetChargedHadronFraction);
    _currentTreePtr->SetBranchAddress("_jetNeutralEmFraction", _jetNeutralEmFraction, &b__jetNeutralEmFraction);
    _currentTreePtr->SetBranchAddress("_jetChargedEmFraction", _jetChargedEmFraction, &b__jetChargedEmFraction);
    _currentTreePtr->SetBranchAddress("_jetHFHadronFraction", _jetHFHadronFraction, &b__jetHFHadronFraction);
    _currentTreePtr->SetBranchAddress("_jetHFEmFraction", _jetHFEmFraction, &b__jetHFEmFraction);
    _currentTreePtr->SetBranchAddress("_met", &_met, &b__met);
    _currentTreePtr->SetBranchAddress("_met_JECDown", &_metJECDown, &b__metJECDown);
    _currentTreePtr->SetBranchAddress("_met_JECUp", &_metJECUp, &b__metJECUp);
    _currentTreePtr->SetBranchAddress("_met_UnclDown", &_metUnclDown, &b__metUnclDown);
    _currentTreePtr->SetBranchAddress("_met_UnclUp", &_metUnclUp, &b__metUnclUp);
    _currentTreePtr->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
    _currentTreePtr->SetBranchAddress("_metPhi_JECDown", &_metPhiJECDown, &b__metPhiJECDown);
    _currentTreePtr->SetBranchAddress("_metPhi_JECUp", &_metPhiJECUp, &b__metPhiJECUp);
    _currentTreePtr->SetBranchAddress("_metPhi_UnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
    _currentTreePtr->SetBranchAddress("_metPhi_UnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
    _currentTreePtr->SetBranchAddress("_metSignificance", &_metSignificance, &b__metSignificance);
    
    if( containsGeneratorInfo() ){
        _currentTreePtr->SetBranchAddress("_weight", &_weight, &b__weight);
        _currentTreePtr->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
        _currentTreePtr->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
        _currentTreePtr->SetBranchAddress("_nPsWeights", &_nPsWeights, &b__nPsWeights);
        _currentTreePtr->SetBranchAddress("_psWeight", _psWeight, &b__psWeight);
        _currentTreePtr->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
        _currentTreePtr->SetBranchAddress("_lheHTIncoming", &_lheHTIncoming, &b__lheHTIncoming);
        _currentTreePtr->SetBranchAddress("_gen_met", &_gen_met, &b__gen_met);
        _currentTreePtr->SetBranchAddress("_gen_metPhi", &_gen_metPhi, &b__gen_metPhi);
        _currentTreePtr->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
        _currentTreePtr->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
        _currentTreePtr->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
        _currentTreePtr->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
        _currentTreePtr->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
        _currentTreePtr->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
        _currentTreePtr->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
        _currentTreePtr->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
        _currentTreePtr->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
        _currentTreePtr->SetBranchAddress("_lIsPrompt", _lIsPrompt, &b__lIsPrompt);
        _currentTreePtr->SetBranchAddress("_lMatchPdgId", _lMatchPdgId, &b__lMatchPdgId);
        _currentTreePtr->SetBranchAddress("_lMatchCharge", _lMatchCharge, &b__lMatchCharge);
        _currentTreePtr->SetBranchAddress("_lMomPdgId",  _lMomPdgId, &b__lMomPdgId);
        _currentTreePtr->SetBranchAddress("_lProvenance", _lProvenance, &b__lProvenance);
        _currentTreePtr->SetBranchAddress("_lProvenanceCompressed", _lProvenanceCompressed, &b__lProvenanceCompressed);
        _currentTreePtr->SetBranchAddress("_lProvenanceConversion", _lProvenanceConversion, &b__lProvenanceConversion);
        _currentTreePtr->SetBranchAddress("_ttgEventType", &_ttgEventType, &b__ttgEventType);
        _currentTreePtr->SetBranchAddress("_zgEventType", &_zgEventType, &b__zgEventType);
        _currentTreePtr->SetBranchAddress("_gen_n"                                      , &_gen_n                                     , &b__gen_n);
        _currentTreePtr->SetBranchAddress("_gen_pt"                                     , _gen_pt                                     , &b__gen_pt);
        _currentTreePtr->SetBranchAddress("_gen_eta"                                    , _gen_eta                                    , &b__gen_eta);
        _currentTreePtr->SetBranchAddress("_gen_phi"                                    , _gen_phi                                    , &b__gen_phi);
        _currentTreePtr->SetBranchAddress("_gen_E"                                      , _gen_E                                      , &b__gen_E);
        _currentTreePtr->SetBranchAddress("_gen_pdgId"                                  , _gen_pdgId                                  , &b__gen_pdgId);
        _currentTreePtr->SetBranchAddress("_gen_charge"                                 , _gen_charge                                 , &b__gen_charge);
        _currentTreePtr->SetBranchAddress("_gen_status"                                 , _gen_status                                 , &b__gen_status);
        _currentTreePtr->SetBranchAddress("_gen_isPromptFinalState"                     , _gen_isPromptFinalState                     , &b__gen_isPromptFinalState);
        _currentTreePtr->SetBranchAddress("_gen_isDirectPromptTauDecayProductFinalState", _gen_isDirectPromptTauDecayProductFinalState, &b__gen_isDirectPromptTauDecayProductFinalState);
        _currentTreePtr->SetBranchAddress("_gen_isLastCopy"                             , _gen_isLastCopy                             , &b__gen_isLastCopy);
        _currentTreePtr->SetBranchAddress("_gen_index"                                  , _gen_index                                  , &b__gen_index);
        _currentTreePtr->SetBranchAddress("_gen_motherIndex"                            , _gen_motherIndex                            , &b__gen_motherIndex);
        _currentTreePtr->SetBranchAddress("_gen_daughter_n"                             , _gen_daughter_n                             , &b__gen_daughter_n);
        _currentTreePtr->SetBranchAddress("_gen_daughterIndex"                          , _gen_daughterIndex                          , &b__gen_daughterIndex);
    } 

    if( !is2018() && isMC() ){
        _currentTreePtr->SetBranchAddress("_prefireWeight", &_prefireWeight, &b__prefireWeight);
        _currentTreePtr->SetBranchAddress("_prefireWeightDown", &_prefireWeightDown, &b__prefireWeightDown);
        _currentTreePtr->SetBranchAddress("_prefireWeightUp", &_prefireWeightUp, &b__prefireWeightUp);
    }

	if( containsSusyMassInfo() ){
		_currentTreePtr->SetBranchAddress("_mChi1", &_mChi1, &b__mChi1);
		_currentTreePtr->SetBranchAddress("_mChi2", &_mChi2, &b__mChi2);
	}

    //add all individually stored triggers 
    //always reset triggers instead of rare case of combining primary datasets to prevent invalidating addresses set by setOutputTree
    if( resetTriggersAndFilters || _triggerMap.empty() ){
        initializeTriggerMap( _currentTreePtr );
    }
    setMapBranchAddresses( _currentTreePtr, _triggerMap, b__triggerMap );

    //add all individually stored MET filters
    //always reset filters instead of rare case of combining primary datasets to prevent invalidating addresses set by setOutputTree
    if( resetTriggersAndFilters || _MetFilterMap.empty() ){
        initializeMetFilterMap( _currentTreePtr );
    }
    setMapBranchAddresses( _currentTreePtr, _MetFilterMap, b__MetFilterMap );
}


void TreeReader::setOutputTree( TTree* outputTree ){
    outputTree->Branch("_runNb",                        &_runNb,                        "_runNb/l");
    outputTree->Branch("_lumiBlock",                    &_lumiBlock,                    "_lumiBlock/l");
    outputTree->Branch("_eventNb",                      &_eventNb,                      "_eventNb/l");
    outputTree->Branch("_nVertex",                      &_nVertex,                      "_nVertex/i");
    outputTree->Branch("_met",                          &_met,                          "_met/D");
    outputTree->Branch("_met_JECDown",                   &_metJECDown,                   "_metJECDown/D");
    outputTree->Branch("_met_JECUp",                     &_metJECUp,                     "_metJECUp/D");
    outputTree->Branch("_met_UnclDown",                  &_metUnclDown,                  "_metUnclDown/D");
    outputTree->Branch("_met_UnclUp",                    &_metUnclUp,                    "_metUnclUp/D");
    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
    outputTree->Branch("_metPhi_JECDown",                &_metPhiJECDown,                "_metPhiJECDown/D");
    outputTree->Branch("_metPhi_JECUp",                  &_metPhiJECUp,                  "_metPhiJECUp/D");
    outputTree->Branch("_metPhi_UnclDown",               &_metPhiUnclDown,               "_metPhiUnclDown/D");
    outputTree->Branch("_metPhi_UnclUp",                 &_metPhiUnclUp,                 "_metPhiUnclUp/D");
    outputTree->Branch("_metSignificance",              &_metSignificance,              "_metSignificance/D");
    outputTree->Branch("_passTrigger_e", &_passTrigger_e, "_passTrigger_e/O");
    outputTree->Branch("_passTrigger_ee", &_passTrigger_ee, "_passTrigger_ee/O");
    outputTree->Branch("_passTrigger_eee", &_passTrigger_eee, "_passTrigger_eee/O");
    outputTree->Branch("_passTrigger_em", &_passTrigger_em, "_passTrigger_em/O");
    outputTree->Branch("_passTrigger_m", &_passTrigger_m, "_passTrigger_m/O");
    outputTree->Branch("_passTrigger_eem", &_passTrigger_eem, "_passTrigger_eem/O");
    outputTree->Branch("_passTrigger_mm", &_passTrigger_mm, "_passTrigger_mm/O");
    outputTree->Branch("_passTrigger_emm", &_passTrigger_emm, "_passTrigger_emm/O");
    outputTree->Branch("_passTrigger_mmm", &_passTrigger_mmm, "_passTrigger_mmm/O");
    outputTree->Branch("_passTrigger_et", &_passTrigger_et, "_passTrigger_et/O");
    outputTree->Branch("_passTrigger_mt", &_passTrigger_mt, "_passTrigger_mt/O");
    outputTree->Branch("_passTrigger_FR", &_passTrigger_FR, "_passTrigger_FR/O");
    outputTree->Branch("_passTrigger_FR_iso", &_passTrigger_FR_iso, "_passTrigger_FR_iso/O");
    outputTree->Branch("_passMETFilters", &_passMETFilters, "_passMETFilters/O");
    outputTree->Branch("_nL",                           &_nL,                           "_nL/i");
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/i");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/i");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/i");
    outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/i");
    outputTree->Branch("_lPt",                          &_lPt,                          "_lPt[_nL]/D");
    outputTree->Branch("_lPtCorr",                      &_lPtCorr,                      "_lPtCorr[_nLight]/D");
    outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
    outputTree->Branch("_lEtaSC",                       &_lEtaSC,                       "_lEtaSC[_nLight]/D");
    outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
    outputTree->Branch("_lE",                           &_lE,                           "_lE[_nL]/D");
    outputTree->Branch("_lECorr",                       &_lECorr,                       "_lECorr[_nLight]/D");
    outputTree->Branch("_lFlavor",                      &_lFlavor,                      "_lFlavor[_nL]/i");
    outputTree->Branch("_lCharge",                      &_lCharge,                      "_lCharge[_nL]/I");
    outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
    outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
    outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
    outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
    outputTree->Branch("_lElectronSummer16MvaGP",       &_lElectronSummer16MvaGP,       "_lElectronSummer16MvaGP[_nLight]/F");
    outputTree->Branch("_lElectronSummer16MvaHZZ",      &_lElectronSummer16MvaHZZ,      "_lElectronSummer16MvaHZZ[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17Iso",        &_lElectronMvaFall17Iso,        "_lElectronMvaFall17Iso[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17NoIso",      &_lElectronMvaFall17NoIso,      "_lElectronMvaFall17NoIso[_nLight]/F");
    outputTree->Branch("_lElectronPassMVAFall17NoIsoWPLoose", &_lElectronPassMVAFall17NoIsoWPLoose, "_lElectronPassMVAFall17NoIsoWPLoose[_nLight]/O");
    outputTree->Branch("_lElectronPassMVAFall17NoIsoWP90", &_lElectronPassMVAFall17NoIsoWP90, "_lElectronPassMVAFall17NoIsoWP90[_nLight]/O");
    outputTree->Branch("_lElectronPassMVAFall17NoIsoWP80", &_lElectronPassMVAFall17NoIsoWP80, "_lElectronPassMVAFall17NoIsoWP80[_nLight]/O");
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_lElectronPassConvVeto",        &_lElectronPassConvVeto,        "_lElectronPassConvVeto[_nLight]/O");
    outputTree->Branch("_lElectronChargeConst",         &_lElectronChargeConst,         "_lElectronChargeConst[_nLight]/O");
    outputTree->Branch("_lElectronMissingHits",         &_lElectronMissingHits,         "_lElectronMissingHits[_nLight]/i");
    outputTree->Branch("_lElectronEInvMinusPInv",       &_lElectronEInvMinusPInv,       "_lElectronEInvMinusPInv[_nLight]/D");
    outputTree->Branch("_lElectronHOverE",              &_lElectronHOverE,              "_lElectronHOverE[_nLight]/D");
    outputTree->Branch("_lElectronSigmaIetaIeta",       &_lElectronSigmaIetaIeta,       "_lElectronSigmaIetaIeta[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH",                 &_leptonMvaTTH,                 "_leptonMvaTTH[_nLight]/D");
    outputTree->Branch("_leptonMvatZq",                 &_leptonMvatZq,                 "_leptonMvatZq[_nLight]/D");
    outputTree->Branch("_leptonMvaTOP",                 &_leptonMvaTOP,                 "_leptonMvaTOP[_nLight]/D");
    outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
    outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
    outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");

    outputTree->Branch("_tauDecayMode",                 &_tauDecayMode,                 "_tauDecayMode[_nL]/i");
    outputTree->Branch("_decayModeFinding",             &_decayModeFinding,             "_decayModeFinding[_nL]/O");
   	outputTree->Branch("_decayModeFindingNew",          &_decayModeFindingNew,          "_decayModeFindingNew[_nL]/O");
    outputTree->Branch("_tauPOGVLoose2015",             &_tauPOGVLoose2015,             "_tauPOGVLoose2015[_nL]/O");
    outputTree->Branch("_tauPOGLoose2015",              &_tauPOGLoose2015,              "_tauPOGLoose2015[_nL]/O");
    outputTree->Branch("_tauPOGMedium2015",             &_tauPOGMedium2015,             "_tauPOGMedium2015[_nL]/O");
    outputTree->Branch("_tauPOGTight2015",              &_tauPOGTight2015,              "_tauPOGTight2015[_nL]/O");
    outputTree->Branch("_tauPOGVTight2015",             &_tauPOGVTight2015,             "_tauPOGVTight2015[_nL]/O");
    outputTree->Branch("_tauVLooseMvaNew2015",          &_tauVLooseMvaNew2015,          "_tauVLooseMvaNew2015[_nL]/O");
    outputTree->Branch("_tauLooseMvaNew2015",           &_tauLooseMvaNew2015,           "_tauLooseMvaNew2015[_nL]/O");
    outputTree->Branch("_tauMediumMvaNew2015",          &_tauMediumMvaNew2015,          "_tauMediumMvaNew2015[_nL]/O");
    outputTree->Branch("_tauTightMvaNew2015",           &_tauTightMvaNew2015,           "_tauTightMvaNew2015[_nL]/O");
    outputTree->Branch("_tauVTightMvaNew2015",          &_tauVTightMvaNew2015,          "_tauVTightMvaNew2015[_nL]/O");
    outputTree->Branch("_tauPOGVVLoose2017v2",          &_tauPOGVVLoose2017v2,          "_tauPOGVVLoose2017v2[_nL]/O");
    outputTree->Branch("_tauPOGVTight2017v2",           &_tauPOGVTight2017v2,           "_tauPOGVTight2017v2[_nL]/O");
    outputTree->Branch("_tauPOGVVTight2017v2",          &_tauPOGVVTight2017v2,          "_tauPOGVVTight2017v2[_nL]/O");
    outputTree->Branch("_tauVLooseMvaNew2017v2",        &_tauVLooseMvaNew2017v2,        "_tauVLooseMvaNew2017v2[_nL]/O");
    outputTree->Branch("_tauLooseMvaNew2017v2",         &_tauLooseMvaNew2017v2,         "_tauLooseMvaNew2017v2[_nL]/O");
    outputTree->Branch("_tauMediumMvaNew2017v2",        &_tauMediumMvaNew2017v2,        "_tauMediumMvaNew2017v2[_nL]/O");
    outputTree->Branch("_tauTightMvaNew2017v2",         &_tauTightMvaNew2017v2,         "_tauTightMvaNew2017v2[_nL]/O");
    outputTree->Branch("_tauVTightMvaNew2017v2",        &_tauVTightMvaNew2017v2,        "_tauVTightMvaNew2017v2[_nL]/O");
	outputTree->Branch("_tauMuonVetoLoose",             &_tauMuonVetoLoose,             "_tauMuonVetoLoose[_nL]/O");
    outputTree->Branch("_tauMuonVetoTight",             &_tauMuonVetoTight,             "_tauMuonVetoTight[_nL]/O");
    outputTree->Branch("_tauEleVetoVLoose",             &_tauEleVetoVLoose,             "_tauEleVetoVLoose[_nL]/O");
    outputTree->Branch("_tauEleVetoLoose",              &_tauEleVetoLoose,              "_tauEleVetoLoose[_nL]/O");
    outputTree->Branch("_tauEleVetoMedium",             &_tauEleVetoMedium,             "_tauEleVetoMedium[_nL]/O");
    outputTree->Branch("_tauEleVetoTight",              &_tauEleVetoTight,              "_tauEleVetoTight[_nL]/O");
    outputTree->Branch("_tauEleVetoVTight",             &_tauEleVetoVTight,             "_tauEleVetoVTight[_nL]/O"); 

    outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
    outputTree->Branch("_relIso0p4",                    &_relIso0p4,                    "_relIso0p4[_nLight]/D");
    outputTree->Branch("_relIso0p4MuDeltaBeta",         &_relIso0p4MuDeltaBeta,         "_relIso0p4MuDeltaBeta[_nMu]/D");
    outputTree->Branch("_miniIso",                      &_miniIso,                      "_miniIso[_nLight]/D");
    outputTree->Branch("_miniIsoCharged",               &_miniIsoCharged,               "_miniIsoCharged[_nLight]/D");
    outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
    outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
    outputTree->Branch("_closestJetCsvV2",              &_closestJetCsvV2,              "_closestJetCsvV2[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_b",          &_closestJetDeepCsv_b,          "_closestJetDeepCsv_b[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_bb",         &_closestJetDeepCsv_bb,         "_closestJetDeepCsv_bb[_nLight]/D");
	outputTree->Branch("_closestJetDeepFlavor_b",       &_closestJetDeepFlavor_b,       "_closestJetDeepFlavor_b[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor_bb",      &_closestJetDeepFlavor_bb,      "_closestJetDeepFlavor_bb[_nLight]/D");
    outputTree->Branch("_closestJetDeepFlavor_lepb",    &_closestJetDeepFlavor_lepb,    "_closestJetDeepFlavor_lepb[_nLight]/D");
    outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
    outputTree->Branch("_lMuonSegComp",                 &_lMuonSegComp,                 "_lMuonSegComp[_nMu]/D");
    outputTree->Branch("_lMuonTrackPt",                 &_lMuonTrackPt,                 "_lMuonTrackPt[_nMu]/D");
    outputTree->Branch("_lMuonTrackPtErr",              &_lMuonTrackPtErr,              "_lMuonTrackPtErr[_nMu]/D");
    outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/i");
    outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
    outputTree->Branch("_jetPt_JECUp",               &_jetPt_JECUp,              "_jetPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetPt_JECDown",             &_jetPt_JECDown,            "_jetPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt",              &_jetSmearedPt,             "_jetSmearedPt[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JECDown",      &_jetSmearedPt_JECDown,     "_jetSmearedPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JECUp",        &_jetSmearedPt_JECUp,       "_jetSmearedPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JERDown",      &_jetSmearedPt_JERDown,     "_jetSmearedPt_JERDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JERUp",        &_jetSmearedPt_JERUp,       "_jetSmearedPt_JERUp[_nJets]/D");
    outputTree->Branch("_jetEta",                    &_jetEta,                   "_jetEta[_nJets]/D");
    outputTree->Branch("_jetPhi",                    &_jetPhi,                   "_jetPhi[_nJets]/D");
    outputTree->Branch("_jetE",                      &_jetE,                     "_jetE[_nJets]/D");
    outputTree->Branch("_jetPt_Uncorrected",         &_jetPt_Uncorrected,        "_jetPt_Uncorrected[_nJets]/D");
    outputTree->Branch("_jetPt_L1",                  &_jetPt_L1,                 "_jetPt_L1[_nJets]/D");
    outputTree->Branch("_jetPt_L2",                  &_jetPt_L2,                 "_jetPt_L2[_nJets]/D");
    outputTree->Branch("_jetPt_L3",                  &_jetPt_L3,                 "_jetPt_L3[_nJets]/D");
    outputTree->Branch("_jetCsvV2",                  &_jetCsvV2,                 "_jetCsvV2[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_udsg",           &_jetDeepCsv_udsg,          "_jetDeepCsv_udsg[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_b",              &_jetDeepCsv_b,             "_jetDeepCsv_b[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_c",              &_jetDeepCsv_c,             "_jetDeepCsv_c[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_bb",             &_jetDeepCsv_bb,            "_jetDeepCsv_bb[_nJets]/D");
	outputTree->Branch("_jetDeepFlavor_b",           &_jetDeepFlavor_b,          "_jetDeepFlavor_b[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_bb",          &_jetDeepFlavor_bb,         "_jetDeepFlavor_bb[_nJets]/D");
    outputTree->Branch("_jetDeepFlavor_lepb",        &_jetDeepFlavor_lepb,       "_jetDeepFlavor_lepb[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/O");
    outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/O");
    outputTree->Branch("_jetNeutralHadronFraction",  &_jetNeutralHadronFraction, "_jetNeutralHadronFraction[_nJets]/D");
    outputTree->Branch("_jetChargedHadronFraction",  &_jetChargedHadronFraction, "_jetChargedHadronFraction[_nJets]/D");
    outputTree->Branch("_jetNeutralEmFraction",      &_jetNeutralEmFraction,     "_jetNeutralEmFraction[_nJets]/D");
    outputTree->Branch("_jetChargedEmFraction",      &_jetChargedEmFraction,     "_jetChargedEmFraction[_nJets]/D");
    outputTree->Branch("_jetHFHadronFraction",       &_jetHFHadronFraction,      "_jetHFHadronFraction[_nJets]/D");
    outputTree->Branch("_jetHFEmFraction",           &_jetHFEmFraction,          "_jetHFEmFraction[_nJets]/D");


    if( containsGeneratorInfo() ){
        outputTree->Branch("_nLheWeights",               &_nLheWeights,               "_nLheWeights/i");
        outputTree->Branch("_lheWeight",                 &_lheWeight,                 "_lheWeight[_nLheWeights]/D");
        outputTree->Branch("_weight",                    &_weight,                    "_weight/D");
        outputTree->Branch("_nPsWeights",                &_nPsWeights,                "_nPsWeights/i");
        outputTree->Branch("_psWeight",                  &_psWeight,                  "_psWeight/D");
        outputTree->Branch("_nTrueInt",                  &_nTrueInt,                  "_nTrueInt/F");
        outputTree->Branch("_lheHTIncoming",             &_lheHTIncoming,             "_lheHTIncoming/D");
        outputTree->Branch("_lIsPrompt",                 &_lIsPrompt,                 "_lIsPrompt[_nL]/O");
        outputTree->Branch("_lMatchPdgId",               &_lMatchPdgId,               "_lMatchPdgId[_nL]/I");
        outputTree->Branch("_lMatchCharge",              &_lMatchCharge,              "_lMatchCharge[_nL]/I");
        outputTree->Branch("_lMomPdgId",                 &_lMomPdgId,                 "_lMomPdgId[_nL]/I");
        outputTree->Branch("_lProvenance",               &_lProvenance,               "_lProvenance[_nL]/i");
        outputTree->Branch("_lProvenanceCompressed",     &_lProvenanceCompressed,     "_lProvenanceCompressed[_nL]/i");
        outputTree->Branch("_lProvenanceConversion",     &_lProvenanceConversion,     "_lProvenanceConversion[_nL]/i");
        outputTree->Branch("_gen_met",                   &_gen_met,                   "_gen_met/D");
        outputTree->Branch("_gen_metPhi",                &_gen_metPhi,                "_gen_metPhi/D");
        outputTree->Branch("_gen_nL",                    &_gen_nL,                    "_gen_nL/i");
        outputTree->Branch("_gen_lPt",                   &_gen_lPt,                   "_gen_lPt[_gen_nL]/D");
        outputTree->Branch("_gen_lEta",                  &_gen_lEta,                  "_gen_lEta[_gen_nL]/D");
        outputTree->Branch("_gen_lPhi",                  &_gen_lPhi,                  "_gen_lPhi[_gen_nL]/D");
        outputTree->Branch("_gen_lE",                    &_gen_lE,                    "_gen_lE[_gen_nL]/D");
        outputTree->Branch("_gen_lFlavor",               &_gen_lFlavor,               "_gen_lFlavor[_gen_nL]/i");
        outputTree->Branch("_gen_lCharge",               &_gen_lCharge,               "_gen_lCharge[_gen_nL]/I");
        outputTree->Branch("_gen_lMomPdg",               &_gen_lMomPdg,               "_gen_lMomPdg[_gen_nL]/I");
        outputTree->Branch("_gen_lIsPrompt",             &_gen_lIsPrompt,             "_gen_lIsPrompt[_gen_nL]/O");
        outputTree->Branch("_ttgEventType",              &_ttgEventType,              "_ttgEventType/i");
        outputTree->Branch("_zgEventType",               &_zgEventType,               "_zgEventType/i");

        outputTree->Branch("_gen_n"                                      , &_gen_n                                      ,"_gen_n/i"); // UInt_t          _gen_n;
        outputTree->Branch("_gen_pt"                                     , &_gen_pt                                     ,"_gen_pt[_gen_n]/D"); // Double_t        _gen_pt[gen_nL_max];
        outputTree->Branch("_gen_eta"                                    , &_gen_eta                                    ,"_gen_eta[_gen_n]/D"); // Double_t        _gen_eta[gen_nL_max];
        outputTree->Branch("_gen_phi"                                    , &_gen_phi                                    ,"_gen_phi[_gen_n]/D"); // Double_t        _gen_phi[gen_nL_max];
        outputTree->Branch("_gen_E"                                      , &_gen_E                                      ,"_gen_E[_gen_n]/D"); // Double_t        _gen_E[gen_nL_max];
        outputTree->Branch("_gen_pdgId"                                  , &_gen_pdgId                                  ,"_gen_pdgId[_gen_n]/I"); // Int_t           _gen_pdgId[gen_nL_max];
        outputTree->Branch("_gen_charge"                                 , &_gen_charge                                 ,"_gen_charge[_gen_n]/I"); // Int_t           _gen_charge[gen_nL_max];
        outputTree->Branch("_gen_status"                                 , &_gen_status                                 ,"_gen_status[_gen_n]/I"); // Int_t           _gen_status[gen_nL_max];
        outputTree->Branch("_gen_isPromptFinalState"                     , &_gen_isPromptFinalState                     ,"_gen_isPromptFinalState[_gen_n]/O"); // Bool_t          _gen_isPromptFinalState[
        outputTree->Branch("_gen_isDirectPromptTauDecayProductFinalState", &_gen_isDirectPromptTauDecayProductFinalState,"_gen_isDirectPromptTauDecayProductFinalState[_gen_n]/O"); // Bool_t          _gen_isDirectPromptTauDe
        outputTree->Branch("_gen_isLastCopy"                             , &_gen_isLastCopy                             ,"_gen_isLastCopy[_gen_n]/O"); // Bool_t          _gen_isLastCopy[gen_nL_m
        outputTree->Branch("_gen_index"                                  , &_gen_index                                  ,"_gen_index[_gen_n]/I"); // Int_t           _gen_index[gen_nL_max];
        outputTree->Branch("_gen_motherIndex"                            , &_gen_motherIndex                            ,"_gen_motherIndex[_gen_n]/I"); // Int_t           _gen_motherIndex[gen_nL_
        outputTree->Branch("_gen_daughter_n"                             , &_gen_daughter_n                             ,"_gen_daughter_n[_gen_n]/I"); // Int_t           _gen_daughter_n[gen_nL_m
        outputTree->Branch("_gen_daughterIndex"                          , &_gen_daughterIndex                          ,"_gen_daughterIndex[_gen_n][10]/I"); // Int_t           _gen_daughterIndex[5000]
    } 

    if( !is2018() && isMC() ){
       	outputTree->Branch("_prefireWeight",             &_prefireWeight,             "_prefireWeight/F");
        outputTree->Branch("_prefireWeightUp",           &_prefireWeightUp,           "_prefireWeightUp/F");
        outputTree->Branch("_prefireWeightDown",         &_prefireWeightDown,         "_prefireWeightDown/F"); 
    }

    if( containsSusyMassInfo() ){
		outputTree->Branch("_mChi1", &_mChi1, "_mChi1/D");
		outputTree->Branch("_mChi2", &_mChi2, "_mChi2/D");
    }

    //write individual trigger decisions to output tree 
    setMapOutputBranches( outputTree, _triggerMap );

    //write individual MET filters to output tree
    setMapOutputBranches( outputTree, _MetFilterMap );
}


//get object from current file 
TObject* TreeReader::getFromCurrentFile( const std::string& name ) const{
	checkCurrentFile();
	return _currentFilePtr->Get( name.c_str() );
}


//Get list of histograms stored in current file
std::vector< std::shared_ptr< TH1 > > TreeReader::getHistogramsFromCurrentFile() const{

	checkCurrentFile();

	//vector containing all histograms in current file
	std::vector< std::shared_ptr< TH1 > > histogramVector;

	//loop over keys in blackJackAndHookers directory
    //this directory gets implicitly deleted by root when the current root file gets deleted. This can NOT be a shared_ptr since this directory will also make the root file inaccessible upon deletion (DIRTY ROOT!!!)
	TDirectory* dir = (TDirectory*) _currentFilePtr->Get("blackJackAndHookers");

    //this is not a memory leak since this object will implicitly be deleted when 'dir' gets deleted (DIRTY ROOT!!!)
    TList* keyList = dir->GetListOfKeys();

    for( const auto objectPtr : *keyList ){

		//try if a dynamic_cast to a histogram works to check if object is histogram
		TH1* histPtr = dynamic_cast< TH1* >( dir->Get( objectPtr->GetName() ) );
		if( histPtr ){

            //make sure histograms don't get deleted by root upon deletion of TDirectory above
            histPtr->SetDirectory( gROOT );
			histogramVector.emplace_back( histPtr );
		}
    }
	return histogramVector;
}


void TreeReader::removeBSMSignalSamples(){
    for( auto it = samples.begin(); it != samples.end(); ){
        if( it->isNewPhysicsSignal() ){
            it = samples.erase( it );
        } else {
            ++it;
        }
    }
}
