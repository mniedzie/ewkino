#ifndef TreeReader_H
#define TreeReader_H

//include ROOT classes
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"

//include other parts of code
#include "../../Tools/interface/Sample.h"


class Event;


class TreeReader {

    public :

        //Constructor
        TreeReader() = default;
        TreeReader( const std::string&, const std::string& ); 

        //Declare leaf types
        static const unsigned nL_max = 20;
        static const unsigned nJets_max = 100;
        static const unsigned gen_nL_max = 20;
        static const unsigned gen_n_max = 1000;
        ULong_t         _runNb;
        ULong_t         _lumiBlock;
        ULong_t         _eventNb;
        UInt_t          _nVertex;
        Float_t         _prefireWeight;
        Float_t         _prefireWeightDown;
        Float_t         _prefireWeightUp;
        Double_t        _weight;
        UInt_t          _nLheWeights;
        Double_t        _lheWeight[148];
        UInt_t          _nPsWeights;
        Double_t        _psWeight[14];
        Float_t         _nTrueInt;
        Double_t        _lheHTIncoming;
        Double_t        _gen_met;
        Double_t        _gen_metPhi;
        UInt_t          _gen_nL;
        Double_t        _gen_lPt[gen_nL_max];   
        Double_t        _gen_lEta[gen_nL_max];   
        Double_t        _gen_lPhi[gen_nL_max];   
        Double_t        _gen_lE[gen_nL_max];   
        UInt_t          _gen_lFlavor[gen_nL_max];   
        Int_t           _gen_lCharge[gen_nL_max];   
        Int_t           _gen_lMomPdg[gen_nL_max];   
        Bool_t          _gen_lIsPrompt[gen_nL_max];   

        UInt_t          _gen_n;
        Double_t        _gen_pt[gen_n_max];
        Double_t        _gen_eta[gen_n_max];
        Double_t        _gen_phi[gen_n_max];
        Double_t        _gen_E[gen_n_max];
        Int_t           _gen_pdgId[gen_n_max];
        Int_t           _gen_charge[gen_n_max];
        Int_t           _gen_status[gen_n_max];
        Bool_t          _gen_isPromptFinalState[gen_n_max];
        Bool_t          _gen_isDirectPromptTauDecayProductFinalState[gen_n_max];
        Bool_t          _gen_isLastCopy[gen_n_max];
        Int_t           _gen_index[gen_n_max];
        Int_t           _gen_motherIndex[gen_n_max];
        Int_t           _gen_daughter_n[gen_n_max];
        Int_t           _gen_daughterIndex[gen_n_max][10];

        UInt_t          _ttgEventType;
        UInt_t          _zgEventType;
        Bool_t          _passTrigger_e;
        Bool_t          _passTrigger_ee;
        Bool_t          _passTrigger_eee;
        Bool_t          _passTrigger_em;
        Bool_t          _passTrigger_m;
        Bool_t          _passTrigger_eem;
        Bool_t          _passTrigger_mm;
        Bool_t          _passTrigger_emm;
        Bool_t          _passTrigger_mmm;
        Bool_t          _passTrigger_et;
        Bool_t          _passTrigger_mt;
        Bool_t          _passTrigger_FR;
        Bool_t          _passTrigger_FR_iso;
        Bool_t          _passMETFilters;
        UInt_t          _nL;
        UInt_t          _nMu;
        UInt_t          _nEle;
        UInt_t          _nLight;
        UInt_t          _nTau;
        Double_t        _lPt[nL_max];   
        Double_t        _lPtCorr[nL_max];
        Double_t        _lEta[nL_max];   
        Double_t        _lEtaSC[nL_max];   
        Double_t        _lPhi[nL_max];   
        Double_t        _lE[nL_max];   
        Double_t        _lECorr[nL_max];
        UInt_t          _lFlavor[nL_max];   
        Int_t           _lCharge[nL_max];   
        Double_t        _dxy[nL_max];   
        Double_t        _dz[nL_max];   
        Double_t        _3dIP[nL_max];   
        Double_t        _3dIPSig[nL_max];   
        Float_t         _lElectronSummer16MvaGP[nL_max];   
        Float_t         _lElectronSummer16MvaHZZ[nL_max];
        Float_t         _lElectronMvaFall17Iso[nL_max];
        Float_t         _lElectronMvaFall17NoIso[nL_max];
        Bool_t          _lElectronPassMVAFall17NoIsoWPLoose[nL_max];
        Bool_t          _lElectronPassMVAFall17NoIsoWP90[nL_max];
        Bool_t          _lElectronPassMVAFall17NoIsoWP80[nL_max];
        Bool_t          _lElectronPassEmu[nL_max];   
        Bool_t          _lElectronPassConvVeto[nL_max];
        Bool_t          _lElectronChargeConst[nL_max];
        UInt_t          _lElectronMissingHits[nL_max];
        Double_t        _lElectronEInvMinusPInv[nL_max];
        Double_t        _lElectronHOverE[nL_max];
        Double_t        _lElectronSigmaIetaIeta[nL_max];
        Double_t        _leptonMvaTTH[nL_max];
        Double_t        _leptonMvatZq[nL_max];
        Double_t        _leptonMvaTOP[nL_max];
        Bool_t          _lPOGVeto[nL_max];   
        Bool_t          _lPOGLoose[nL_max];   
        Bool_t          _lPOGMedium[nL_max];   
        Bool_t          _lPOGTight[nL_max];   

        UInt_t          _tauDecayMode[nL_max];
        Bool_t          _decayModeFinding[nL_max];
        Bool_t          _decayModeFindingNew[nL_max];   
        Bool_t          _tauMuonVetoLoose[nL_max];
        Bool_t          _tauMuonVetoTight[nL_max];
        Bool_t          _tauEleVetoVLoose[nL_max];
        Bool_t          _tauEleVetoLoose[nL_max];
        Bool_t          _tauEleVetoMedium[nL_max];
        Bool_t          _tauEleVetoTight[nL_max];
        Bool_t          _tauEleVetoVTight[nL_max];
		Bool_t 			_tauPOGVLoose2015[nL_max];
		Bool_t			_tauPOGLoose2015[nL_max];
		Bool_t			_tauPOGMedium2015[nL_max];
		Bool_t			_tauPOGTight2015[nL_max];
		Bool_t			_tauPOGVTight2015[nL_max];
		Bool_t			_tauVLooseMvaNew2015[nL_max];
		Bool_t			_tauLooseMvaNew2015[nL_max];
		Bool_t			_tauMediumMvaNew2015[nL_max];
		Bool_t			_tauTightMvaNew2015[nL_max];
		Bool_t			_tauVTightMvaNew2015[nL_max];
		Bool_t			_tauPOGVVLoose2017v2[nL_max];
		Bool_t			_tauPOGVTight2017v2[nL_max];
		Bool_t			_tauPOGVVTight2017v2[nL_max];
		Bool_t			_tauVLooseMvaNew2017v2[nL_max];
		Bool_t			_tauLooseMvaNew2017v2[nL_max];
		Bool_t			_tauMediumMvaNew2017v2[nL_max];
		Bool_t			_tauTightMvaNew2017v2[nL_max];
		Bool_t			_tauVTightMvaNew2017v2[nL_max];

        Double_t        _relIso[nL_max];   
        Double_t        _relIso0p4[nL_max];
        Double_t        _relIso0p4MuDeltaBeta[nL_max];
        Double_t        _miniIso[nL_max];   
        Double_t        _miniIsoCharged[nL_max];   
        Double_t        _ptRel[nL_max];   
        Double_t        _ptRatio[nL_max];   
        Double_t        _closestJetCsvV2[nL_max];
        Double_t        _closestJetDeepCsv_b[nL_max];
        Double_t        _closestJetDeepCsv_bb[nL_max];
        Double_t        _closestJetDeepFlavor_b[nL_max];
        Double_t        _closestJetDeepFlavor_bb[nL_max];
        Double_t        _closestJetDeepFlavor_lepb[nL_max];
        UInt_t          _selectedTrackMult[nL_max];   
        Double_t        _lMuonSegComp[nL_max];   
        Double_t        _lMuonTrackPt[nL_max];
        Double_t        _lMuonTrackPtErr[nL_max];
        Bool_t          _lIsPrompt[nL_max];   
        Int_t           _lMatchPdgId[nL_max];   
        Int_t           _lMatchCharge[nL_max];
        Int_t           _lMomPdgId[nL_max];
        UInt_t          _lProvenance[nL_max];
        UInt_t          _lProvenanceCompressed[nL_max];
        UInt_t          _lProvenanceConversion[nL_max];
        UInt_t          _nJets;
        Double_t        _jetPt[nJets_max];   
        Double_t        _jetPt_JECUp[nJets_max];   
        Double_t        _jetPt_JECDown[nJets_max];   
        Double_t        _jetSmearedPt[nJets_max];
        Double_t        _jetSmearedPt_JECDown[nJets_max];
        Double_t        _jetSmearedPt_JECUp[nJets_max];
        Double_t        _jetSmearedPt_JERDown[nJets_max];
        Double_t        _jetSmearedPt_JERUp[nJets_max];
        Double_t        _jetPt_Uncorrected[nJets_max];
        Double_t        _jetPt_L1[nJets_max];
        Double_t        _jetPt_L2[nJets_max];
        Double_t        _jetPt_L3[nJets_max];
        Double_t        _jetEta[nJets_max];   
        Double_t        _jetPhi[nJets_max];   
        Double_t        _jetE[nJets_max];   
        Double_t        _jetCsvV2[nJets_max];   
        Double_t        _jetDeepCsv_udsg[nJets_max];   
        Double_t        _jetDeepCsv_b[nJets_max];   
        Double_t        _jetDeepCsv_c[nJets_max];   
        Double_t        _jetDeepCsv_bb[nJets_max];
        Double_t        _jetDeepFlavor_b[nJets_max];
        Double_t        _jetDeepFlavor_bb[nJets_max];
        Double_t        _jetDeepFlavor_lepb[nJets_max];
        UInt_t          _jetHadronFlavor[nJets_max];   
        Bool_t          _jetIsTight[nJets_max];
        Bool_t          _jetIsTightLepVeto[nJets_max];
        Double_t        _jetNeutralHadronFraction[nJets_max];
        Double_t        _jetChargedHadronFraction[nJets_max];
        Double_t        _jetNeutralEmFraction[nJets_max];
        Double_t        _jetChargedEmFraction[nJets_max];
        Double_t        _jetHFHadronFraction[nJets_max];
        Double_t        _jetHFEmFraction[nJets_max];
        Double_t        _met;
        Double_t        _metJECDown;
        Double_t        _metJECUp;
        Double_t        _metUnclDown;
        Double_t        _metUnclUp;
        Double_t        _metPhi;
        Double_t        _metPhiJECDown;
        Double_t        _metPhiJECUp;
        Double_t        _metPhiUnclDown;
        Double_t        _metPhiUnclUp;       
        Double_t        _metSignificance;
        Double_t        _mChi1;
        Double_t        _mChi2;

        std::map< std::string, bool > _triggerMap;
        std::map< std::string, bool > _MetFilterMap;

        //weight including cross section scaling 
        double          _scaledWeight;

        //set up tree for reading and writing
        //always reset triggers instead of rare case of combining primary datasets!
        void initTree( const bool resetTriggersAndFilters = true);
        void setOutputTree( TTree* );

        //void combinePD(std::vector<std::string>& datasets, const bool is2017, std::string outputDirectory = "");

        //initialize the next sample
        void initSample();
        void initSample(const Sample&);  

        //read sample list from text file
        void readSamples2016(const std::string&, const std::string&);
        void readSamples2017(const std::string&, const std::string&);
        void readSamples(const std::string& list, const std::string& directory);

        //initialize the current sample directly from a root file
        //always reset triggers instead of rare case of combining primary datasets to prevent invalidating addresses set by setOutputTree
        void initSampleFromFile( const std::string& pathToFile, const bool is2017, const bool is2018, const bool resetTriggersAndFilters = true );
        void initSampleFromFile( const std::string& pathToFile, const bool resetTriggersAndFilters = true );

        //Get entry from Tree, should not be used except for test purposes
        void GetEntry(const Sample&, long unsigned );
        void GetEntry(long unsigned );

        //Build event (this will implicitly use GetEntry )
        //Use these functions in analysis code 
        Event buildEvent( const Sample&, long unsigned , const bool readIndividualTriggers = false, const bool readIndividualMetFilters = false );
        Event buildEvent( long unsigned, const bool readIndividualTriggers = false, const bool readIndividualMetFilters = false );

        //check whether generator info is present in current tree
        bool containsGeneratorInfo() const;

        //check whether SUSY mass info is present in the current sample ( this is the case for SUSY signal scans )
        bool containsSusyMassInfo() const;

        //check whether a particular trigger is present 
        bool containsTriggerInfo( const std::string& triggerPath ) const;

        //check which year the current sample belongs to
        bool is2016() const;
        bool is2017() const;
        bool is2018() const;

        //check whether the current sample is data or MC, and or signal
        bool isData() const;
        bool isMC() const;
        bool isSMSignal() const;
        bool isNewPhysicsSignal() const;
        bool isSusy() const{ return _isSusy; }

        //access number of samples and current sample
        const Sample& currentSample() const{ return *_currentSamplePtr; }
        const Sample* currentSamplePtr() const{ return _currentSamplePtr.get(); }
        std::vector< Sample >::size_type numberOfSamples() const{ return samples.size(); }
        std::vector< Sample > sampleVector() const{ return samples; }
        void removeBSMSignalSamples();

        //access current file and tree 
        std::shared_ptr< TFile > currentFilePtr(){ return _currentFilePtr; }

        //get object from current file 
        TObject* getFromCurrentFile( const std::string& name ) const;//{ return currentFilePtr->Get( name.c_str() ); }

        //Get list of histograms stored in current file
        std::vector< std::shared_ptr< TH1 > > getHistogramsFromCurrentFile() const;

        //functions for event selection
        /*
        void orderByPt(std::vector<unsigned>&, const double*, const unsigned) const;
        unsigned dilFlavorComb(const std::vector<unsigned>&) const;
        double coneCorr(const unsigned) const;
        void applyConeCorrection();
        bool lepIsLoose(const unsigned) const;
        bool lepIsGood(const unsigned) const;
        bool lepIsTight(const unsigned) const;
        bool lepFromMEExtConversion(const unsigned) const;
        bool eleIsClean(const unsigned) const;
        double closestJetDeepCSV(const unsigned) const;

        unsigned selectLep(std::vector<unsigned>&) const;
        unsigned selectLepConeCorr(std::vector<unsigned>&);
        unsigned tightLepCount(const std::vector<unsigned>&, const unsigned) const;

        bool passPtCuts(const std::vector<unsigned>&) const;
        bool jetIsClean(const unsigned) const;
        bool jetIsGood(const unsigned, const double ptCut = 25., const unsigned unc = 0, const bool clean = true, const bool allowForward = false) const;
        unsigned nJets(const unsigned unc = 0, const bool clean = true, const bool allowForward = false) const;                                   //without jet pt ordering
        unsigned nJets(std::vector<unsigned>& jetInd, const unsigned unc = 0, const bool clean = true, const bool allowForward = false) const;    //with jet pt ordering
        double deepCSV(const unsigned) const;
        bool bTagged(const unsigned ind, const unsigned wp = 1, const bool deepCSV = true) const;
        unsigned nBJets(const unsigned unc = 0, const bool deepCSV = true, const bool clean = true, const unsigned wp = 1) const;
        unsigned nBJets(std::vector<unsigned>& bJetInd, const unsigned unc = 0, const bool deepCSV = true, const bool clean = true, const unsigned wp = 1) const;

        //baseline selection for leptonMva training
        bool lepPassBaseline(const unsigned) const;

        //trigger decitions
        bool passSingleLeptonTriggers() const;
        bool passDileptonTriggers() const;
        bool passTrileptonTriggers() const;
        bool passTriggerCocktail() const;
        bool passMETFilters() const;

        //overlap removal between samples
        bool photonOverlap(const bool mcNonprompt = true) const;                                            //sample overlap due to photons
        bool photonOverlap(const Sample&, const bool mcNonprompt = true) const;
        bool htOverlap() const;                                                                             //sample overlap due to HT binning
        bool htOverlap(const Sample&) const;

        //check if leptons are prompt in MC
        bool promptLeptons() const;

        //compute b-tagging efficiency
        void computeBTagEff(const std::string& analysis, const bool clean, const bool deepCSV, const bool is2016);

        //event weights
        //pileup reweighting
        double puWeight(const unsigned unc = 0) const;

        //b-tag reweighting
        double bTagWeight_cut_singleJet(const unsigned jetIndex, const unsigned unc = 0) const;
        double bTagWeight_reshaping_singleJet(const unsigned jetIndex, const unsigned unc = 0) const;
        double bTagWeight_base(const unsigned jetFlavor, const unsigned unc, double (TreeReader::*jetWeight)(const unsigned, const unsigned) const ) const;
        double bTagWeight_cut( const unsigned jetFlavor, const unsigned unc = 0) const;
        double bTagWeight_reshaping( const unsigned jetFlavor, const unsigned unc = 0) const;
        double bTagWeight(const unsigned jetFlavor, const unsigned unc = 0) const;
        double bTagWeight(const std::vector<unsigned>& jetInd, const unsigned jetFlavor, const unsigned unc = 0) const; //more efficient version if jets were already selected 
        double bTagWeight_udsg(const unsigned unc = 0) const;
        double bTagWeight_c(const unsigned unc = 0) const;
        double bTagWeight_b(const unsigned unc = 0) const;
        double bTagWeight(const unsigned unc = 0) const;

        //lepton reweighting
        double leptonWeight(const std::string& unc = "") const;

        //fake-rate
        double fakeRateWeight(const unsigned unc = 0);
        double sfWeight();
        double jetPrefiringWeight(const unsigned unc = 0) const;
        */

        unsigned long numberOfEntries() const;

    private:

        //list of samples to loop over 
        std::vector< Sample > samples;
        std::vector< Sample > samples2016;
        std::vector< Sample > samples2017;
        std::vector< Sample > samples2018;

        //current sample
        std::shared_ptr< const Sample > _currentSamplePtr;

        //TFile associated to current sample
        std::shared_ptr< TFile > _currentFilePtr;

        //TTree associated to current sample 
        TTree* _currentTreePtr = nullptr;

        //cache whether current sample is SUSY to avoid having to check the branch names for each event
        bool _isSusy = false;

        //check whether current sample is initialized, throw an error if it is not 
        void checkCurrentSample() const;

        //check whether current Tree is initialized, throw an error if it is not 
        void checkCurrentTree() const;

        //check whether current File is initialized, throw an error if it is not
        void checkCurrentFile() const;

        //current index in samples vector
        int currentSampleIndex = -1;

        //luminosity scaling
        double scale = 0;

        //std::shared_ptr<Reweighter> reweighter;                                 //instance of reweighter class used for reweighting functions
        //initialize SF weights
        //void initializeWeights();


        /*
        //check lepton flavors 
        bool isElectron(const unsigned leptonIndex) const { return (_lFlavor[leptonIndex] == 0); }
        bool isMuon(const unsigned leptonIndex) const { return (_lFlavor[leptonIndex] == 1); }
        bool isTau(const unsigned leptonIndex) const { return (_lFlavor[leptonIndex] == 2); }

        //era-specific event selection functions
        bool lepIsLooseBase(const unsigned) const;
        bool lepIsLoose2016(const unsigned) const;
        bool lepIsLoose2017(const unsigned) const;

        bool lepIsGoodBase(const unsigned) const;
        bool lepIsGood2016(const unsigned) const;
        bool lepIsGood2017(const unsigned) const;

        bool lepIsTightBase(const unsigned) const;
        bool lepIsTight2016(const unsigned) const;
        bool lepIsTight2017(const unsigned) const;

        bool eleIsCleanBase(const unsigned, bool (TreeReader::*looseMuon)(const unsigned) const) const;
        bool eleIsClean2016(const unsigned) const;
        bool eleIsClean2017(const unsigned) const;

        bool jetIsCleanBase(const unsigned, bool (TreeReader::*leptonIsFO)(const unsigned) const) const;

        bool bTaggedDeepCSVBase(const unsigned, const unsigned wp, const double cuts[3]) const;
        bool bTaggedDeepCSV2016(const unsigned, const unsigned wp = 1) const;
        bool bTaggedDeepCSV2017(const unsigned, const unsigned wp = 1) const;
        bool bTaggedDeepCSV(const unsigned, const unsigned wp = 1) const;

        bool bTaggedCSVv2Base(const unsigned, const unsigned wp, const double cuts[3]) const;
        bool bTaggedCSVv22016(const unsigned, const unsigned wp = 1) const;
        bool bTaggedCSVv22017(const unsigned, const unsigned wp = 1) const;
        bool bTaggedCSVv2(const unsigned, const unsigned wp = 1) const;


        //lepton selection for different parts of tZq and ttV analysis (to compute bTag efficiencies for everyone)
        bool lepIsGoodtZq(const unsigned) const;

        bool lepIsGoodttZ3l2016(const unsigned) const;
        bool lepIsGoodttZ3l2017(const unsigned) const;
        bool lepIsGoodttZ3l(const unsigned) const;

        bool lepIsGoodttZ4l2016(const unsigned) const;
        bool lepIsGoodttZ4l2017(const unsigned) const;
        bool lepIsGoodttZ4l(const unsigned) const;

        bool lepIsGoodttW2016(const unsigned) const;
        bool lepIsGoodttW2017(const unsigned) const;
        bool lepIsGoodttW(const unsigned) const;

        bool lepIsGoodMultiAnalysis(const std::string&, const unsigned) const;

        */


        //some safety-checks for errors 
        void checkSampleEraConsistency() const;  //make sure a sample is not is2016() AND 2017() 
        void checkEraOrthogonality() const;        //make sure no sample from the wrong era is being used (i.e. no 2016 sample in the list of 2017 samples) 

        //debugging prints
        //void printLeptonContent( std::ostream& os = std::cout ) const;
        //void printLeptonPairing( std::ostream& os = std::cout ) const;

        //general function to read a list of samples
        void readSamples(const std::string&, const std::string&, std::vector<Sample>&);

        //initialize triggerMap
        void initializeTriggerMap( TTree* );
        void initializeMetFilterMap( TTree* );

        //list of branches
        TBranch        *b__runNb;   
        TBranch        *b__lumiBlock;   
        TBranch        *b__eventNb;   
        TBranch        *b__nVertex;   
        TBranch        *b__prefireWeight;
        TBranch        *b__prefireWeightDown;
        TBranch        *b__prefireWeightUp;
        TBranch        *b__weight;   
        TBranch        *b__nLheWeights;   
        TBranch        *b__nPsWeights;
        TBranch        *b__psWeight;
        TBranch        *b__lheWeight;   
        TBranch        *b__nTrueInt;   
        TBranch        *b__lheHTIncoming;
        TBranch        *b__gen_met;   
        TBranch        *b__gen_metPhi;   
        TBranch        *b__gen_nL;   
        TBranch        *b__gen_lPt;   
        TBranch        *b__gen_lEta;   
        TBranch        *b__gen_lPhi;   
        TBranch        *b__gen_lE;   
        TBranch        *b__gen_lFlavor;   
        TBranch        *b__gen_lCharge;   
        TBranch        *b__gen_lMomPdg;   
        TBranch        *b__gen_lIsPrompt;   
        TBranch        *b__ttgEventType;
        TBranch        *b__zgEventType;

        TBranch        *b__gen_n;
        TBranch        *b__gen_pt;
        TBranch        *b__gen_eta;
        TBranch        *b__gen_phi;
        TBranch        *b__gen_E;
        TBranch        *b__gen_pdgId;
        TBranch        *b__gen_charge;
        TBranch        *b__gen_status;
        TBranch        *b__gen_isPromptFinalState;
        TBranch        *b__gen_isDirectPromptTauDecayProductFinalState;
        TBranch        *b__gen_isLastCopy;
        TBranch        *b__gen_index;
        TBranch        *b__gen_motherIndex;
        TBranch        *b__gen_daughter_n;
        TBranch        *b__gen_daughterIndex;

        TBranch        *b__passTrigger_e;   
        TBranch        *b__passTrigger_ee;   
        TBranch        *b__passTrigger_eee;   
        TBranch        *b__passTrigger_em;   
        TBranch        *b__passTrigger_m;   
        TBranch        *b__passTrigger_eem;   
        TBranch        *b__passTrigger_mm;   
        TBranch        *b__passTrigger_emm;   
        TBranch        *b__passTrigger_mmm;   
        TBranch        *b__passTrigger_et;
        TBranch        *b__passTrigger_mt;
        TBranch        *b__passTrigger_FR;
        TBranch        *b__passTrigger_FR_iso;
        TBranch        *b__passMETFilters;   
        TBranch        *b__nL;   
        TBranch        *b__nMu;   
        TBranch        *b__nEle;   
        TBranch        *b__nLight;   
        TBranch        *b__nTau;   
        TBranch        *b__lPt;   
        TBranch        *b__lPtCorr;
        TBranch        *b__lEta;   
        TBranch        *b__lEtaSC;   
        TBranch        *b__lPhi;   
        TBranch        *b__lE;   
        TBranch        *b__lECorr;
        TBranch        *b__lFlavor;   
        TBranch        *b__lCharge;   
        TBranch        *b__dxy;   
        TBranch        *b__dz;   
        TBranch        *b__3dIP;   
        TBranch        *b__3dIPSig;   
        TBranch        *b__lElectronSummer16MvaGP;   
        TBranch        *b__lElectronSummer16MvaHZZ;
        TBranch        *b__lElectronMvaFall17Iso;
        TBranch        *b__lElectronMvaFall17NoIso;
        TBranch        *b__lElectronPassMVAFall17NoIsoWPLoose;
        TBranch        *b__lElectronPassMVAFall17NoIsoWP90;
        TBranch        *b__lElectronPassMVAFall17NoIsoWP80;
        TBranch        *b__lElectronPassEmu;   
        TBranch        *b__lElectronPassConvVeto;
        TBranch        *b__lElectronChargeConst;
        TBranch        *b__lElectronMissingHits;
        TBranch        *b__lElectronEInvMinusPInv;
        TBranch        *b__lElectronHOverE;
        TBranch        *b__lElectronSigmaIetaIeta;
        TBranch        *b__leptonMvaTTH;
        TBranch        *b__leptonMvatZq;
        TBranch        *b__leptonMvaTOP;
        TBranch        *b__lPOGVeto;   
        TBranch        *b__lPOGLoose;   
        TBranch        *b__lPOGMedium;   
        TBranch        *b__lPOGTight;   
		TBranch		   *b__tauDecayMode;
        TBranch		   *b__decayModeFinding;
        TBranch		   *b__decayModeFindingNew;
        TBranch		   *b__tauMuonVetoLoose;
        TBranch		   *b__tauMuonVetoTight;
        TBranch		   *b__tauEleVetoVLoose;
        TBranch		   *b__tauEleVetoLoose;
        TBranch		   *b__tauEleVetoMedium;
        TBranch		   *b__tauEleVetoTight;
        TBranch		   *b__tauEleVetoVTight;
        TBranch		   *b__tauPOGVLoose2015;
        TBranch		   *b__tauPOGLoose2015;
        TBranch		   *b__tauPOGMedium2015;
        TBranch		   *b__tauPOGTight2015;
        TBranch		   *b__tauPOGVTight2015;
        TBranch		   *b__tauVLooseMvaNew2015;
        TBranch		   *b__tauLooseMvaNew2015;
        TBranch		   *b__tauMediumMvaNew2015;
        TBranch		   *b__tauTightMvaNew2015;
        TBranch		   *b__tauVTightMvaNew2015;
        TBranch		   *b__tauPOGVVLoose2017v2;
        TBranch		   *b__tauPOGVTight2017v2;
        TBranch		   *b__tauPOGVVTight2017v2;
        TBranch		   *b__tauVLooseMvaNew2017v2;
        TBranch		   *b__tauLooseMvaNew2017v2;
        TBranch		   *b__tauMediumMvaNew2017v2;
        TBranch		   *b__tauTightMvaNew2017v2;
        TBranch		   *b__tauVTightMvaNew2017v2;	
        TBranch        *b__relIso;   
        TBranch        *b__relIso0p4;
        TBranch        *b__relIso0p4MuDeltaBeta;
        TBranch        *b__miniIso;   
        TBranch        *b__miniIsoCharged;   
        TBranch        *b__ptRel;   
        TBranch        *b__ptRatio;   
        TBranch        *b__closestJetCsvV2;
        TBranch        *b__closestJetDeepCsv_b;
        TBranch        *b__closestJetDeepCsv_bb;
        TBranch        *b__closestJetDeepFlavor_b;
        TBranch        *b__closestJetDeepFlavor_bb;
        TBranch        *b__closestJetDeepFlavor_lepb;
        TBranch        *b__selectedTrackMult;   
        TBranch        *b__lMuonSegComp;   
        TBranch        *b__lMuonTrackPt;
        TBranch        *b__lMuonTrackPtErr;
        TBranch        *b__lIsPrompt;   
        TBranch        *b__lMatchPdgId;
        TBranch        *b__lMatchCharge;
        TBranch        *b__lMomPdgId;
        TBranch        *b__lProvenance;
        TBranch        *b__lProvenanceCompressed;
        TBranch        *b__lProvenanceConversion;
        TBranch        *b__nJets;   
        TBranch        *b__jetPt;   
        TBranch        *b__jetPt_JECUp;   
        TBranch        *b__jetPt_JECDown;   
        TBranch        *b__jetSmearedPt;
        TBranch        *b__jetSmearedPt_JECDown;
        TBranch        *b__jetSmearedPt_JECUp;
        TBranch        *b__jetSmearedPt_JERDown;
        TBranch        *b__jetSmearedPt_JERUp;
        TBranch        *b__jetPt_Uncorrected;
        TBranch        *b__jetPt_L1;
        TBranch        *b__jetPt_L2;
        TBranch        *b__jetPt_L3;
        TBranch        *b__jetEta;   
        TBranch        *b__jetPhi;   
        TBranch        *b__jetE;   
        TBranch        *b__jetCsvV2;   
        TBranch        *b__jetDeepCsv_udsg;   
        TBranch        *b__jetDeepCsv_b;   
        TBranch        *b__jetDeepCsv_c;   
        TBranch        *b__jetDeepCsv_bb;   
        TBranch        *b__jetDeepFlavor_b;
        TBranch        *b__jetDeepFlavor_bb;
        TBranch        *b__jetDeepFlavor_lepb;
        TBranch        *b__jetHadronFlavor;   
        TBranch        *b__jetId;   
        TBranch        *b__jetIsTight;
        TBranch        *b__jetIsTightLepVeto;
        TBranch        *b__jetNeutralHadronFraction;
        TBranch        *b__jetChargedHadronFraction;
        TBranch        *b__jetNeutralEmFraction;
        TBranch        *b__jetChargedEmFraction;
        TBranch        *b__jetHFHadronFraction;
        TBranch        *b__jetHFEmFraction;
        TBranch        *b__met;   
        TBranch        *b__metJECDown;   
        TBranch        *b__metJECUp;   
        TBranch        *b__metUnclDown;   
        TBranch        *b__metUnclUp;   
        TBranch        *b__metPhi;   
        TBranch        *b__metPhiJECDown;   
        TBranch        *b__metPhiJECUp;   
        TBranch        *b__metPhiUnclDown;   
        TBranch        *b__metPhiUnclUp;   
        TBranch        *b__metSignificance;
        TBranch        *b__mChi1;
        TBranch        *b__mChi2;

        std::map< std::string, TBranch* > b__triggerMap;
        std::map< std::string, TBranch* > b__MetFilterMap; 
};

#endif
