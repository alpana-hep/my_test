#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PFCalibration/PFChargedHadronAnalyzer/plugins/alpsEdAnalyzer.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
//#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"                                                                                                             
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include <TROOT.h>
#include <TVector3.h>
#include <TH1.h>




using namespace std;
using namespace edm;
using namespace reco;

// define parameters you want to use 
alpsEdAnalyzer::alpsEdAnalyzer(const edm::ParameterSet& iConfig) {
  nCh = std::vector<unsigned int>(10,static_cast<unsigned int>(0));
  nEv = std::vector<unsigned int>(2,static_cast<unsigned int>(0));

  inputTagPFCandidates_
    = iConfig.getParameter<InputTag>("PFCandidates");
  tokenPFCandidates_ = consumes<reco::PFCandidateCollection>(inputTagPFCandidates_);
  
  inputTagPFSimParticles_
    = iConfig.getParameter<InputTag>("PFSimParticles");
  tokenPFSimParticles_ = consumes<reco::PFSimParticleCollection>(inputTagPFSimParticles_);
  
  
  inputTagEcalPFClusters_
    = iConfig.getParameter<InputTag>("EcalPFClusters");
  tokenEcalPFClusters_ = consumes<reco::PFClusterCollection>(inputTagEcalPFClusters_);
  isMBMC_ = iConfig.getUntrackedParameter<bool>("isMinBiasMC",false);
  
  verbose_ =
    iConfig.getUntrackedParameter<bool>("verbose",false);
  
  LogDebug("alpsEdAnalyzer")
    <<" input collection : "<<inputTagPFCandidates_ ;
  
  
  //    edm::Service<TFileService> fs;
  
  // The root tuple                                                                                                                                                       
  outputfile_ = iConfig.getParameter<std::string>("rootOutputFile");
  tf1 = new TFile(outputfile_.c_str(), "RECREATE");
  //tf1 = new TFile("alps.root","recreate");
  //    s->Branch("photon_pt",&photon_pt,"photon_pt");
  
  s = new TTree("s"," PFCalibration");
  
  s->Branch("true",&true_,"true/F");
  //  s->Branch("photon_pt",&photon_pt,"photon_pt");
  s->Branch("photon_pt",&photon_pt);
  s->Branch("photon_eta",&photon_eta);
  s->Branch("photon_phi",&photon_phi);

  s->Branch("neutral_pt",&neutral_pt);
  s->Branch("neutral_eta",&neutral_eta);
  s->Branch("neutral_phi",&neutral_phi);
  s->Branch("charged_pt",&charged_pt);
  s->Branch("charged_eta",&charged_eta);
  s->Branch("charged_phi",&charged_phi);
  s->Branch("electron_pt",&electron_pt); 
  s->Branch("electron_eta",&electron_eta);
  s->Branch("electron_phi",&electron_phi);
  s->Branch("muon_pt",&muon_pt); 
  s->Branch("muon_eta",&muon_eta);
  s->Branch("muon_phi",&muon_phi);
  s->Branch("pfc_pt",&pfc_pt);
  s->Branch("pfc_eta",&pfc_eta);
  s->Branch("pfc_phi",&pfc_phi);
  //  Event *event = new Event();
  //  s->Branch("photon",&photon_,"ptt/F:etaa/F:phii/F");
  //  s->Branch("electron",&electron_,"pt/F:eta:phii");
  //s->Branch("muon",&muon_,"pt/F:eta:phii");
  
  // TH1F* photon_pt = new TH1F("photon_pt","pt of photons",100,200,500);
  // TH1F* photon_eta = new TH1F("photon_eta","eta of photons",10,-10,10);
  // TH1F* photon_phi = new TH1F("photon_phi","phi of photons",10,-10,10);
  // TH1F* electron_pt = new TH1F("electron_pt","pt of electrons",100,200,500);
  // TH1F* electron_eta = new TH1F("electron_eta","eta of electrons",10,-10,10);
  // TH1F* electron_phi = new TH1F("electron_phi","phi of electrons",10,-10,10);
  // TH1F* muon_pt = new TH1F("muon_pt","pt of muons",100,200,500);
  // TH1F* muon_eta = new TH1F("muon_eta","eta of muons",10,-10,10);
  // TH1F* muon_phi = new TH1F("muon_phi","phi of muons",10,-10,10);


}
alpsEdAnalyzer::~alpsEdAnalyzer() {

  std::cout << "Total number of events .............. " << nEv[0] << std::endl;
  std::cout << "Total number of events with atleast 1 sim particle .............. " << nEv[1] << std::endl;
  std::cout<<"total number of pf candidates .... "<<nCh[0]<<std::endl;
  std::cout<<"total number of  photons .... "<<nCh[1]<<std::endl;
  std::cout<<"total number of  neutral hadrons .... "<<nCh[2]<<std::endl;
  std::cout<<"total number of charged hadrons .... "<<nCh[3]<<std::endl;
  std::cout<<"total number of electrons .... "<<nCh[4]<<std::endl;
  std::cout<<"total number of  muons .... "<<nCh[5]<<std::endl;
  

  tf1->cd();
  s->Write();


  tf1->Write();
  tf1->Close();



}
//need to figure out what it is
void
alpsEdAnalyzer::beginRun(const edm::Run& run,
                                  const edm::EventSetup & es) { }


// loop per event
void
alpsEdAnalyzer::analyze(const Event& iEvent,
                                 const EventSetup& iSetup) {

  LogDebug("alpsEdAnalyzer")<<"START event: "<<iEvent.id().event()
				     <<" in run "<<iEvent.id().run()<<endl;




  //get PF candidate
  Handle<PFCandidateCollection> pfCandidates;
  iEvent.getByToken(tokenPFCandidates_, pfCandidates);

  Handle<PFSimParticleCollection> trueParticles;
  bool truesim = iEvent.getByToken(tokenPFSimParticles_, trueParticles);
  //  if(isMBMC_)
  pfcsID.clear(); 

  photon_pt.clear();
  photon_eta.clear();
  photon_phi.clear();
  photon_phi.clear();
  electron_pt.clear();
  electron_eta.clear();
  electron_phi.clear();
  muon_pt.clear();
  muon_eta.clear();
  muon_phi.clear();
  neutral_pt.clear();
  neutral_eta.clear();
  neutral_phi.clear();
  charged_pt.clear();
  charged_eta.clear();
  charged_phi.clear();
  pfc_pt.clear();
  pfc_eta.clear();
  pfc_phi.clear();


  if(isMBMC_)  
  truesim=false;
  if ( truesim ) {
    nEv[0]++;
    if ( (*trueParticles).size() != 1 ) return;
    nEv[1]++;
  
  } 
  // nEv[0]++;

   for(auto sim = pfCandidates->begin(); sim != pfCandidates->end(); ++sim){
     nCh[0]++;
     const reco::PFCandidate& pfc = *sim;
     pfcsID.push_back(pfc.particleId());
     pfc_pt.push_back(pfc.pt());
     pfc_eta.push_back(pfc.eta());
     pfc_phi.push_back(pfc.phi());
     //std::cout<<"pt of  photons .... "<<pfc.pt()<<std::endl;
     //std::cout<<"eta of   photons .... "<<pfc.eta()<<std::endl;

     if ( pfc.particleId() == 4 ) // for photons
     {
       nCh[1]++;
       photon_pt.push_back(pfc.pt());
       photon_eta.push_back(pfc.eta());
       photon_phi.push_back(pfc.phi());
       //       std::cout<<"pt of  photons .... "<<pfc.pt()<<std::endl;
       //       std::cout<<"eta of   photons .... "<<pfc.eta()<<std::endl;
     }
     else if (pfc.particleId() == 5)// for neutral hadrons
     {
       nCh[2]++;
       neutral_pt.push_back(pfc.pt());
       neutral_eta.push_back(pfc.eta());
       neutral_phi.push_back(pfc.phi());
     }
     else if (pfc.particleId() == 1)// for charged hadrons
     {
       nCh[3]++;
       charged_pt.push_back(pfc.pt());
       charged_eta.push_back(pfc.eta());
       charged_phi.push_back(pfc.phi());

     }
     else if (pfc.particleId() == 2)//for electrons
     {
       nCh[4]++;
       electron_pt.push_back(pfc.pt());
       electron_eta.push_back(pfc.eta());
       electron_phi.push_back(pfc.phi());
     }
     else if (pfc.particleId() == 3) // for muons
     {
       nCh[5]++;
       muon_pt.push_back(pfc.pt());
       muon_eta.push_back(pfc.eta());
       muon_phi.push_back(pfc.phi());
     }
   }
     s->Fill();
        return;
   






      
}

  DEFINE_FWK_MODULE(alpsEdAnalyzer);

