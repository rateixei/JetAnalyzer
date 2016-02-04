// -*- C++ -*-
//
// Package:    bbgg/GenAnalyzer
// Class:      GenAnalyzer
// 
/**\class GenAnalyzer GenAnalyzer.cc bbgg/GenAnalyzer/plugins/GenAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rafael Teixeira De Lima
//         Created:  Tue, 01 Dec 2015 09:18:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"


#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <vector>

//
// class declaration
//

class GenAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenAnalyzer(const edm::ParameterSet&);
      ~GenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      typedef math::XYZTLorentzVector LorentzVector;

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      double DeltaR(LorentzVector v1, LorentzVector v2);

      edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken_;
      edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
      edm::EDGetTokenT<edm::View<reco::GenJet> > genjetToken_;

      std::string outFile_;

      TFile* outFile;
      TTree* outTree;

      double pho1_reco_pt;
      double pho1_reco_eta;
      double pho1_reco_phi;
      double pho2_reco_pt;
      double pho2_reco_eta;
      double pho2_reco_phi;
      double pho_reco_mass;
      double pho1_gen_pt;
      double pho1_gen_eta;
      double pho1_gen_phi;
      double pho2_gen_pt;
      double pho2_gen_eta;
      double pho2_gen_phi;
      double pho_gen_mass;
      double bjet1_reco_pt;
      double bjet1_reco_eta;
      double bjet1_reco_phi;
      double bjet2_reco_pt;
      double bjet2_reco_eta;
      double bjet2_reco_phi;
      double bjet_reco_mass;
      double bjet1_gen_pt;
      double bjet1_gen_eta;
      double bjet1_gen_phi;
      double bjet2_gen_pt;
      double bjet2_gen_eta;
      double bjet2_gen_phi;
      double bjet_gen_mass;
      double gjet1_reco_pt;
      double gjet1_reco_eta;
      double gjet1_reco_phi;
      double gjet2_reco_pt;
      double gjet2_reco_eta;
      double gjet2_reco_phi;
      double gjet_reco_mass;
      double gjet1_gen_pt;
      double gjet1_gen_eta;
      double gjet1_gen_phi;
      double gjet2_gen_pt;
      double gjet2_gen_eta;
      double gjet2_gen_phi;
      double gjet_gen_mass;
      double bgjet1_reco_pt;
      double bgjet1_reco_eta;
      double bgjet1_reco_phi;
      double bgjet2_reco_pt;
      double bgjet2_reco_eta;
      double bgjet2_reco_phi;
      double bgjet_reco_mass;
      double bgjet1_gen_pt;
      double bgjet1_gen_eta;
      double bgjet1_gen_phi;
      double bgjet2_gen_pt;
      double bgjet2_gen_eta;
      double bgjet2_gen_phi;
      double bgjet_gen_mass;



      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

double GenAnalyzer::DeltaR(LorentzVector v1, LorentzVector v2)
{
	double DR2 = (v1.eta() - v2.eta())*(v1.eta() - v2.eta()) + (v1.phi() - v2.phi())*(v1.phi() - v2.phi());
	return sqrt(DR2);
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig) :
        photonToken_( consumes<edm::View<pat::Photon> >( iConfig.getParameter<edm::InputTag> ( "photonTag" ) ) ),
        genToken_( consumes<edm::View<reco::GenParticle> >( iConfig.getParameter<edm::InputTag>( "genTag" ) ) ),
        jetToken_( consumes<edm::View<pat::Jet> >( iConfig.getParameter<edm::InputTag>( "jetTag" ) ) ),
        genjetToken_( consumes<edm::View<reco::GenJet> >( iConfig.getParameter<edm::InputTag>( "genjetTag" ) ) )
{

    outFile_ = iConfig.getUntrackedParameter<std::string>( "outFile");

   //now do what ever initialization is needed
   outFile = new TFile(TString(outFile_), "RECREATE");
   outTree = new TTree("genTree", "genTree");

      outTree->Branch("pho1_reco_pt", &pho1_reco_pt, "pho1_reco_pt/D");
      outTree->Branch("pho1_reco_eta", &pho1_reco_eta, "pho1_reco_eta/D");
      outTree->Branch("pho1_reco_phi", &pho1_reco_phi, "pho1_reco_phi/D");
      outTree->Branch("pho2_reco_pt", &pho2_reco_pt, "pho2_reco_pt/D");
      outTree->Branch("pho2_reco_eta", &pho2_reco_eta, "pho2_reco_eta/D");
      outTree->Branch("pho2_reco_phi", &pho2_reco_phi, "pho2_reco_phi/D");
      outTree->Branch("pho_reco_mass", &pho_reco_mass, "pho_reco_mass/D");
      outTree->Branch("pho1_gen_pt", &pho1_gen_pt, "pho1_gen_pt/D");
      outTree->Branch("pho1_gen_eta", &pho1_gen_eta, "pho1_gen_eta/D");
      outTree->Branch("pho1_gen_phi", &pho1_gen_phi, "pho1_gen_phi/D");
      outTree->Branch("pho2_gen_pt", &pho2_gen_pt, "pho2_gen_pt/D");
      outTree->Branch("pho2_gen_eta", &pho2_gen_eta, "pho2_gen_eta/D");
      outTree->Branch("pho2_gen_phi", &pho2_gen_phi, "pho2_gen_phi/D");
      outTree->Branch("pho_gen_mass", &pho_gen_mass, "pho_gen_mass/D");
      outTree->Branch("bjet1_reco_pt", &bjet1_reco_pt, "bjet1_reco_pt/D");
      outTree->Branch("bjet1_reco_eta", &bjet1_reco_eta, "bjet1_reco_eta/D");
      outTree->Branch("bjet1_reco_phi", &bjet1_reco_phi, "bjet1_reco_phi/D");
      outTree->Branch("bjet2_reco_pt", &bjet2_reco_pt, "bjet2_reco_pt/D");
      outTree->Branch("bjet2_reco_eta", &bjet2_reco_eta, "bjet2_reco_eta/D");
      outTree->Branch("bjet2_reco_phi", &bjet2_reco_phi, "bjet2_reco_phi/D");
      outTree->Branch("bjet_reco_mass", &bjet_reco_mass, "bjet_reco_mass/D");
      outTree->Branch("bjet1_gen_pt", &bjet1_gen_pt, "bjet1_gen_pt/D");
      outTree->Branch("bjet1_gen_eta", &bjet1_gen_eta, "bjet1_gen_eta/D");
      outTree->Branch("bjet1_gen_phi", &bjet1_gen_phi, "bjet1_gen_phi/D");
      outTree->Branch("bjet2_gen_pt", &bjet2_gen_pt, "bjet2_gen_pt/D");
      outTree->Branch("bjet2_gen_eta", &bjet2_gen_eta, "bjet2_gen_eta/D");
      outTree->Branch("bjet2_gen_phi", &bjet2_gen_phi, "bjet2_gen_phi/D");
      outTree->Branch("bjet_gen_mass", &bjet_gen_mass, "bjet_gen_mass/D");
      outTree->Branch("gjet1_reco_pt", &gjet1_reco_pt, "gjet1_reco_pt/D");
      outTree->Branch("gjet1_reco_eta", &gjet1_reco_eta, "gjet1_reco_eta/D");
      outTree->Branch("gjet1_reco_phi", &gjet1_reco_phi, "gjet1_reco_phi/D");
      outTree->Branch("gjet2_reco_pt", &gjet2_reco_pt, "gjet2_reco_pt/D");
      outTree->Branch("gjet2_reco_eta", &gjet2_reco_eta, "gjet2_reco_eta/D");
      outTree->Branch("gjet2_reco_phi", &gjet2_reco_phi, "gjet2_reco_phi/D");
      outTree->Branch("gjet_reco_mass", &gjet_reco_mass, "gjet_reco_mass/D");
      outTree->Branch("gjet1_gen_pt", &gjet1_gen_pt, "gjet1_gen_pt/D");
      outTree->Branch("gjet1_gen_eta", &gjet1_gen_eta, "gjet1_gen_eta/D");
      outTree->Branch("gjet1_gen_phi", &gjet1_gen_phi, "gjet1_gen_phi/D");
      outTree->Branch("gjet2_gen_pt", &gjet2_gen_pt, "gjet2_gen_pt/D");
      outTree->Branch("gjet2_gen_eta", &gjet2_gen_eta, "gjet2_gen_eta/D");
      outTree->Branch("gjet2_gen_phi", &gjet2_gen_phi, "gjet2_gen_phi/D");
      outTree->Branch("gjet_gen_mass", &gjet_gen_mass, "gjet_gen_mass/D");
      outTree->Branch("bgjet1_reco_pt", &bgjet1_reco_pt, "bgjet1_reco_pt/D");
      outTree->Branch("bgjet1_reco_eta", &bgjet1_reco_eta, "bgjet1_reco_eta/D");
      outTree->Branch("bgjet1_reco_phi", &bgjet1_reco_phi, "bgjet1_reco_phi/D");
      outTree->Branch("bgjet2_reco_pt", &bgjet2_reco_pt, "bgjet2_reco_pt/D");
      outTree->Branch("bgjet2_reco_eta", &bgjet2_reco_eta, "bgjet2_reco_eta/D");
      outTree->Branch("bgjet2_reco_phi", &bgjet2_reco_phi, "bgjet2_reco_phi/D");
      outTree->Branch("bgjet_reco_mass", &bgjet_reco_mass, "bgjet_reco_mass/D");
      outTree->Branch("bgjet1_gen_pt", &bgjet1_gen_pt, "bgjet1_gen_pt/D");
      outTree->Branch("bgjet1_gen_eta", &bgjet1_gen_eta, "bgjet1_gen_eta/D");
      outTree->Branch("bgjet1_gen_phi", &bgjet1_gen_phi, "bgjet1_gen_phi/D");
      outTree->Branch("bgjet2_gen_pt", &bgjet2_gen_pt, "bgjet2_gen_pt/D");
      outTree->Branch("bgjet2_gen_eta", &bgjet2_gen_eta, "bgjet2_gen_eta/D");
      outTree->Branch("bgjet2_gen_phi", &bgjet2_gen_phi, "bgjet2_gen_phi/D");
      outTree->Branch("bgjet_gen_mass", &bgjet_gen_mass, "bgjet_gen_mass/D");

}


GenAnalyzer::~GenAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

      pho1_reco_pt = -1;
      pho1_reco_eta = -1;
      pho1_reco_phi = -1;
      pho2_reco_pt = -1;
      pho2_reco_eta = -1;
      pho2_reco_phi = -1;
      pho_reco_mass = -1;
      pho1_gen_pt = -1;
      pho1_gen_eta = -1;
      pho1_gen_phi = -1;
      pho2_gen_pt = -1;
      pho2_gen_eta = -1;
      pho2_gen_phi = -1;
      pho_gen_mass = -1;
      bjet1_reco_pt = -1;
      bjet1_reco_eta = -1;
      bjet1_reco_phi = -1;
      bjet2_reco_pt = -1;
      bjet2_reco_eta = -1;
      bjet2_reco_phi = -1;
      bjet_reco_mass = -1;
      bjet1_gen_pt = -1;
      bjet1_gen_eta = -1;
      bjet1_gen_phi = -1;
      bjet2_gen_pt = -1;
      bjet2_gen_eta = -1;
      bjet2_gen_phi = -1;
      bjet_gen_mass = -1;
      gjet1_reco_pt = -1;
      gjet1_reco_eta = -1;
      gjet1_reco_phi = -1;
      gjet2_reco_pt = -1;
      gjet2_reco_eta = -1;
      gjet2_reco_phi = -1;
      gjet_reco_mass = -1;
      gjet1_gen_pt = -1;
      gjet1_gen_eta = -1;
      gjet1_gen_phi = -1;
      gjet2_gen_pt = -1;
      gjet2_gen_eta = -1;
      gjet2_gen_phi = -1;
      gjet_gen_mass = -1;
      bgjet1_reco_pt = -1;
      bgjet1_reco_eta = -1;
      bgjet1_reco_phi = -1;
      bgjet2_reco_pt = -1;
      bgjet2_reco_eta = -1;
      bgjet2_reco_phi = -1;
      bgjet_reco_mass = -1;
      bgjet1_gen_pt = -1;
      bgjet1_gen_eta = -1;
      bgjet1_gen_phi = -1;
      bgjet2_gen_pt = -1;
      bgjet2_gen_eta = -1;
      bgjet2_gen_phi = -1;
      bgjet_gen_mass = -1;

   using namespace edm;
   using namespace std;

	Handle<View<pat::Photon> > photons;
        iEvent.getByToken( photonToken_, photons );

	Handle<View<reco::GenParticle> > gens;
        iEvent.getByToken( genToken_, gens );

	Handle<View<pat::Jet> > jets;
        iEvent.getByToken( jetToken_, jets );

	Handle<View<reco::GenJet> > genJets;
        iEvent.getByToken( genjetToken_, genJets );
//        iEvent.getByLabel( "genJets", genJets );


	vector<LorentzVector> genPhotons;
	vector<LorentzVector> genBs;
	vector<LorentzVector> genMatchedRecoPhotons;
	vector<unsigned int> genMatchedRecoPhotons_genIndex;
	vector<LorentzVector> genMatchedRecoJets;
	vector<unsigned int> genMatchedRecoJets_genIndex;
	vector<LorentzVector> genMatchedGenJets;
	vector<unsigned int> genMatchedGenJets_genIndex;
	vector<LorentzVector> genJetsMatchedRecoJets;
	vector<unsigned int> genJetsMatchedRecoJets_genIndex;
	

	for(unsigned int i = 0; i < gens->size(); i++)
	{
		Ptr<reco::GenParticle> gen = gens->ptrAt(i);
		int motherId = (gen->mother()) ? gen->mother()->pdgId() : -1;
//		int mothermotherId = -1;
//		if(motherId)
//			mothermotherId = (gen->mother()->mother()) ? gen->mother()->mother()->pdgId() : -1;

//		if(gen->pdgId() == 5 || gen->pdgId() == -5)
//			std::cout << "Found b jet! Mother: " << motherId << " MotherMother: " << mothermotherId  << " Status: " << gen->status() << std::endl;

		if(motherId != 25)
			continue;
		else{
			if(gen->pdgId() == 22) genPhotons.push_back(gen->p4());
			if(gen->pdgId() == 5) genBs.push_back(gen->p4());
			if(gen->pdgId() == -5) genBs.push_back(gen->p4());
		}
	}

//	std::cout << "Number of b's: " << genBs.size() << std::endl;

	for(unsigned int i = 0; i < genPhotons.size(); i++)
	{
		for(unsigned int p = 0; p < photons->size(); p++)
		{
			Ptr<pat::Photon> pho = photons->ptrAt(p);
			LorentzVector v_pho = pho->p4();
			if( GenAnalyzer::DeltaR(v_pho, genPhotons[i]) < 0.3){
				genMatchedRecoPhotons.push_back(v_pho);
				genMatchedRecoPhotons_genIndex.push_back(i);
				continue;
			}
		}
	}

	if(genMatchedRecoPhotons.size() == 1){
		pho1_reco_pt = genMatchedRecoPhotons[0].pt();
		pho1_reco_eta = genMatchedRecoPhotons[0].eta();
		pho1_reco_phi = genMatchedRecoPhotons[0].phi();
		pho1_gen_pt = genPhotons[ genMatchedRecoPhotons_genIndex[0] ].pt();
		pho1_gen_eta = genPhotons[ genMatchedRecoPhotons_genIndex[0] ].eta();
		pho1_gen_phi = genPhotons[ genMatchedRecoPhotons_genIndex[0] ].phi();
	}	
	if(genMatchedRecoPhotons.size() == 2){
		unsigned int p1 = ( genMatchedRecoPhotons[0].pt() > genMatchedRecoPhotons[1].pt()) ? 0 : 1;
		unsigned int p2 = ( genMatchedRecoPhotons[0].pt() > genMatchedRecoPhotons[1].pt()) ? 1 : 0;

		pho1_reco_pt = genMatchedRecoPhotons[p1].pt();
		pho1_reco_eta = genMatchedRecoPhotons[p1].eta();
		pho1_reco_phi = genMatchedRecoPhotons[p1].phi();
		pho1_gen_pt = genPhotons[ genMatchedRecoPhotons_genIndex[p1] ].pt();
		pho1_gen_eta = genPhotons[ genMatchedRecoPhotons_genIndex[p1] ].eta();
		pho1_gen_phi = genPhotons[ genMatchedRecoPhotons_genIndex[p1] ].phi();
		pho2_reco_pt = genMatchedRecoPhotons[p2].pt();
		pho2_reco_eta = genMatchedRecoPhotons[p2].eta();
		pho2_reco_phi = genMatchedRecoPhotons[p2].phi();
		pho2_gen_pt = genPhotons[ genMatchedRecoPhotons_genIndex[p2] ].pt();
		pho2_gen_eta = genPhotons[ genMatchedRecoPhotons_genIndex[p2] ].eta();
		pho2_gen_phi = genPhotons[ genMatchedRecoPhotons_genIndex[p2] ].phi();
		pho_reco_mass = ( genMatchedRecoPhotons[p1] + genMatchedRecoPhotons[p2]).M();
		pho_gen_mass = ( genPhotons[ genMatchedRecoPhotons_genIndex[p2] ] + genPhotons[ genMatchedRecoPhotons_genIndex[p1] ] ).M();
	}	


	for(unsigned int i = 0; i < genBs.size(); i++)
	{
		for(unsigned int p = 0; p < jets->size(); p++)
		{
			Ptr<pat::Jet> jet = jets->ptrAt(p);

            std::vector<std::string> userFloats = jet->userFloatNames();
            for (unsigned int uf = 0; uf < userFloats.size(); uf++)
                std::cout << " ## " << userFloats[uf] << " = " << jet->userFloat(userFloats[uf]) << " // ";

			LorentzVector v_jet = jet->p4();
			if( GenAnalyzer::DeltaR(v_jet, genBs[i]) < 0.3){
				genMatchedRecoJets.push_back(v_jet);
				genMatchedRecoJets_genIndex.push_back(i);
				continue;
			}
		}
	}

        if(genMatchedRecoJets.size() == 1){
                bjet1_reco_pt = genMatchedRecoJets[0].pt();
                bjet1_reco_eta = genMatchedRecoJets[0].eta();
                bjet1_reco_phi = genMatchedRecoJets[0].phi();
                bjet1_gen_pt = genBs[ genMatchedRecoJets_genIndex[0] ].pt();
                bjet1_gen_eta = genBs[ genMatchedRecoJets_genIndex[0] ].eta();
                bjet1_gen_phi = genBs[ genMatchedRecoJets_genIndex[0] ].phi();
        }
        if(genMatchedRecoJets.size() == 2){
                unsigned int p1 = ( genMatchedRecoJets[0].pt() > genMatchedRecoJets[1].pt()) ? 0 : 1;
                unsigned int p2 = ( genMatchedRecoJets[0].pt() > genMatchedRecoJets[1].pt()) ? 1 : 0;

                bjet1_reco_pt = genMatchedRecoJets[p1].pt();
                bjet1_reco_eta = genMatchedRecoJets[p1].eta();
                bjet1_reco_phi = genMatchedRecoJets[p1].phi();
                bjet1_gen_pt = genBs[ genMatchedRecoJets_genIndex[p1] ].pt();
                bjet1_gen_eta = genBs[ genMatchedRecoJets_genIndex[p1] ].eta();
                bjet1_gen_phi = genBs[ genMatchedRecoJets_genIndex[p1] ].phi();
                bjet2_reco_pt = genMatchedRecoJets[p2].pt();
                bjet2_reco_eta = genMatchedRecoJets[p2].eta();
                bjet2_reco_phi = genMatchedRecoJets[p2].phi();
                bjet2_gen_pt = genBs[ genMatchedRecoJets_genIndex[p2] ].pt();
                bjet2_gen_eta = genBs[ genMatchedRecoJets_genIndex[p2] ].eta();
                bjet2_gen_phi = genBs[ genMatchedRecoJets_genIndex[p2] ].phi();
                bjet_reco_mass = ( genMatchedRecoJets[p1] + genMatchedRecoJets[p2]).M();
                bjet_gen_mass = ( genBs[ genMatchedRecoJets_genIndex[p2] ] + genBs[ genMatchedRecoJets_genIndex[p1] ] ).M();
        }



	for(unsigned int i = 0; i < genBs.size(); i++)
	{
		for(unsigned int p = 0; p < genJets->size(); p++)
		{
			Ptr<reco::GenJet> gjet = genJets->ptrAt(p);
			LorentzVector v_gjet = gjet->p4();
			if( GenAnalyzer::DeltaR(v_gjet, genBs[i]) < 0.3){
				genMatchedGenJets.push_back(v_gjet);
				genMatchedGenJets_genIndex.push_back(i);
				continue;
			}
		}
	}


        if(genMatchedGenJets.size() == 1){
                gjet1_reco_pt = genMatchedGenJets[0].pt();
                gjet1_reco_eta = genMatchedGenJets[0].eta();
                gjet1_reco_phi = genMatchedGenJets[0].phi();
                gjet1_gen_pt = genBs[ genMatchedGenJets_genIndex[0] ].pt();
                gjet1_gen_eta = genBs[ genMatchedGenJets_genIndex[0] ].eta();
                gjet1_gen_phi = genBs[ genMatchedGenJets_genIndex[0] ].phi();
        }
        if(genMatchedGenJets.size() == 2){
                unsigned int p1 = ( genMatchedGenJets[0].pt() > genMatchedGenJets[1].pt()) ? 0 : 1;
                unsigned int p2 = ( genMatchedGenJets[0].pt() > genMatchedGenJets[1].pt()) ? 1 : 0;

                gjet1_reco_pt = genMatchedGenJets[p1].pt();
                gjet1_reco_eta = genMatchedGenJets[p1].eta();
                gjet1_reco_phi = genMatchedGenJets[p1].phi();
                gjet1_gen_pt = genBs[ genMatchedGenJets_genIndex[p1] ].pt();
                gjet1_gen_eta = genBs[ genMatchedGenJets_genIndex[p1] ].eta();
                gjet1_gen_phi = genBs[ genMatchedGenJets_genIndex[p1] ].phi();
                gjet2_reco_pt = genMatchedGenJets[p2].pt();
                gjet2_reco_eta = genMatchedGenJets[p2].eta();
                gjet2_reco_phi = genMatchedGenJets[p2].phi();
                gjet2_gen_pt = genBs[ genMatchedGenJets_genIndex[p2] ].pt();
                gjet2_gen_eta = genBs[ genMatchedGenJets_genIndex[p2] ].eta();
                gjet2_gen_phi = genBs[ genMatchedGenJets_genIndex[p2] ].phi();
                gjet_reco_mass = ( genMatchedGenJets[p1] + genMatchedGenJets[p2]).M();
                gjet_gen_mass = ( genBs[ genMatchedGenJets_genIndex[p2] ] + genBs[ genMatchedGenJets_genIndex[p1] ] ).M();
        }


	for(unsigned int i = 0; i < genMatchedGenJets.size(); i++)
	{
		for(unsigned int p = 0; p < jets->size(); p++)
		{
			Ptr<pat::Jet> jet = jets->ptrAt(p);
			LorentzVector v_jet = jet->p4();
			if( GenAnalyzer::DeltaR(v_jet, genMatchedGenJets[i]) < 0.3){
				genJetsMatchedRecoJets.push_back(v_jet);
				genJetsMatchedRecoJets_genIndex.push_back(i);
				continue;
			}
		}

	}

        if(genJetsMatchedRecoJets.size() == 1){
                bgjet1_reco_pt = genJetsMatchedRecoJets[0].pt();
                bgjet1_reco_eta = genJetsMatchedRecoJets[0].eta();
                bgjet1_reco_phi = genJetsMatchedRecoJets[0].phi();
                bgjet1_gen_pt = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[0] ].pt();
                bgjet1_gen_eta = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[0] ].eta();
                bgjet1_gen_phi = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[0] ].phi();
        }
        if(genJetsMatchedRecoJets.size() == 2){
                unsigned int p1 = ( genJetsMatchedRecoJets[0].pt() > genJetsMatchedRecoJets[1].pt()) ? 0 : 1;
                unsigned int p2 = ( genJetsMatchedRecoJets[0].pt() > genJetsMatchedRecoJets[1].pt()) ? 1 : 0;

                bgjet1_reco_pt = genJetsMatchedRecoJets[p1].pt();
                bgjet1_reco_eta = genJetsMatchedRecoJets[p1].eta();
                bgjet1_reco_phi = genJetsMatchedRecoJets[p1].phi();
                bgjet1_gen_pt = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[p1] ].pt();
                bgjet1_gen_eta = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[p1] ].eta();
                bgjet1_gen_phi = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[p1] ].phi();
                bgjet2_reco_pt = genJetsMatchedRecoJets[p2].pt();
                bgjet2_reco_eta = genJetsMatchedRecoJets[p2].eta();
                bgjet2_reco_phi = genJetsMatchedRecoJets[p2].phi();
                bgjet2_gen_pt = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[p2] ].pt();
                bgjet2_gen_eta = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[p2] ].eta();
                bgjet2_gen_phi = genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[p2] ].phi();
                bgjet_reco_mass = ( genJetsMatchedRecoJets[p1] + genJetsMatchedRecoJets[p2]).M();
                bgjet_gen_mass = ( genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[p2] ] + genMatchedGenJets[ genJetsMatchedRecoJets_genIndex[p1] ] ).M();
        }

	outTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
GenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenAnalyzer::endJob() 
{

	outFile->cd();
	outTree->Write();
	outFile->Write();
	outFile->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
GenAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
GenAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
GenAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
GenAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
