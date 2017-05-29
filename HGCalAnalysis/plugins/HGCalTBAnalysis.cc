// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"



#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
//#include "RecoLocalCalo/HGCalRecAlgos/src/HGCalDepthPreClusterer.cc"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

//#include "RecoNtuples/HGCalTBAnalysis/interface/AEvent.h"
//#include "RecoNtuples/HGCalTBAnalysis/interface/AObData.h"

#include "RecoNtuples/HGCalAnalysis/interface/UsefulClasses.h"
//#include "../utils/UsefulClasses.cpp"

#include <string>
#include <map>

#include "TH3F.h"

class HGCalTBAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalTBAnalysis(const edm::ParameterSet&);
  ~HGCalTBAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<reco::CaloClusterCollection> _clusters;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > _multiClusters;

  std::string                detector;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  bool                       rawRecHits;
  hgcal::RecHitTools         recHitTools;

  //  std::unique_ptr<hgcal::ClusterTools> clusterTools;

  //  int CMSSW_cellLayer_1[8] = {10, 11, 12, 13, 14, 16, 17, 19};   //cfg1 
  int CMSSW_cellLayer[8] = {8, 12, 16, 19, 22, 23, 25, 28};   // cfg2


  TH2F* h3_bary_yz;
  TH2F* h3_bary_xz;
  TH2F* h3_seed_yz;
  TH2F* h3_seed_xz;
  TH2F* h3_2dCl_yz;
  TH2F* h3_2dCl_xz;
  //  TH3* h3_2dCl_xyz;


  TH2F* h2_YvsX_rh[8];
  TH2F* h2_YvsX_2dCl[8];
  TH2F* h2_YvsX_rh2dCl[8];


  TH3F* ZvsYvsX_rh;
  TH3F* ZvsYvsX_2Drh;
  TH3F* ZvsYvsX_2Dcl;

  std::vector<int> nHits_layer;
  std::vector<int> nHitsWithTime_layer;

  TH2F* h_EnergyMC_vsN2D;
  TH2F* h_EnergyFracMC_vsN2D;
  TH1F* h_EfracMCl_wrt_Best;

  // TH2F* h2_dR_SeedCl2_multiCl_vsEta;
  // TProfile* tp_dR_SeedCl2_multiCl_vsEta;
  // TH1F* h_dR_SeedCl2_multiCl;

  TH2F* h_nMCvsnum2dClBestMC;
  TH1F* h_n2dClPerMC;
  TH1F* h_n2dClPerAllMC;
  TH1F* h_nMCPerEvt;

  TH1F* h_n2dClinMCPerEvt;
  TH1F* h_n2dClPerEvt;

  TProfile* n2DClinM_vsL;
  TProfile* n2DCl_vsL;

  TH1F* h_dR_2DCl_multiCl;
  TH1F* h_dR_2DCl;
  //TH1F* h_2DClSumE_overMCE;
  TH1F* h_2DClSumE;
  TProfile* tp_2DClSumE_vsLayer;
  TProfile* tp_2DClSumE_vsRadius_wrtMCl;


  TH1F* h_dR_2DClAllM;
  TH1F* h_2DClAllMSumE;
  TProfile* tp_2DClAllMSumE_vsLayer;

  TH1F* h_dR_2DClAllSelM;
  TH1F* h_2DClAllSelMSumE;
  TProfile* tp_2DClAllSelMSumE_vsLayer;

  TH1F* h_dR_2DClAll;
  TH1F* h_2DClAllSumE;
  TProfile* tp_2DClAllSumE_vsLayer;

  TH1F* h_dR_recHitAll;
  TH1F* h_recHitAllSumE;
  TProfile* tp_recHitAllSumE_vsLayer;

  TH1F* h_mtcRHSumE;
  TProfile* tp_mtcRHSumE_vsLayer;


  TH2F* h2_dPhivsdEta_rhAxis;



};  



HGCalTBAnalysis::HGCalTBAnalysis(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  rawRecHits(iConfig.getParameter<bool>("rawRecHits"))
{

  //  std::cout << " >>> costruttore " << std::endl;

  //now do what ever initialization is needed
  usesResource("TFileService");

  if(detector=="EE") {
    _recHitsEE = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"));
    algo = 2;
  }
  _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("hgcalLayerClusters"));
  _multiClusters = consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("hgcalLayerClusters"));


  //  auto sumes = consumesCollector();
  //  clusterTools = std::make_unique<hgcal::ClusterTools>(iConfig,sumes);

  edm::Service<TFileService> fs;

  for(int ij=0; ij<8; ++ij){
    h2_YvsX_rh[ij] = fs->make<TH2F>(Form("h2_YvsX_rh_layer%d", ij+1), "", 100, -100., -85., 100, 85., 105.);
    h2_YvsX_2dCl[ij] = fs->make<TH2F>(Form("h2_YvsX_2dCl_layer%d", ij+1), "", 100, -100., -85., 100, 85., 105.);
    h2_YvsX_rh2dCl[ij] = fs->make<TH2F>(Form("h2_YvsX_rh2dCl_layer%d", ij+1), "", 100, -100., -85., 100, 85., 105.);
  }


  ZvsYvsX_rh = fs->make<TH3F>("ZvsYvsX_rh", " ", 10, 0., 10., 100, -100., -85., 100, 85., 105.);
  ZvsYvsX_2Drh = fs->make<TH3F>("ZvsYvsX_2Drh", " ", 10, 0., 10., 100, -100., -85., 100, 85., 105.);
  ZvsYvsX_2Dcl = fs->make<TH3F>("ZvsYvsX_2Dcl", " ", 10, 0., 10., 100, -100., -85., 100, 85., 105.);

  h3_bary_xz = fs->make<TH2F>("h3_bary_xz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_bary_yz = fs->make<TH2F>("h3_bary_yz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_seed_xz = fs->make<TH2F>("h3_seed_xz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_seed_yz = fs->make<TH2F>("h3_seed_yz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_2dCl_xz = fs->make<TH2F>("h3_2dCl_xz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_2dCl_yz = fs->make<TH2F>("h3_2dCl_yz", "",  1000, 250., 650., 1000, -100., 100.);
  //  h3_2dCl_xyz = fs->make<TH3>("h3_2dCl_xyz", "",  1000, 200., 500., 1000, -100., 100., 1000, -100., 100.);

  h_EfracMCl_wrt_Best = fs->make<TH1F>("h_EfracMCl_wrt_Best", "", 500, 0., 1.1);
  h_EnergyFracMC_vsN2D = fs->make<TH2F>("h_EnergyFracMC_vsN2D", "", 50, 0., 50., 500, 0., 1.1);
  h_EnergyMC_vsN2D = fs->make<TH2F>("h_EnergyMC_vsN2D", "", 50, 0., 50., 300, 0., 150.);

  //  h2_dR_SeedCl2_multiCl_vsEta = fs->make<TH2F>("h2_dR_SeedCl2_multiCl_vsEta", "", 500, 1.5, 3.5, 1000, 0., 5.);
  //  tp_dR_SeedCl2_multiCl_vsEta = fs->make<TProfile>("tp_dR_SeedCl2_multiCl_vsEta", "", 500, 1.5, 3.5);
  //  h_dR_SeedCl2_multiCl = fs->make<TH1F>("h_dR_SeedCl2_multiCl", "", 1000, 0., 3.);

  h_nMCvsnum2dClBestMC = fs->make<TH2F>("h_nMCvsnum2dClBestMC", "", 50, 0., 50., 100, 0., 100.);
  h_n2dClPerMC = fs->make<TH1F>("h_n2dClPerMC", "", 50, 0., 50.);
  h_n2dClPerAllMC = fs->make<TH1F>("h_n2dClPerAllMC", "", 50, 0., 50.);
  h_nMCPerEvt = fs->make<TH1F>("h_nMCPerEvt", "", 100, 0., 100.);

  h_n2dClinMCPerEvt = fs->make<TH1F>("h_n2dClinMCPerEvt", "", 100, 0., 100.);
  h_n2dClPerEvt = fs->make<TH1F>("h_n2dClPerEvt", "", 100, 0., 100.);

  n2DClinM_vsL = fs->make<TProfile>("n2DClinM_vsL", "", 30, 0., 30.);
  n2DCl_vsL = fs->make<TProfile>("n2DCl_vsL", "", 30, 0., 30.);


  h_dR_2DCl_multiCl = fs->make<TH1F>("h_dR_2DCl_multiCl", "", 1000, 0., 3.);
  h_dR_2DCl = fs->make<TH1F>("h_dR_2DCl", "", 1500, 0., 15.);
  //  h_2DClSumE_overMCE = fs->make<TH1F>("h_2DClSumE_overMCE", "", 500, 0., 2.);  

  h_2DClSumE = fs->make<TH1F>("h_2DClSumE", "", 500, 0., 250.);
  tp_2DClSumE_vsLayer = fs->make<TProfile>("tp_2DClSumE_vsLayer", "", 1000, 0., 30);
  tp_2DClSumE_vsRadius_wrtMCl = fs->make<TProfile>("tp_2DClSumE_vsRadius_wrtMCl", "", 5, 0., 5.);

  h_dR_2DClAllM = fs->make<TH1F>("h_dR_2DClAllM", "", 1500, 0., 15.);
  h_2DClAllMSumE = fs->make<TH1F>("h_2DClAllMSumE", "", 500, 0., 250.);
  tp_2DClAllMSumE_vsLayer = fs->make<TProfile>("tp_2DClAllMSumE_vsLayer", "", 1000, 0., 30);

  h_dR_recHitAll = fs->make<TH1F>("h_dR_recHitAll", "", 1500, 0., 15.);
  h_recHitAllSumE = fs->make<TH1F>("h_recHitAllSumE", "", 500, 0., 250.);
  tp_recHitAllSumE_vsLayer = fs->make<TProfile>("tp_recHitAllSumE_vsLayer", "", 1000, 0., 30);

  h_dR_2DClAllSelM = fs->make<TH1F>("h_dR_2DClAllSelM", "", 1500, 0., 15.);
  h_2DClAllSelMSumE = fs->make<TH1F>("h_2DClAllSelMSumE", "", 500, 0., 250.);
  tp_2DClAllSelMSumE_vsLayer = fs->make<TProfile>("tp_2DClAllSelMSumE_vsLayer", "", 1000, 0., 30);

  h_dR_2DClAll = fs->make<TH1F>("h_dR_2DClAll", "", 1500, 0., 15.);
  h_2DClAllSumE = fs->make<TH1F>("h_2DClAllSumE", "", 500, 0., 250.);
  tp_2DClAllSumE_vsLayer = fs->make<TProfile>("tp_2DClAllSumE_vsLayer", "", 1000, 0., 30);

  //sel

  h_mtcRHSumE = fs->make<TH1F>("h_mtcRHSumE", "", 500, 0., 250.);
  tp_mtcRHSumE_vsLayer = fs->make<TProfile>("tp_mtcRHSumE_vsLayer", "", 1000, 0., 30);


  h2_dPhivsdEta_rhAxis = fs->make<TH2F>("h2_dPhivsdEta_rhAxis", "", 1000, -2., 2., 1000, -2., 2.);
}

HGCalTBAnalysis::~HGCalTBAnalysis()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
HGCalTBAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << " >>> analyzer " << std::endl;
  using namespace edm;

  recHitTools.getEventSetup(iSetup);

  Handle<HGCRecHitCollection> recHitHandleEE;

  Handle<reco::CaloClusterCollection> clusterHandle;
  iEvent.getByToken(_clusters,clusterHandle);
  const reco::CaloClusterCollection &clusters = *clusterHandle;

  Handle<std::vector<reco::HGCalMultiCluster> > multiClusterHandle;
  iEvent.getByToken(_multiClusters, multiClusterHandle);
  const std::vector<reco::HGCalMultiCluster>& multiClusters = *multiClusterHandle;


  //make a map detid-rechit
  std::map<DetId,const HGCRecHit*> hitmap;
  switch(algo){
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      for(unsigned int i = 0; i < rechitsEE.size(); i++){
	hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      break;
    }
  default:
    break;
  }

  const HGCRecHitCollection& rechitsEE = *recHitHandleEE;

  // std::cout << " >>> hitmap.size() =  " << hitmap.size() << std::endl;
  // for(std::map<DetId,const HGCRecHit*>::iterator iop=hitmap.begin(); iop != hitmap.end(); ++iop){
  //   std::cout <<" energy " << iop. << " energy = " << iop->second << std::endl; 
  // }

  ////////////////////

  //  std::cout << " >>> now caloparticles " << std::endl;
  /*
  float seedEnergy = 0;
  float scondEnergy = 0;
  unsigned int seedIndex = 0;
  unsigned int secondIndex = 0;
  for(unsigned int i = 0; i < clusters.size(); i++){
    if(clusters.at(i).eta() < 0) continue;
    h3_2dCl_xz->Fill(clusters.at(i).z(), clusters.at(i).x(),  clusters.at(i).energy());
    h3_2dCl_yz->Fill(clusters.at(i).z(), clusters.at(i).y(), clusters.at(i).energy());
    //    h3_2dCl_xyz->Fill(clusters.at(i).z(), clusters.at(i).y(), clusters.at(i).x(), clusters.at(i).energy());
    if(clusters.at(i).energy() > seedEnergy){
      seedEnergy = clusters.at(i).energy();
      seedIndex = i;
    }
  }
    
  for(unsigned int i = 0; i < clusters.size(); i++){
    if(clusters.at(i).eta() < 0) continue;
    if(clusters.at(i).energy() > scondEnergy && i != seedIndex){
      scondEnergy = clusters.at(i).energy();
      secondIndex = i;
    }
  }

  h3_bary_xz->Fill(clusters.at(secondIndex).z(), clusters.at(secondIndex).x());
  h3_bary_yz->Fill(clusters.at(secondIndex).z(), clusters.at(secondIndex).y());
  h3_bary_xz->Fill(clusters.at(seedIndex).z(), clusters.at(seedIndex).x());
  h3_bary_yz->Fill(clusters.at(seedIndex).z(), clusters.at(seedIndex).y());
  h3_seed_xz->Fill(clusters.at(seedIndex).z(), clusters.at(seedIndex).x());
  h3_seed_yz->Fill(clusters.at(seedIndex).z(), clusters.at(seedIndex).y());
  h3_seed_xz->Fill(0., 0.);
  h3_seed_yz->Fill(0., 0.);
  */
  ////////////////


  std::array<double,3> vtx{ {-92.78, 96.42, 0.}};
      
  int nmclusCount = 0;
  int num2DClinBestMC = 0;
  int num2DClinMC = 0;
  int num2DClinMC_vsL[30];
  for(int ij=0; ij<30; ++ij) num2DClinMC_vsL[ij] = 0;
  unsigned int MC_best_index = 0;
  float MC_best_energy = 0;
  //    std::cout << " >>> now MCL " << std::endl;

  UsefulClasses utilsMet = UsefulClasses(0, 0); 

  for(unsigned int i = 0; i < multiClusters.size(); i++){
      
    if(multiClusters[i].energy() > MC_best_energy){
      MC_best_index = i;
      MC_best_energy = multiClusters[i].energy();
    }
  } // found MC seed = best
  
  float sum2DClusterEnergy = 0;
  float sum2DClusterAllMEnergy = 0.;
  float sum2DClusterAllSelMEnergy = 0.;
  float sumRecHitAllEnergy = 0.;
  float sum2DClusterEnergy_vsRadius_wrtMCl[5];
  float energy2D_vsLayer[30];
  float energy2DAllM_vsLayer[30];
  float energy2DAllSelM_vsLayer[30];
  float energyRecHitAll_vsLayer[30];
  for(int ij=0; ij<30; ++ij){
    energy2D_vsLayer[ij] = 0.;
    energy2DAllM_vsLayer[ij] = 0.;
    energy2DAllSelM_vsLayer[ij] = 0.;
    energyRecHitAll_vsLayer[ij] = 0.;
    if(ij > 4) continue;
    sum2DClusterEnergy_vsRadius_wrtMCl[ij] = 0;
  }

  //      std::cout << " >>> MC energy = " << multiClusters[i].energy() << " i  = " << i << std::endl;
  for(unsigned int i = 0; i < multiClusters.size(); i++){


    std::array<double,3> fromMC{ {multiClusters[i].x(),multiClusters[i].y(),multiClusters[i].z()} };
    //    std::cout << "MC posX = " << fromMC[0] << " MC posY = " << fromMC[1] << " MC posZ = " <<  fromMC[2] << std::endl;

    float etaMC = multiClusters[i].eta();
    float phiMC = multiClusters[i].phi();

    h_EnergyMC_vsN2D->Fill(multiClusters[i].size(), multiClusters[i].energy());
    h_EnergyFracMC_vsN2D->Fill(multiClusters[i].size(), multiClusters[i].energy() / multiClusters[MC_best_index].energy());
    h_EfracMCl_wrt_Best->Fill(multiClusters[i].energy() / multiClusters[MC_best_index].energy());
    if(multiClusters[i].energy() >= 0.10 * multiClusters[MC_best_index].energy()) ++nmclusCount;

    h_n2dClPerAllMC->Fill(multiClusters[i].size());

    if(MC_best_index == i){
      h_n2dClPerMC->Fill(multiClusters[i].size());
      num2DClinBestMC = multiClusters[i].size();
      float MC_pT = multiClusters[i].energy()/cosh(multiClusters[i].eta());
	
      int cl2dSeed = 0;
      unsigned int layerSeedRH = 0;
	//	std::cout << " >>> nuovo ciclo " << std::endl;
	//	std::cout << " >>> now 2D cluster " << std::endl;
      std::vector<float> Cluster2DZ;
      std::vector<int> Cluster2DIt;
      Cluster2DZ.clear();
      Cluster2DIt.clear();
      for(reco::HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
	  it!=multiClusters[i].end(); it++){
	
	if((*it)->energy() > (*(multiClusters[i].begin()+cl2dSeed))->energy()){
	  cl2dSeed = it - multiClusters[i].begin();
	  //	    std::cout << " Ra = " << (*it)->energy() << std::endl;
	}	


	//not used loop over rechits		
	//	int rhSeed = 0;
	//	  std::cout << " >>> 1st rechit loop " << std::endl;
	/*
	const std::vector< std::pair<DetId, float> > &hf = (*it)->hitsAndFractions();
	for(unsigned int j = 0; j < hf.size(); j++){
	  //here we loop over detid/fraction pairs
	  const HGCRecHit *hit = hitmap[hf[j].first];
	  float rhX = recHitTools.getPosition(hf[j].first).x();
	  float rhY = recHitTools.getPosition(hf[j].first).y();
	  float rhZ = recHitTools.getPosition(hf[j].first).z();
	  float rhEta = recHitTools.getEta(recHitTools.getPosition(hf[j].first));

	  int cmsswLayer = utilsMet.Ztolayer(rhZ, rhEta);
	  int tbLayer = -1;
	  for(int tr=0;tr<8; ++tr){
	    if (CMSSW_cellLayer[tr] == cmsswLayer){
	      tbLayer = tr;
	      break;
	    }
	  }
	  h2_YvsX_rh2dCl[tbLayer]->Fill(rhX, rhY, hit->energy()*hf[j].second);
	  ZvsYvsX_2Drh->Fill(tbLayer+1, rhX, rhY, hit->energy()*hf[j].second);
	}
	*/
	float eta2DCl = (*it)->eta();
	float phi2DCl = (*it)->phi();

	//	std::cout << "2D cluster Z = " << (*it)->z() << " eta = " << eta2DCl << " layer = " << tbLayer <<  " cmsswLayer = " << cmsswLayer << std::endl;

	sum2DClusterEnergy += (*it)->energy();
	energy2D_vsLayer[utilsMet.Ztolayer((*it)->z(), eta2DCl)] += (*it)->energy();


	// Cluster2DZ.push_back(utilsMet.Ztolayer((*it)->z(), eta2DCl) );
	// Cluster2DIt.push_back(it - multiClusters[i].begin());


	float x2DCl = (*it)->x();
	float y2DCl = (*it)->y();
	float z2DCl = (*it)->z();

	float Radius_2DClMCl = utilsMet.dsGenRecoObj(etaMC, phiMC, z2DCl, x2DCl, y2DCl);
	  
	if(int(Radius_2DClMCl) < 5) sum2DClusterEnergy_vsRadius_wrtMCl[int(Radius_2DClMCl)] += (*it)->energy();

	h_dR_2DCl_multiCl->Fill(Radius_2DClMCl);

	std::array<double,3> to{ {0., 0., z2DCl} };
        utilsMet.layerIntersection(to, fromMC, vtx);
	float deltaR = sqrt(pow(to[0]-x2DCl, 2) + pow(to[1]-y2DCl, 2));
	//	std::cout << "best  to[0] = " << to[0] << " x2DCl = " << x2DCl << " to[1] = " << to[1] << " y2DCl = " << y2DCl << std::endl;      
	h_dR_2DCl->Fill(deltaR);
	
      } // loop over 2D cluster
      
	// FIXME LAYER dai recHit
	//    math::XYZPoint mclPosition = clusterTools->getMultiClusterPosition(multiClusters[i]);
	// std::cout << " mclPosition.eta = " << mclPosition.eta() << " mclPosition.phi = " << mclPosition.phi() 
	// 	      << " mcl.eta() = " << multiClusters[i].eta() << " phi = " << multiClusters[i].phi() << std::endl;
      


      h3_seed_yz->Fill((*(multiClusters[i].begin() + cl2dSeed))->z(), (*(multiClusters[i].begin() + cl2dSeed))->y());
      h3_seed_xz->Fill((*(multiClusters[i].begin() + cl2dSeed))->z(), (*(multiClusters[i].begin() + cl2dSeed))->x());

      float eta2dS = (*(multiClusters[i].begin() + cl2dSeed))->eta();
      float phi2dS = (*(multiClusters[i].begin() + cl2dSeed))->phi();

      /*
      for(unsigned int is=0; is<Cluster2DZ.size(); ++is){
	if(Cluster2DIt.at(is) == cl2dSeed){ 
	  //	    std::cout << " delta layer seed = " << Cluster2DZ.at(is) - utilsMet.Ztolayer((*(multiClusters[i].begin() + cl2dSeed))->z(), eta2dS) << std::endl; 
	  continue;
	}
	//	  std::cout << " delta layer altro = " << Cluster2DZ.at(is) - utilsMet.Ztolayer((*(multiClusters[i].begin() + cl2dSeed))->z(), eta2dS) << std::endl; 
	h_deltaLayer_2dClPerMC->Fill(Cluster2DZ.at(is) - utilsMet.Ztolayer((*(multiClusters[i].begin() + cl2dSeed))->z(), eta2dS));
	h2_deltaLayer_n2dClPerMC_vs_2dClPerMC->Fill(Cluster2DZ.at(is) - utilsMet.Ztolayer((*(multiClusters[i].begin() + cl2dSeed))->z(), eta2dS), num2DClinBestMC);
	h2_deltaLayer_n2dClPerMC_vs_2dClEnergy->Fill(Cluster2DZ.at(is) - utilsMet.Ztolayer((*(multiClusters[i].begin() + cl2dSeed))->z(), eta2dS), 
						     (*(multiClusters[i].begin() + Cluster2DIt.at(is)))->energy() / (*(multiClusters[i].begin() + cl2dSeed))->energy() );
      }
      */

      //      h2_dR_SeedCl2_multiCl_vsEta->Fill(std::abs(etaMC), reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
      //      tp_dR_SeedCl2_multiCl_vsEta->Fill(std::abs(etaMC), reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
      //	h2_dR_SeedCl2_multiCl_vsLayerS2d->Fill(layerSeedRH, reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
      //      h_dR_SeedCl2_multiCl->Fill(reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));




      //      h_2DClSumE_overMCE->Fill(sum2DClusterEnergy/multiClusters[i].energy());
      h_2DClSumE->Fill(sum2DClusterEnergy);
	

      for(int kk=0; kk<29; ++kk){
	tp_2DClSumE_vsLayer->Fill(kk+1, energy2D_vsLayer[kk+1]);
	  if(kk > 4) continue;
	  tp_2DClSumE_vsRadius_wrtMCl->Fill(kk, sum2DClusterEnergy_vsRadius_wrtMCl[kk]);
      }
      //	std::cout << " reco::deltaR = " << reco::deltaR(etaMC, phiMC, eta2dS, phi2dS) << " ok = " << sqrt(pow(etaMC-eta2dS, 2 ) + pow(reco::deltaPhi(phiMC,phi2dS), 2) )<< std::endl;


      ///////////////////////////loop over all recHits and get the ones in 2cm around best MC axis
      for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end(); ++it_hit) {
	const HGCalDetId hitid = it_hit->detid();
	bool found = false;
	float rhEnergy = 0;
	std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmap.find(hitid);
	if(trovatore == hitmap.end()){
	  continue;
	}
	else{
	  const HGCRecHit *hit = hitmap[hitid];
	  rhEnergy = hit->energy();
	  found = true;
	}
	if(found){
	  float rhX = recHitTools.getPosition(hitid).x();
	  float rhY = recHitTools.getPosition(hitid).y();
	  int rhL = recHitTools.getLayerWithOffset(hitid);
	  float rhEta = recHitTools.getEta(recHitTools.getPosition(hitid));
	  float rhZ = utilsMet.layerToZ(rhL, rhEta);
	  
	  std::array<double,3> to{ {0., 0., rhZ} };
	  utilsMet.layerIntersection(to, fromMC, vtx);
	  float deltaR = sqrt(pow(to[0] - rhX, 2) + pow(to[1] - rhY, 2));
	  h_dR_recHitAll->Fill(deltaR);
	  if(deltaR < 2){
	    sumRecHitAllEnergy += rhEnergy;
	    energyRecHitAll_vsLayer[rhL] += rhEnergy;
	  }
	}
      }
      ///////////////////////////////////////////////////////////////////////////////////////
      h_recHitAllSumE->Fill(sumRecHitAllEnergy);
      
      for(int ij=0; ij<30; ++ij){
	tp_recHitAllSumE_vsLayer->Fill(ij, energyRecHitAll_vsLayer[ij]);
      }
      
           
    } // MC best
    
    ///////////////////////////////////////////////////////////////
    for(reco::HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
	it!=multiClusters[i].end(); it++){
      num2DClinMC += multiClusters[i].size();	

      float eta2DCl = (*it)->eta();
      float phi2DCl = (*it)->phi();


      int cmsswLayer = utilsMet.Ztolayer((*it)->z(), eta2DCl);
      int tbLayer = -1;
      for(int tr=0;tr<8; ++tr){
	if (CMSSW_cellLayer[tr] == cmsswLayer){
	  tbLayer= tr;
	  break;
	}
      }
      h2_YvsX_2dCl[tbLayer]->Fill((*it)->x(), (*it)->y(), (*it)->energy());	
      ZvsYvsX_2Dcl->Fill(tbLayer+1, (*it)->x(), (*it)->y(), (*it)->energy());
      
      
      sum2DClusterAllMEnergy += (*it)->energy();
      energy2DAllM_vsLayer[utilsMet.Ztolayer((*it)->z(), eta2DCl)] += (*it)->energy();
      num2DClinMC_vsL[utilsMet.Ztolayer((*it)->z(), eta2DCl)] += multiClusters[i].size();
      

      const std::vector< std::pair<DetId, float> > &hf = (*it)->hitsAndFractions();
      for(unsigned int j = 0; j < hf.size(); j++){
	//here we loop over detid/fraction pairs                                                                                                                                
	const HGCRecHit *hit = hitmap[hf[j].first];
	float rhX = recHitTools.getPosition(hf[j].first).x();
	float rhY = recHitTools.getPosition(hf[j].first).y();
	float rhZ = recHitTools.getPosition(hf[j].first).z();
	float rhEta = recHitTools.getEta(recHitTools.getPosition(hf[j].first));

	int cmsswLayer = utilsMet.Ztolayer(rhZ, rhEta);
	int tbLayer = -1;
	for(int tr=0;tr<8; ++tr){
	  if (CMSSW_cellLayer[tr] == cmsswLayer){
	    tbLayer = tr;
	    break;
	  }
	}
	h2_YvsX_rh2dCl[tbLayer]->Fill(rhX, rhY, hit->energy()*hf[j].second);
	ZvsYvsX_2Drh->Fill(tbLayer+1, rhX, rhY, hit->energy()*hf[j].second);
      }

      std::array<double,3> to{ {0., 0., (*it)->z()} };
      utilsMet.layerIntersection(to, fromMC, vtx);
      float deltaR = sqrt(pow(to[0] - (*it)->x(), 2) + pow(to[1] - (*it)->y(), 2));
      //      std::cout << " 2D in M to[0] = " << to[0] << " x2DCl = " << (*it)->x() << " to[1] = " << to[1] << " y2DCl = " << (*it)->y() << std::endl;      
      h_dR_2DClAllM->Fill(deltaR);


      if(MC_best_index == i || multiClusters[i].energy() >= 0.10 * multiClusters[MC_best_index].energy()){
	sum2DClusterAllSelMEnergy += (*it)->energy();
	energy2DAllSelM_vsLayer[utilsMet.Ztolayer((*it)->z(), eta2DCl)] += (*it)->energy();
	h_dR_2DClAllSelM->Fill(deltaR);
      }
    } // loop over 2D cluster
  
    h_2DClAllMSumE->Fill(sum2DClusterAllMEnergy);
    h_2DClAllSelMSumE->Fill(sum2DClusterAllSelMEnergy);
	
    for(int kk=0; kk<30; ++kk){
      tp_2DClAllMSumE_vsLayer->Fill(kk, energy2DAllM_vsLayer[kk]);
      tp_2DClAllSelMSumE_vsLayer->Fill(kk, energy2DAllSelM_vsLayer[kk]);
    }
    /////////////////////////////////////////////////////////////////


    //    std::cout << " multicluster i = " << i << " energy = " << multiClusters[i].energy() << " num 2D = " << multiClusters[i].size() << " layer = "
    //	      << utilsMet.Ztolayer(multiClusters[i].z(), multiClusters[i].eta()) << std::endl;
  } // loop MC

  h_nMCvsnum2dClBestMC->Fill(num2DClinBestMC, nmclusCount);  
  h_nMCPerEvt->Fill(nmclusCount);
  h_n2dClinMCPerEvt->Fill(num2DClinMC);
  for(int ij=0; ij<30; ++ij)  n2DClinM_vsL->Fill(ij, num2DClinMC_vsL[ij]);


  //shower axis by recHits
  float axisX = 0;
  float axisY = 0;
  float axisZ = 0;
  float sumEnergyToNorm = 0;
  GlobalPoint showerAxis;
  
  float sumMtcRHEnergy = 0;
  float sumMtcRHEnergy_vsLayer[30];
  for(int ij=0; ij<30; ++ij){
    sumMtcRHEnergy_vsLayer[ij] = 0;
  }
  
  //  const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
  // loop over EE RecHits                                                                                                                                                   
  for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end(); ++it_hit) {
    const HGCalDetId hitid = it_hit->detid();
    
    std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmap.find(hitid);
    if(trovatore == hitmap.end()){
      //	  std::cout << " fine " << std::endl;
      continue;
    }
    else{
      const HGCRecHit *hit = hitmap[hitid];
      float rhEnergy = hit->energy();
      sumEnergyToNorm += rhEnergy;
      axisX += rhEnergy * recHitTools.getPosition(hitid).x();
      axisY += rhEnergy * recHitTools.getPosition(hitid).y();
      axisZ += rhEnergy * recHitTools.getPosition(hitid).z();
      }
  }
  showerAxis = GlobalPoint(axisX/sumEnergyToNorm, axisY/sumEnergyToNorm, (axisZ - 0)/sumEnergyToNorm);
  
  float axEta = showerAxis.eta();
  float axPhi = showerAxis.phi();
  float axX = showerAxis.x();
  float axY = showerAxis.y();
  float axZ = showerAxis.z();
      
  
  
  // std::cout << " axEta = " << axEta << std::endl;
  // std::cout << " axPhi = " << axPhi << std::endl;
  
  //      unsigned int allRHSeedLayer = 0; 
  float allRHSeedEnergy = 0;   
  //loop over hits
  for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end(); ++it_hit) {
    const HGCalDetId hitid = it_hit->detid();
    
    bool found = false;
    float rhEnergy = 0;
    //	float rhTime = -1;
    
    std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmap.find(hitid);
    if(trovatore == hitmap.end()){
      continue;
    }
    else{
      const HGCRecHit *hit = hitmap[hitid];
      rhEnergy = hit->energy();
      //	  rhTime = hit->time() - 1;
      found = true;
      /*
	if(it_caloPart->eta() > 0){
	h_dEta_pos->Fill(recHitTools.getEta(hitid) - it_caloPart->eta());
	h_dPhi_pos->Fill(reco::deltaPhi(recHitTools.getPhi(hitid), it_caloPart->phi()) );
	}
	if(it_caloPart->eta() < 0){
	h_dEta_neg->Fill(recHitTools.getEta(hitid) - it_caloPart->eta());
	h_dPhi_neg->Fill(reco::deltaPhi(recHitTools.getPhi(hitid), it_caloPart->phi()) );
	}
      */
      
      //here look for seed
      if(hit->energy()> allRHSeedEnergy){
	allRHSeedEnergy = hit->energy();
	//	    allRHSeedLayer = recHitTools.getLayerWithOffset(hitid);
      }
    }
    if(found){
      float rhEta = recHitTools.getEta(hitid);
      float rhPhi = recHitTools.getPhi(hitid);
      
      sumMtcRHEnergy += rhEnergy;      
      
      h2_dPhivsdEta_rhAxis->Fill(reco::deltaPhi(rhPhi, axPhi), rhEta - axEta);
      
      float rhX = recHitTools.getPosition(hitid).x();
      float rhY = recHitTools.getPosition(hitid).y();
      float rhL = recHitTools.getLayerWithOffset(hitid);
      float rhZ = utilsMet.layerToZ(rhL, rhEta);  
      sumMtcRHEnergy_vsLayer[ utilsMet.Ztolayer(rhZ, rhEta)] += rhEnergy;


      int cmsswLayer = utilsMet.Ztolayer(rhZ, rhEta) - 1;
      int tbLayer = -1;
      for(int tr=0;tr<8; ++tr){
	if (CMSSW_cellLayer[tr] == cmsswLayer){
	  tbLayer= tr;
	  break;
	}
      }
      if( rhEnergy > 0){
	h2_YvsX_rh[tbLayer]->Fill(rhX, rhY, rhEnergy);
	ZvsYvsX_rh->Fill(tbLayer+1, rhX, rhY, rhEnergy);
      }

      //      std::cout << "rh all  Z = " << rhZ << " eta = " << rhEta << " layer = " << tbLayer <<  " cmsswLayer = " << cmsswLayer << " energy = " << rhEnergy << " detID = " << hitid.rawId() << std::endl;
    }
  }// first loop over rechits
  
  h_mtcRHSumE -> Fill(sumMtcRHEnergy);
  for(int ij=0; ij<30; ++ij){
    tp_mtcRHSumE_vsLayer->Fill(ij, sumMtcRHEnergy_vsLayer[ij]);
  }


  h_n2dClPerEvt->Fill(clusters.size());
  float sum2DClusterAllEnergy = 0.;
  float energy2DAll_vsLayer[30];
  int num2DCl_vsL[30];
  for(int ij=0; ij<30; ++ij) {
    energy2DAll_vsLayer[ij] = 0.;
    num2DCl_vsL[ij] = 0;
  }

  //loop over 2D cluster
  int cl2dSeed = 0; 
  for (reco::CaloClusterCollection::const_iterator it_cl = clusters.begin(); it_cl < clusters.end(); ++it_cl) {
    float eta2DCl = it_cl->eta();
    float phi2DCl = it_cl->phi();
    
    if(it_cl->energy() > (clusters.begin()+cl2dSeed)->energy()){
      cl2dSeed = it_cl - clusters.begin();
      //        std::cout << " Ra = " << (*it)->energy() << std::endl;                                                                                                                                  
    }
  }

  if(clusters.size() == 0) return;

  std::array<double,3> fromCl{ {(clusters.begin() + cl2dSeed)->x(), (clusters.begin() + cl2dSeed)->y(), (clusters.begin() + cl2dSeed)->z()} };
  //  std::cout << "Cl posX = " << fromCl[0] << " Cl posY = " << fromCl[1] << " Cl posZ = " <<  fromCl[2] << std::endl;
    
  for (reco::CaloClusterCollection::const_iterator it_cl = clusters.begin(); it_cl < clusters.end(); ++it_cl) {
    float eta2DCl = it_cl->eta();
    float phi2DCl = it_cl->phi();
    

    h3_2dCl_yz->Fill(it_cl->z(), it_cl->y());
    h3_2dCl_xz->Fill(it_cl->z(), it_cl->x());

    sum2DClusterAllEnergy += it_cl->energy();
    energy2DAll_vsLayer[utilsMet.Ztolayer(it_cl->z(), eta2DCl)] += it_cl->energy();      

    num2DCl_vsL[utilsMet.Ztolayer(it_cl->z(), eta2DCl)] += 1;


    std::array<double,3> to{ {0., 0., it_cl->z()} };
    utilsMet.layerIntersection(to, fromCl, vtx);
    float deltaR = sqrt(pow(to[0] - it_cl->x(), 2) + pow(to[1] - it_cl->y(), 2));
    //    std::cout << "all 2Dcl to[0] = " << to[0] << " x2DCl = " << it_cl->x() << " to[1] = " << to[1] << " y2DCl = " << it_cl->y() << std::endl;      
    if(cl2dSeed != it_cl - clusters.begin())    h_dR_2DClAll->Fill(deltaR);

  } // loop over 2D cluster
    
  h_2DClAllSumE->Fill(sum2DClusterAllEnergy);
	
  for(int kk=0; kk<30; ++kk){
    tp_2DClAllSumE_vsLayer->Fill(kk, energy2DAll_vsLayer[kk]);
    n2DCl_vsL->Fill(kk, num2DCl_vsL[kk]);
  }

  
    
}

void
HGCalTBAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalTBAnalysis::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalTBAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalTBAnalysis);
