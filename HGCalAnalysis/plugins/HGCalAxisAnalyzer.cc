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

//#include "RecoNtuples/HGCalAxisAnalyzer/interface/AEvent.h"
//#include "RecoNtuples/HGCalAxisAnalyzer/interface/AObData.h"

#include "RecoNtuples/HGCalAnalysis/interface/UsefulClasses.h"
//#include "../utils/UsefulClasses.cpp"

#include <string>
#include <map>

class HGCalAxisAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalAxisAnalyzer(const edm::ParameterSet&);
  ~HGCalAxisAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsFH;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsBH;
  edm::EDGetTokenT<std::vector<TrackingVertex> > _vtx;
  edm::EDGetTokenT<std::vector<TrackingParticle> > _part;
  //  edm::EDGetTokenT<std::vector<SimCluster> > _simClusters;
  edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;
  edm::EDGetTokenT<reco::CaloClusterCollection> _clusters;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > _multiClusters;

  std::string                detector;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  bool                       rawRecHits;
  hgcal::RecHitTools         recHitTools;

  std::unique_ptr<hgcal::ClusterTools> clusterTools;
  //  UsefulClasses utilsMet;

  TH1F* h_Vtx_x;
  TH1F* h_Vtx_y;
  TH1F* h_Vtx_z;

  TH1F* h_Vtx_dvx;
  TH1F* h_Vtx_dvy;
  TH1F* h_Vtx_dvz;

  TH2F* h3_bary_yz;
  TH2F* h3_bary_xz;
  TH2F* h3_seed_yz;
  TH2F* h3_seed_xz;
  TH2F* h3_2dCl_yz;
  TH2F* h3_2dCl_xz;
  //  TH3* h3_2dCl_xyz;

  std::vector<int> nHits_layer;
  std::vector<int> nHitsWithTime_layer;

  //  TH1F* h_numberOfMC;
  TH1F* h_energyFractionInMC;
  TH1F* h_energyFractionInBestMC;
  TH1F* h_energyFractionInSelectedMC;
  TH1F* h_EfracMCl_wrt_Best;
  TH1F* h_EfracSelectedMCl_wrt_Best;

  TH2F* h2_dPhiMCtoGen_vsPt;
  TH1F* h_dPhiMCtoGen;

  TH2F* h2_dR_SeedCl2_multiCl_vsEta;
  TProfile* tp_dR_SeedCl2_multiCl_vsEta;
  //  TH2F* h2_dR_SeedCl2_multiCl_vsLayerS2d;
  TH1F* h_dR_SeedCl2_multiCl;

  TH2F* h_nMCvsnum2dClBestMC;
  TH1F* h_n2dClPerMC;
  TH1F* h_deltaLayer_2dClPerMC;
  TH2F* h2_deltaLayer_n2dClPerMC_vs_2dClPerMC;
  TH2F* h2_deltaLayer_n2dClPerMC_vs_2dClEnergy;

  TH1F* h_nMCPerEvt;

  TH2F* h2_dPhi2DCltoGen_vsPt;
  TH1F* h_dPhi2DCltoGen;
  TH1F* h_dR_2DCl_multiCl;
  TH1F* h_dR_2DCl_Gen;

  TH1F* h_2DClSumE_overMCE;
  TH1F* h_2DClSumE;
  TH2F* h2_2DClSumE_vsVtxZ;
  TProfile* tp_2DClSumE_vsVtxZ;

  TProfile* tp_2DClSumE_vsDR_wrtGen;
  TProfile* tp_2DClSumE_vsDR_wrtMCl;
  TProfile* tp_2DClSumE_vsRadius_wrtGen;
  TProfile* tp_2DClSumE_vsRadius_wrtMCl;

  //selected
  TH2F* h2_dPhiSelMCtoGen_vsPt;
  TH1F* h_dPhiSelMCtoGen;

  TH2F* h2_dR_SeedCl2_SelmultiCl_vsEta;
  TProfile* tp_dR_SeedCl2_SelmultiCl_vsEta;
  //  TH2F* h2_dR_SeedCl2_SelmultiCl_vsLayerS2d;
  TH1F* h_dR_SeedCl2_SelmultiCl;

  TH2F* h_nSelMCvsnum2dClSelMC;
  TH1F* h_n2dClPerSelMC;
  TH1F* h_nSelMCPerEvt;

  TH2F* h2_dPhi2DClSeltoGen_vsPt;
  TH1F* h_dPhi2DClSeltoGen;
  TH1F* h_dR_2DClSel_multiCl;
  TH1F* h_dR_2DClSel_Gen;

  //  TH1F* h_2DClSelSumE_overMCE;
  TH1F* h_2DClSelSumE;
  TProfile* tp_2DClSelSumE_vsDR_wrtGen;
  TProfile* tp_2DClSelSumE_vsDR_wrtMCl;
  TProfile* tp_2DClSelSumE_vsRadius_wrtGen;
  TProfile* tp_2DClSelSumE_vsRadius_wrtMCl;
  //
  TH1F* h_2DClSelConsSumE;
  TProfile* tp_2DClSelConsSumE_vsRadius_wrtGen;
  TProfile* tp_2DClSelConsSumE_vsRadius_wrtMCl;

  TH2F* h2_dPhiMatchHittoGen_vsPt;
  TH1F* h_dPhiMatchHittoGen;

  TH1F* h_dR_MatchHit_Gen;

  TH1F* h_mtcRHSumE;
  TProfile* tp_mtcRHSumE_vsDR_wrtGen;
  TProfile* tp_mtcRHSumE_vsRadius_wrtGen;

  TH2F* h2_dPhivsdEta_rhGen;
  TH2F* h2_dPhivsdEta_rhGen_fixL;
  TH2F* h2_dPhivsdEta_rhAxis;
  TH2F* h2_dPhivsdEta_GenAxis;


  /*  
  TH1F* h_dEta_pos;
  TH1F* h_dPhi_pos;
  TH1F* h_dEta_neg;
  TH1F* h_dPhi_neg;
  */
};  



HGCalAxisAnalyzer::HGCalAxisAnalyzer(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  rawRecHits(iConfig.getParameter<bool>("rawRecHits"))
{

  //  std::cout << " >>> costruttore " << std::endl;

  //now do what ever initialization is needed
  usesResource("TFileService");

  if(detector=="all") {
    _recHitsEE = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"));
    _recHitsFH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCFHInput"));
    _recHitsBH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCBHInput"));
    algo = 1;
  }else if(detector=="EM") {
    _recHitsEE = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"));
    algo = 2;
  }else if(detector=="HAD") {
    _recHitsFH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCFHInput"));
    _recHitsBH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCBHInput"));
    algo = 3;
  }
  _vtx = consumes<std::vector<TrackingVertex> >(edm::InputTag("mix","MergedTrackTruth"));
  _part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));
  //  _simClusters = consumes<std::vector<SimCluster> >(edm::InputTag("mix","MergedCaloTruth"));
  _caloParticles = consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));
  _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("hgcalLayerClusters"));
  _multiClusters = consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("hgcalLayerClusters"));


  auto sumes = consumesCollector();
  clusterTools = std::make_unique<hgcal::ClusterTools>(iConfig,sumes);

  //  utilsMet = UsefulClasses();

  edm::Service<TFileService> fs;

  h_Vtx_x = fs->make<TH1F>("h_Vtx_x", "", 1000, -10., 10.);
  h_Vtx_y = fs->make<TH1F>("h_Vtx_y", "", 1000, -10., 10.);
  h_Vtx_z = fs->make<TH1F>("h_Vtx_z", "", 1000, -10., 10.);

  h_Vtx_dvx = fs->make<TH1F>("h_Vtx_dvx", "", 1000, -10., 10.);
  h_Vtx_dvy = fs->make<TH1F>("h_Vtx_dvy", "", 1000, -10., 10.);
  h_Vtx_dvz = fs->make<TH1F>("h_Vtx_dvz", "", 1000, -10., 10.);

  h3_bary_xz = fs->make<TH2F>("h3_bary_xz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_bary_yz = fs->make<TH2F>("h3_bary_yz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_seed_xz = fs->make<TH2F>("h3_seed_xz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_seed_yz = fs->make<TH2F>("h3_seed_yz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_2dCl_xz = fs->make<TH2F>("h3_2dCl_xz", "",  1000, 250., 650., 1000, -100., 100.);
  h3_2dCl_yz = fs->make<TH2F>("h3_2dCl_yz", "",  1000, 250., 650., 1000, -100., 100.);
  //  h3_2dCl_xyz = fs->make<TH3>("h3_2dCl_xyz", "",  1000, 200., 500., 1000, -100., 100., 1000, -100., 100.);

  h_energyFractionInMC = fs->make<TH1F>("h_energyFractionInMC", "", 1000, 0., 2.);
  h_energyFractionInBestMC = fs->make<TH1F>("h_energyFractionInBestMC", "", 1000, 0., 2.);
  h_energyFractionInSelectedMC = fs->make<TH1F>("h_energyFractionInSelectedMC", "", 1000, 0., 2.);
  h_EfracMCl_wrt_Best = fs->make<TH1F>("h_EfracMCl_wrt_Best", "", 500, 0., 1.1);
  h_EfracSelectedMCl_wrt_Best = fs->make<TH1F>("h_EfracSelectedMCl_wrt_Best", "", 500, 0., 1.1);

  h2_dPhiMCtoGen_vsPt = fs->make<TH2F>("h2_dPhiMCtoGen_vsPt", "", 300, 0., 300., 600, -0.3, 0.3);
  h_dPhiMCtoGen = fs->make<TH1F>("h_dPhiMCtoGen", "", 600, -3., 3.);

  h_dR_SeedCl2_multiCl = fs->make<TH1F>("h_dR_SeedCl2_multiCl", "", 1000, 0., 3.);
  //  h2_dR_SeedCl2_multiCl_vsLayerS2d = fs->make<TH2F>("h2_dR_SeedCl2_multiCl_vsLayerS2d", "", 100, 0., 100., 1000, 0., 5.);
  h2_dR_SeedCl2_multiCl_vsEta = fs->make<TH2F>("h2_dR_SeedCl2_multiCl_vsEta", "", 500, 1.5, 3.5, 1000, 0., 5.);
  tp_dR_SeedCl2_multiCl_vsEta = fs->make<TProfile>("tp_dR_SeedCl2_multiCl_vsEta", "", 500, 1.5, 3.5);

  h_nMCvsnum2dClBestMC = fs->make<TH2F>("h_nMCvsnum2dClBestMC", "", 50, 0., 50., 100, 0., 100.);
  h_n2dClPerMC = fs->make<TH1F>("h_n2dClPerMC", "", 50, 0., 50.);
  h_nMCPerEvt = fs->make<TH1F>("h_nMCPerEvt", "", 100, 0., 100.);
  h_deltaLayer_2dClPerMC = fs->make<TH1F>("h_deltaLayer_2dClPerMC", "", 50, -25., 25.);
  h2_deltaLayer_n2dClPerMC_vs_2dClPerMC = fs->make<TH2F>("h2_deltaLayer_n2dClPerMC_vs_2dClPerMC", "", 50, -25., 25., 50, 0., 50.);
  h2_deltaLayer_n2dClPerMC_vs_2dClEnergy = fs->make<TH2F>("h2_deltaLayer_n2dClPerMC_vs_2dClEnergy", "", 50, -25., 25., 100, 0., 1.);


  h2_dPhi2DCltoGen_vsPt = fs->make<TH2F>("h2_dPhi2DCltoGen_vsPt", "", 300, 0., 300., 600, -0.3, 0.3);
  h_dPhi2DCltoGen = fs->make<TH1F>("h_dPhi2DCltoGen", "", 600, -3., 3.);

  h_dR_2DCl_multiCl = fs->make<TH1F>("h_dR_2DCl_multiCl", "", 1000, 0., 3.);
  h_dR_2DCl_Gen = fs->make<TH1F>("h_dR_2DCl_Gen", "", 1000, 0., 3.);

  h_2DClSumE_overMCE = fs->make<TH1F>("h_2DClSumE_overMCE", "", 500, 0., 2.);  
  h_2DClSumE = fs->make<TH1F>("h_2DClSumE", "", 500, 0., 2.);
  h2_2DClSumE_vsVtxZ = fs->make<TH2F>("h2_2DClSumE_vsVtxZ", "", 1000, -10., 10, 500, 0., 2.);
  tp_2DClSumE_vsVtxZ = fs->make<TProfile>("tp_2DClSumE_vsVtxZ", "", 1000, -10., 10);


  tp_2DClSumE_vsDR_wrtGen = fs->make<TProfile>("tp_2DClSumE_vsDR_wrtGen", "", 5, 0., 5);
  tp_2DClSumE_vsDR_wrtMCl = fs->make<TProfile>("tp_2DClSumE_vsDR_wrtMCl", "", 5, 0., 5);
  tp_2DClSumE_vsRadius_wrtGen = fs->make<TProfile>("tp_2DClSumE_vsRadius_wrtGen", "", 5, 0., 5);
  tp_2DClSumE_vsRadius_wrtMCl = fs->make<TProfile>("tp_2DClSumE_vsRadius_wrtMCl", "", 5, 0., 5);


  //sel
  h2_dPhiSelMCtoGen_vsPt = fs->make<TH2F>("h2_dPhiSelMCtoGen_vsPt", "", 300, 0., 300., 600, -0.3, 0.3);
  h_dPhiSelMCtoGen = fs->make<TH1F>("h_dPhiSelMCtoGen", "", 600, -3., 3.);

  h_dR_SeedCl2_SelmultiCl = fs->make<TH1F>("h_dR_SeedCl2_SelmultiCl", "", 1000, 0., 3.);
  //  h2_dR_SeedCl2_SelmultiCl_vsLayerS2d = fs->make<TH2F>("h2_dR_SeedCl2_SelmultiCl_vsLayerS2d", "", 100, 0., 100., 1000, 0., 5.);
  h2_dR_SeedCl2_SelmultiCl_vsEta = fs->make<TH2F>("h2_dR_SeedCl2_SelmultiCl_vsEta", "", 500, 1.5, 3.5, 1000, 0., 5.);
  tp_dR_SeedCl2_SelmultiCl_vsEta = fs->make<TProfile>("tp_dR_SeedCl2_SelmultiCl_vsEta", "", 500, 1.5, 3.5);

  h_nSelMCvsnum2dClSelMC = fs->make<TH2F>("h_nSelMCvsnum2dClSelMC", "", 50, 0., 50., 100, 0., 100.);
  h_n2dClPerSelMC = fs->make<TH1F>("h_n2dClPerSelMC", "", 50, 0., 50.);
  h_nSelMCPerEvt = fs->make<TH1F>("h_nSelMCPerEvt", "", 100, 0., 100.);

  h2_dPhi2DClSeltoGen_vsPt = fs->make<TH2F>("h2_dPhi2DClSeltoGen_vsPt", "", 300, 0., 300., 600, -0.3, 0.3);
  h_dPhi2DClSeltoGen = fs->make<TH1F>("h_dPhi2DClSeltoGen", "", 600, -3., 3.);

  h_dR_2DClSel_multiCl = fs->make<TH1F>("h_dR_2DClSel_multiCl", "", 1000, 0., 3.);
  h_dR_2DClSel_Gen = fs->make<TH1F>("h_dR_2DClSel_Gen", "", 1000, 0., 3.);

  //  h_2DClSelSumE_overMCE = fs->make<TH1F>("h_2DClSelSumE_overMCE", "", 500, 0., 2.);  
  h_2DClSelSumE = fs->make<TH1F>("h_2DClSelSumE", "", 500, 0., 2.);
  tp_2DClSelSumE_vsDR_wrtGen = fs->make<TProfile>("tp_2DClSelSumE_vsDR_wrtGen", "", 5, 0., 5);
  tp_2DClSelSumE_vsDR_wrtMCl = fs->make<TProfile>("tp_2DClSelSumE_vsDR_wrtMCl", "", 5, 0., 5);
  tp_2DClSelSumE_vsRadius_wrtGen = fs->make<TProfile>("tp_2DClSelSumE_vsRadius_wrtGen", "", 5, 0., 5);
  tp_2DClSelSumE_vsRadius_wrtMCl = fs->make<TProfile>("tp_2DClSelSumE_vsRadius_wrtMCl", "", 5, 0., 5);
  //

  h_2DClSelConsSumE = fs->make<TH1F>("h_2DClSelConsSumE", "", 500, 0., 2.);
  tp_2DClSelConsSumE_vsRadius_wrtGen = fs->make<TProfile>("tp_2DClSelConsSumE_vsRadius_wrtGen", "", 5, 0., 5);
  tp_2DClSelConsSumE_vsRadius_wrtMCl = fs->make<TProfile>("tp_2DClSelConsSumE_vsRadius_wrtMCl", "", 5, 0., 5);


  h2_dPhiMatchHittoGen_vsPt = fs->make<TH2F>("h2_dPhiMatchHittoGen_vsPt", "", 300, 0., 300., 600, -0.3, 0.3);
  h_dPhiMatchHittoGen = fs->make<TH1F>("h_dPhiMatchHittoGen", "", 600, -3., 3.);
  h_dR_MatchHit_Gen = fs->make<TH1F>("h_dR_MatchHit_Gen", "", 1000, 0., 3.);


  h_mtcRHSumE = fs->make<TH1F>("h_mtcRHSumE", "", 500, 0., 2.);
  tp_mtcRHSumE_vsDR_wrtGen = fs->make<TProfile>("tp_mtcRHSumE_vsDR_wrtGen", "", 5, 0., 5);
  tp_mtcRHSumE_vsRadius_wrtGen = fs->make<TProfile>("tp_mtcRHSumE_vsRadius_wrtGen", "", 5, 0., 5);


  h2_dPhivsdEta_rhGen = fs->make<TH2F>("h2_dPhivsdEta_rhGen", "", 500, -0.2, 0.2, 500, -0.2, 0.2);
  h2_dPhivsdEta_rhGen_fixL = fs->make<TH2F>("h2_dPhivsdEta_rhGen_fixL", "", 500, -0.2, 0.2, 500, -0.2, 0.2);
  h2_dPhivsdEta_rhAxis = fs->make<TH2F>("h2_dPhivsdEta_rhAxis", "", 1000, -2., 2., 1000, -2., 2.);
  h2_dPhivsdEta_GenAxis = fs->make<TH2F>("h2_dPhivsdEta_GenAxis", "", 1000, -2., 2., 1000, -2., 2.);

  /*
    h_dEta_pos = fs->make<TH1F>("h_dEta_pos", "", 500, -0.1, 0.1);
    h_dPhi_pos = fs->make<TH1F>("h_dPhi_pos", "", 500, -0.1, 0.1);
    h_dEta_neg = fs->make<TH1F>("h_dEta_neg", "", 500, -0.1, 0.1);
    h_dPhi_neg = fs->make<TH1F>("h_dPhi_neg", "", 500, -0.1, 0.1);
  */

}

HGCalAxisAnalyzer::~HGCalAxisAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
HGCalAxisAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  std::cout << " >>> analyzer " << std::endl;
  using namespace edm;

  recHitTools.getEventSetup(iSetup);

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;

  Handle<std::vector<TrackingVertex> > vtxHandle;
  Handle<std::vector<TrackingParticle> > partHandle;
  iEvent.getByToken(_vtx,vtxHandle);
  iEvent.getByToken(_part,partHandle);
  const std::vector<TrackingVertex>& vtxs = *vtxHandle;
  const std::vector<TrackingParticle>& part = *partHandle;

  // Handle<std::vector<SimCluster> > simClusterHandle;
  // iEvent.getByToken(_simClusters, simClusterHandle);
  // const std::vector<SimCluster>& simClusters = *simClusterHandle;

  Handle<std::vector<CaloParticle> > caloParticleHandle;
  iEvent.getByToken(_caloParticles, caloParticleHandle);
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;

  Handle<reco::CaloClusterCollection> clusterHandle;
  iEvent.getByToken(_clusters,clusterHandle);
  const reco::CaloClusterCollection &clusters = *clusterHandle;
  /*
  edm::PtrVector<reco::CaloCluster> clusterPtrs;
  for( unsigned i = 0; i < clusterHandle->size(); ++i ){
    edm::Ptr<reco::CaloCluster> ptr(clusterHandle,i);
    clusterPtrs.push_back(ptr);
  }
  */
  Handle<std::vector<reco::HGCalMultiCluster> > multiClusterHandle;
  iEvent.getByToken(_multiClusters, multiClusterHandle);
  const std::vector<reco::HGCalMultiCluster>& multiClusters = *multiClusterHandle;

  float vx = 0.;
  float vy = 0.;
  float vz = 0.;
  if(vtxs.size()!=0){
    vx = vtxs[0].position().x();
    vy = vtxs[0].position().y();
    vz = vtxs[0].position().z();
  }

  h_Vtx_x->Fill(vx);
  h_Vtx_y->Fill(vy);
  h_Vtx_z->Fill(vz);

  // TODO: should fall back to beam spot if no vertex
  //  npart = part.size();
  for(unsigned int i=0;i<part.size();++i){
    if(part[i].parentVertex()->nGenVertices()>0){
      float dvx=0.;
      float dvy=0.;
      float dvz=0.;
      if(part[i].decayVertices().size()==1){
	 dvx=part[i].decayVertices()[0]->position().x();
	 dvy=part[i].decayVertices()[0]->position().y();
	 dvz=part[i].decayVertices()[0]->position().z();

	 h_Vtx_dvx->Fill(dvx);
	 h_Vtx_dvy->Fill(dvy);
	 h_Vtx_dvz->Fill(dvz);
      }
      //      agpc->push_back(AGenPart(part[i].eta(),part[i].phi(),part[i].pt(),part[i].energy(),dvx,dvy,dvz,part[i].pdgId()));
    }
  }

  //make a map detid-rechit
  std::map<DetId,const HGCRecHit*> hitmap;
  switch(algo){
  case 1:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);
      const auto& rechitsEE = *recHitHandleEE;
      const auto& rechitsFH = *recHitHandleFH;
      const auto& rechitsBH = *recHitHandleBH;
      for(unsigned int i = 0; i < rechitsEE.size(); ++i){
	hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      for(unsigned int i = 0; i < rechitsFH.size(); ++i){
	hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for(unsigned int i = 0; i < rechitsBH.size(); ++i){
	hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      for(unsigned int i = 0; i < rechitsEE.size(); i++){
	hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      break;
    }
  case 3:
    {
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);
      const auto& rechitsFH = *recHitHandleFH;
      const auto& rechitsBH = *recHitHandleBH;
      for(unsigned int i = 0; i < rechitsFH.size(); i++){
	hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for(unsigned int i = 0; i < rechitsBH.size(); i++){
	hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
  default:
    break;
  }

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

  // loop over caloParticles
  for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){

    float etaGen = it_caloPart->eta();
    float phiGen = it_caloPart->phi();

    UsefulClasses utilsMet = UsefulClasses(etaGen, phiGen); 
      
    int nmclusCount = 0;
    int nSelmclusCount = 0;
    int num2DClinBestMC = 0;
    int num2DClinSelMC = 0;
    unsigned int MC_best_index = 0;
    float MC_best_energy = 0;
    //    std::cout << " >>> now MCL " << std::endl;
    for(unsigned int i = 0; i < multiClusters.size(); i++){
      if(multiClusters[i].eta() * it_caloPart->eta() < 0) continue;

      if(multiClusters[i].energy() > MC_best_energy){
	MC_best_index = i;
	MC_best_energy = multiClusters[i].energy();
      }
    } // found MC seed = best

    float sum2DClusterEnergy = 0;
    float sum2DClusterEnergy_vsDR_wrtGen[5];
    float sum2DClusterEnergy_vsDR_wrtMCl[5];
    float sum2DClusterEnergy_vsRadius_wrtGen[5];
    float sum2DClusterEnergy_vsRadius_wrtMCl[5];
    float sum2DClSelEnergy = 0;
    float sum2DClSelEnergy_vsDR_wrtGen[5];
    float sum2DClSelEnergy_vsDR_wrtMCl[5];
    float sum2DClSelEnergy_vsRadius_wrtGen[5];
    float sum2DClSelEnergy_vsRadius_wrtMCl[5];
    float sum2DClSelConsEnergy = 0;
    float sum2DClSelConsEnergy_vsRadius_wrtGen[5];
    float sum2DClSelConsEnergy_vsRadius_wrtMCl[5];
    for(int ij=0; ij<5; ++ij){
      sum2DClusterEnergy_vsDR_wrtGen[ij] = 0;
      sum2DClusterEnergy_vsRadius_wrtGen[ij] = 0;
      sum2DClusterEnergy_vsDR_wrtMCl[ij] = 0;
      sum2DClusterEnergy_vsRadius_wrtMCl[ij] = 0;
      sum2DClSelEnergy_vsDR_wrtGen[ij] = 0;
      sum2DClSelEnergy_vsRadius_wrtGen[ij] = 0;
      sum2DClSelEnergy_vsDR_wrtMCl[ij] = 0;
      sum2DClSelEnergy_vsRadius_wrtMCl[ij] = 0;
      sum2DClSelConsEnergy_vsRadius_wrtGen[ij] = 0;
      sum2DClSelConsEnergy_vsRadius_wrtMCl[ij] = 0;
    }
    //      std::cout << " >>> MC energy = " << multiClusters[i].energy() << " i  = " << i << std::endl;
    for(unsigned int i = 0; i < multiClusters.size(); i++){
      if(multiClusters[i].eta() * it_caloPart->eta() < 0) continue;

      float etaMC = multiClusters[i].eta();
      float phiMC = multiClusters[i].phi();

      h_EfracMCl_wrt_Best->Fill(multiClusters[i].energy() / multiClusters[MC_best_index].energy());
      h_energyFractionInMC->Fill(multiClusters[i].energy()/it_caloPart->energy()); 
      if(multiClusters[i].energy() >= 0.10 * multiClusters[MC_best_index].energy()) ++nmclusCount;


      if(multiClusters[i].size() >= 3 || MC_best_index == i) {
	nSelmclusCount += 1;
	h_EfracSelectedMCl_wrt_Best->Fill(multiClusters[i].energy() / multiClusters[MC_best_index].energy());
	h_energyFractionInSelectedMC->Fill(multiClusters[i].energy()/it_caloPart->energy()); 

	num2DClinSelMC += multiClusters[i].size(); 

	h2_dPhiSelMCtoGen_vsPt->Fill(it_caloPart->pt(), reco::deltaPhi(multiClusters[i].phi(), it_caloPart->phi()));
        h_dPhiSelMCtoGen->Fill(reco::deltaPhi(multiClusters[i].phi(), it_caloPart->phi()));

	bool areConsecutive = false;
	int layerCL[52];

	//look for MC 2DCseed
	for(int il=0; il<52; ++il) layerCL[il] = 0;
	int cl2dSeed = 0;
        unsigned int layerSeedRH = 0;
        for(reco::HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
            it!=multiClusters[i].end(); it++){
	  if((*it)->energy() > (*(multiClusters[i].begin()+cl2dSeed))->energy()){
            cl2dSeed = it - multiClusters[i].begin();
	  }
  
	  sum2DClSelEnergy += (*it)->energy();

          float eta2DCl = (*it)->eta();
          float phi2DCl = (*it)->phi();

          float dR_2DClSelGen = reco::deltaR(etaGen, eta2DCl, phiGen, phi2DCl);
          float dR_2DClSelMCl = reco::deltaR(multiClusters[MC_best_index].eta(), eta2DCl, multiClusters[MC_best_index].phi(), phi2DCl);

          if(int(dR_2DClSelMCl) < 5) sum2DClSelEnergy_vsDR_wrtMCl[int(dR_2DClSelMCl)] += (*it)->energy();
          if(int(dR_2DClSelGen) < 5) sum2DClSelEnergy_vsDR_wrtGen[int(dR_2DClSelGen)] += (*it)->energy();


          float x2DCl = (*it)->x();
          float y2DCl = (*it)->y();
          float z2DCl = (*it)->z();

	  layerCL[utilsMet.Ztolayer(z2DCl, eta2DCl)] += 1;

          float Radius_2DClSelGen = utilsMet.dsGenRecoObj(etaGen, phiGen, z2DCl, x2DCl, y2DCl);
          float Radius_2DClSelMCl = utilsMet.dsGenRecoObj(multiClusters[MC_best_index].eta(), multiClusters[MC_best_index].phi(), z2DCl, x2DCl, y2DCl);

          if(int(Radius_2DClSelGen) < 5) sum2DClSelEnergy_vsRadius_wrtGen[int(Radius_2DClSelGen)] += (*it)->energy();
          if(int(Radius_2DClSelMCl) < 5) sum2DClSelEnergy_vsRadius_wrtMCl[int(Radius_2DClSelMCl)] += (*it)->energy();


          h2_dPhi2DClSeltoGen_vsPt->Fill(it_caloPart->pt(), reco::deltaPhi(phiGen, phi2DCl));
          h_dPhi2DClSeltoGen->Fill(reco::deltaPhi(phiGen, phi2DCl));

          h_dR_2DClSel_multiCl->Fill(Radius_2DClSelMCl);
          h_dR_2DClSel_Gen->Fill(Radius_2DClSelGen);

	}// end loop1 to find seed

	float eta2dS = (*(multiClusters[i].begin() + cl2dSeed))->eta();
        float phi2dS = (*(multiClusters[i].begin() + cl2dSeed))->phi();

        h2_dR_SeedCl2_SelmultiCl_vsEta->Fill(std::abs(multiClusters[MC_best_index].eta()), 
					     reco::deltaR(multiClusters[MC_best_index].eta(), multiClusters[MC_best_index].phi(), eta2dS, phi2dS));
        tp_dR_SeedCl2_SelmultiCl_vsEta->Fill(std::abs(multiClusters[MC_best_index].eta()), 
					     reco::deltaR(multiClusters[MC_best_index].eta(), multiClusters[MC_best_index].phi(), eta2dS, phi2dS));
	//        h2_dR_SeedCl2_SelmultiCl_vsLayerS2d->Fill(layerSeedRH, reco::deltaR(multiClusters[MC_best_index].eta(), multiClusters[MC_best_index].phi(), eta2dS, phi2dS));
        h_dR_SeedCl2_SelmultiCl->Fill(reco::deltaR(multiClusters[MC_best_index].eta(), multiClusters[MC_best_index].phi(), eta2dS, phi2dS));



	//	h_2DClSelSumE_overMCE->Fill(sum2DClSelEnergy/multiClusters[i].energy());
        h_2DClSelSumE->Fill(sum2DClSelEnergy/it_caloPart->energy());

        for(int kk=0; kk<5; ++kk){
          tp_2DClSelSumE_vsDR_wrtGen->Fill(kk, sum2DClSelEnergy_vsDR_wrtGen[kk] / it_caloPart->energy());
          tp_2DClSelSumE_vsDR_wrtMCl->Fill(kk, sum2DClSelEnergy_vsDR_wrtMCl[kk] / it_caloPart->energy());
          tp_2DClSelSumE_vsRadius_wrtGen->Fill(kk, sum2DClSelEnergy_vsRadius_wrtGen[kk] / it_caloPart->energy());
          tp_2DClSelSumE_vsRadius_wrtMCl->Fill(kk, sum2DClSelEnergy_vsRadius_wrtMCl[kk] / it_caloPart->energy());
        }


	for(int il=0; il<52-2; ++il){
	  std::cout << " layer = " << il << " counts = " << layerCL[il] << std::endl;
	  if(layerCL[il] > 0 && layerCL[il+1] > 0 && layerCL[il+2] > 0){
	    areConsecutive = true;
	    break;
	  }
	}
	std::cout << " areConsecutive =  " << areConsecutive  << std::endl;
	if(areConsecutive){
	  for(reco::HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
	      it!=multiClusters[i].end(); it++){
	    sum2DClSelConsEnergy += (*it)->energy();
	    float Radius_2DClSelConsGen = utilsMet.dsGenRecoObj(etaGen, phiGen, (*it)->z(), (*it)->x(), (*it)->y());
	    float Radius_2DClSelConsMCl = utilsMet.dsGenRecoObj(multiClusters[MC_best_index].eta(), multiClusters[MC_best_index].phi(), (*it)->z(), (*it)->x(), (*it)->y());

	    if(int(Radius_2DClSelConsGen) < 5) sum2DClSelConsEnergy_vsRadius_wrtGen[int(Radius_2DClSelConsGen)] += (*it)->energy();
	    if(int(Radius_2DClSelConsMCl) < 5) sum2DClSelConsEnergy_vsRadius_wrtMCl[int(Radius_2DClSelConsMCl)] += (*it)->energy();
	  }
	}
	
      }

      if(MC_best_index == i){
	h_energyFractionInBestMC->Fill(multiClusters[i].energy()/it_caloPart->energy()); 
	h_n2dClPerMC->Fill(multiClusters[i].size());
	num2DClinBestMC = multiClusters[i].size();
	float MC_pT = multiClusters[i].energy()/cosh(multiClusters[i].eta());
	//	h2_dPhiMCtoGen_vsPt->Fill(MC_pT, reco::deltaPhi(multiClusters[i].phi(), it_caloPart->phi()));
	h2_dPhiMCtoGen_vsPt->Fill(it_caloPart->pt(), reco::deltaPhi(multiClusters[i].phi(), it_caloPart->phi()));
	h_dPhiMCtoGen->Fill(reco::deltaPhi(multiClusters[i].phi(), it_caloPart->phi()));
	
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
	  /*  
	  int rhSeed = 0;
	  //	  std::cout << " >>> 1st rechit loop " << std::endl;
	  const std::vector< std::pair<DetId, float> > &hf = (*it)->hitsAndFractions();
	  for(unsigned int j = 0; j < hf.size(); j++){
	    //here we loop over detid/fraction pairs
	    const HGCRecHit *hit = hitmap[hf[j].first];
	    //	    std::cout << " >>>  found " << std::endl;
	    if(hit->energy()*hf[j].second > hitmap[hf[rhSeed].first]->energy() * hf[rhSeed].second){
	      rhSeed = j;
	      layerSeedRH = recHitTools.getLayerWithOffset(hf[j].first);
	    }
	    //	    std::cout << " >>>  now time  " << std::endl;
	  }
	  //	  std::cout << " >>> 2nd rechit loop seedL = " << layerSeedRH << std::endl;

	  // 2nd loop to stick on the seed layer	
	  for(unsigned int j = 0; j < hf.size(); j++){
            const HGCRecHit *hit = hitmap[hf[j].first];
	    if(recHitTools.getLayerWithOffset(hf[j].first) == layerSeedRH){
	    
	    }	  
	  }
	  */
	

	  sum2DClusterEnergy += (*it)->energy();

	  float eta2DCl = (*it)->eta();
	  float phi2DCl = (*it)->phi();

	  Cluster2DZ.push_back(utilsMet.Ztolayer((*it)->z(), eta2DCl) );
	  Cluster2DIt.push_back(it - multiClusters[i].begin());

	  float dR_2DClGen = reco::deltaR(etaGen, eta2DCl, phiGen, phi2DCl);
	  float dR_2DClMCl = reco::deltaR(etaMC, eta2DCl, phiMC, phi2DCl);
	  
	  if(int(dR_2DClMCl) < 5) sum2DClusterEnergy_vsDR_wrtMCl[int(dR_2DClMCl)] += (*it)->energy();
	  if(int(dR_2DClGen) < 5) sum2DClusterEnergy_vsDR_wrtGen[int(dR_2DClGen)] += (*it)->energy();


	  float x2DCl = (*it)->x();
	  float y2DCl = (*it)->y();
	  float z2DCl = (*it)->z();

	  float Radius_2DClGen = utilsMet.dsGenRecoObj(etaGen, phiGen, z2DCl, x2DCl, y2DCl);
	  float Radius_2DClMCl = utilsMet.dsGenRecoObj(etaMC, phiMC, z2DCl, x2DCl, y2DCl);
	  
	  if(int(Radius_2DClGen) < 5) sum2DClusterEnergy_vsRadius_wrtGen[int(Radius_2DClGen)] += (*it)->energy();
	  if(int(Radius_2DClMCl) < 5) sum2DClusterEnergy_vsRadius_wrtMCl[int(Radius_2DClMCl)] += (*it)->energy();


	  h2_dPhi2DCltoGen_vsPt->Fill(it_caloPart->pt(), reco::deltaPhi(phiGen, phi2DCl));
	  h_dPhi2DCltoGen->Fill(reco::deltaPhi(phiGen, phi2DCl));

	  h_dR_2DCl_multiCl->Fill(Radius_2DClMCl);
	  h_dR_2DCl_Gen->Fill(Radius_2DClGen);


	} // loop over 2D cluster

	// FIXME LAYER dai recHit
	//    math::XYZPoint mclPosition = clusterTools->getMultiClusterPosition(multiClusters[i]);
	// std::cout << " mclPosition.eta = " << mclPosition.eta() << " mclPosition.phi = " << mclPosition.phi() 
	// 	      << " mcl.eta() = " << multiClusters[i].eta() << " phi = " << multiClusters[i].phi() << std::endl;


	float eta2dS = (*(multiClusters[i].begin() + cl2dSeed))->eta();
	float phi2dS = (*(multiClusters[i].begin() + cl2dSeed))->phi();


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


	h2_dR_SeedCl2_multiCl_vsEta->Fill(std::abs(etaMC), reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
	tp_dR_SeedCl2_multiCl_vsEta->Fill(std::abs(etaMC), reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
	//	h2_dR_SeedCl2_multiCl_vsLayerS2d->Fill(layerSeedRH, reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
	h_dR_SeedCl2_multiCl->Fill(reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));




	h_2DClSumE_overMCE->Fill(sum2DClusterEnergy/multiClusters[i].energy());
	h_2DClSumE->Fill(sum2DClusterEnergy/it_caloPart->energy());
	
	h2_2DClSumE_vsVtxZ->Fill(vz, sum2DClusterEnergy/it_caloPart->energy());
	tp_2DClSumE_vsVtxZ->Fill(vz, sum2DClusterEnergy/it_caloPart->energy()); 

	for(int kk=0; kk<5; ++kk){
	  tp_2DClSumE_vsDR_wrtGen->Fill(kk, sum2DClusterEnergy_vsDR_wrtGen[kk] / it_caloPart->energy());
	  tp_2DClSumE_vsDR_wrtMCl->Fill(kk, sum2DClusterEnergy_vsDR_wrtMCl[kk] / it_caloPart->energy());
	  tp_2DClSumE_vsRadius_wrtGen->Fill(kk, sum2DClusterEnergy_vsRadius_wrtGen[kk] / it_caloPart->energy());
	  tp_2DClSumE_vsRadius_wrtMCl->Fill(kk, sum2DClusterEnergy_vsRadius_wrtMCl[kk] / it_caloPart->energy());
	}
	//	std::cout << " reco::deltaR = " << reco::deltaR(etaMC, phiMC, eta2dS, phi2dS) << " ok = " << sqrt(pow(etaMC-eta2dS, 2 ) + pow(reco::deltaPhi(phiMC,phi2dS), 2) )<< std::endl;
      } // MC best
    } // loop MC

    h_nMCvsnum2dClBestMC->Fill(num2DClinBestMC, nmclusCount);  
    h_nMCPerEvt->Fill(nmclusCount);

    h_nSelMCvsnum2dClSelMC->Fill(num2DClinSelMC, nSelmclusCount);  
    h_nSelMCPerEvt->Fill(nSelmclusCount);

    //shower axis by recHits
    float axisX = 0;
    float axisY = 0;
    float axisZ = 0;
    float sumEnergyToNorm = 0;
    GlobalPoint showerAxis;

    float sumMtcRHEnergy = 0;
    float sumMtcRHEnergy_vsDR_wrtGen[5];
    float sumMtcRHEnergy_vsRadius_wrtGen[5];
    for(int ij=0; ij<5; ++ij){
      sumMtcRHEnergy_vsDR_wrtGen[ij] = 0;
      sumMtcRHEnergy_vsRadius_wrtGen[ij] = 0;
    }
    //loop on rechit - matched to gen
    const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
    for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin(); it_sc != simClusterRefVector.end(); ++it_sc) {
      const SimCluster simCluster = (*(*it_sc));
      const std::vector<std::pair<uint32_t,float> > hits_and_fractions = simCluster.hits_and_fractions();
      
      //first loop get shower axis from rechits
      for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
	//	unsigned int hitlayer = recHitTools.getLayerWithOffset(it_haf->first);
	DetId hitid = (it_haf->first);

	float rhEnergy = 0;

	std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmap.find(hitid);
	if(trovatore == hitmap.end()){
	  //	  std::cout << " fine " << std::endl;
	  continue;
	}
	else{
	  const HGCRecHit *hit = hitmap[hitid];
          rhEnergy = hit->energy();
	  sumEnergyToNorm += rhEnergy*it_haf->second;
	  axisX += rhEnergy*it_haf->second * recHitTools.getPosition(hitid).x();
          axisY += rhEnergy*it_haf->second * recHitTools.getPosition(hitid).y();
          axisZ += rhEnergy*it_haf->second * recHitTools.getPosition(hitid).z();
	}
      }
      showerAxis = GlobalPoint(axisX/sumEnergyToNorm, axisY/sumEnergyToNorm, (axisZ - vz)/sumEnergyToNorm);

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
      for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
	DetId hitid = (it_haf->first);

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
	  if(hit->energy()*it_haf->second > allRHSeedEnergy){
	    allRHSeedEnergy = hit->energy()*it_haf->second;
	    //	    allRHSeedLayer = recHitTools.getLayerWithOffset(hitid);
	  }
	}
	if(found){
	  sumMtcRHEnergy += rhEnergy*it_haf->second;
	  float rhEta = recHitTools.getEta(hitid);
	  float rhPhi = recHitTools.getPhi(hitid);
	  
	  h2_dPhivsdEta_rhGen->Fill(reco::deltaPhi(rhPhi, it_caloPart->phi()), rhEta - it_caloPart->eta());
	  h2_dPhivsdEta_rhAxis->Fill(reco::deltaPhi(rhPhi, axPhi), rhEta - axEta);
	  h2_dPhivsdEta_GenAxis->Fill(reco::deltaPhi(it_caloPart->phi(), axPhi), it_caloPart->eta() - axEta);

	  float rhX = recHitTools.getPosition(hitid).x();
	  float rhY = recHitTools.getPosition(hitid).y();
	  float rhL = recHitTools.getLayerWithOffset(hitid);
	  float rhZ = utilsMet.layerToZ(rhL, rhEta);

	  float dR_2DClGen = sqrt( pow(rhX - axX, 2) + pow(rhY - axY, 2) + pow(rhZ - axZ, 2) );
	  if(int(dR_2DClGen) < 5) sumMtcRHEnergy_vsDR_wrtGen[int(dR_2DClGen)] += rhEnergy*it_haf->second;
	  
	  //	  float dEP_2DClGen = sqrt( pow(rhEta - axEta, 2) + pow(reco::deltaPhi(rhPhi, axPhi), 2));
	  
	  float Radius_2DClGen = utilsMet.dsGenRecHit(etaGen, phiGen, rhL, rhX, rhY);
	  if(int(Radius_2DClGen) < 5) sumMtcRHEnergy_vsRadius_wrtGen[int(Radius_2DClGen)] += rhEnergy*it_haf->second;	  


	  h2_dPhiMatchHittoGen_vsPt->Fill(it_caloPart->pt(), reco::deltaPhi(rhPhi, it_caloPart->phi()));
	  h_dPhiMatchHittoGen->Fill(reco::deltaPhi(rhPhi, it_caloPart->phi()));
	  h_dR_MatchHit_Gen->Fill(Radius_2DClGen);
	}
      }// first loop over rechits


      //second to select seed layer => not needed
      
    }

    h_mtcRHSumE -> Fill(sumMtcRHEnergy/it_caloPart->energy());
    for(int ij=0; ij<5; ++ij){
      tp_mtcRHSumE_vsDR_wrtGen->Fill(ij, sumMtcRHEnergy_vsDR_wrtGen[ij]/it_caloPart->energy());
      tp_mtcRHSumE_vsRadius_wrtGen->Fill(ij, sumMtcRHEnergy_vsRadius_wrtGen[ij]/it_caloPart->energy());
    }

  }//caloparticle

}

void
HGCalAxisAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalAxisAnalyzer::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalAxisAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalAxisAnalyzer);
