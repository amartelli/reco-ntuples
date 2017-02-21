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

#include "HGCalCalibration/HitValidation/interface/UsefulClasses.h"
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
  //  edm::EDGetTokenT<reco::CaloClusterCollection> _clusters;
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


  std::vector<int> nHits_layer;
  std::vector<int> nHitsWithTime_layer;

  //  TH1F* h_numberOfMC;
  TH1F* h_energyFractionInMC;
  //  TH2F* h2_energyFractionInMC_vs_VtxZ;
  TH2F* h2_dPhiMCtoGen_vsPt;
  TH1F* h_dPhiMCtoGen;
  TH1F* h_dR_SeedCl2_multiCl;
  TH2F* h2_dR_SeedCl2_multiCl_vsLayerS2d;
  TH2F* h2_dR_SeedCl2_multiCl_vsEta;
  TProfile* tp_dR_SeedCl2_multiCl_vsEta;
  TH1F* h_n2dClPerMC;
  TH1F* h_nMCPerEvt;

  TH2F* h2_cellRadius_vsThick;

  TH1F* h_EfracMCl_wrt_Best;

  //  TH1F* h_allhitFractionRHWithTimeinRing_wrtGen[6];
  TH1F* h_allhitAverageTimeinRing_wrtGen[6];
  TH1F* h_allhitTimeinRing_wrtGen[6];
  TH1F* h_allhitTimeinRing_wrtGen_fixL[6];

  TH1F* h_hitTimeinRing_wrtGen[6];
  TH1F* h_hitTimeinRing_wrtMCl[6];
  TH2F* h2_hitTimeinRing_wrtGen_vsL;

  TH1F* h_hitTimeinRing_wrtGen_fixL[6];
  TH1F* h_hitTimeinRing_wrtMCl_fixL[6];
  TH2F* h2_hitTimeinRing_wrtGen_vsL_fixL;

  TH1F* h_allRH_TimesOfRadius;
  TH1F* h_allRH_TimesOfRadius_fixL;

  TH2F* h2_allRH_TimesVsDRaxis;
  TProfile* tp_allRH_TimesVsDRaxis;
  TH2F* h2_allRH_TimesVsDEPaxis;
  TProfile* tp_allRH_TimesVsDEPaxis;
  TH2F* h2_allRH_TimesVsDRaxis_fixL;
  TProfile* tp_allRH_TimesVsDRaxis_fixL;


  TProfile2D* dPhivsDEta_Time;
  TProfile2D* dRvsDZ_Time;

  TH2F* h2_allRH_TimesVsRadiusaxis;
  TProfile* tp_allRH_TimesVsRadiusaxis;
  TH2F* h2_allRH_TimesVsRadiusaxis_fixL;
  TProfile* tp_allRH_TimesVsRadiusaxis_fixL;

  TH1F* h_2DClSumE;
  TProfile* tp_2DClSumE_vsDR_wrtGen;
  TProfile* tp_2DClSumE_vsDR_wrtMCl;
  TProfile* tp_2DClSumE_vsRadius_wrtGen;
  TProfile* tp_2DClSumE_vsRadius_wrtMCl;

  TH1F* h_mtcRHSumE;
  TProfile* tp_mtcRHSumE_vsDR_wrtGen;
  //  TProfile* tp_mtcRHSumE_vsDR_wrtMCl;
  TProfile* tp_mtcRHSumE_vsRadius_wrtGen;
  //  TProfile* tp_mtcRHSumE_vsRadius_wrtMCl;

  TH1F* h_dEta_pos;
  TH1F* h_dPhi_pos;
  TH1F* h_dEta_neg;
  TH1F* h_dPhi_neg;

  TH1F* h_numRhPerEvt_vsdR;
  TH1F* h_rHGen_dR;
  TH2F* h2_rHGen_dR_vsEta;
  TH1F* h_rHGen_dDist;
  TH2F* h2_rHGen_dDist_vsEta;
  TProfile* tp_rHGen_dDist_vsEta;
  TH2F* h2_Energia_vsEta;
  TH2F* h2_Pt_vsEta;
  TH1F* h_rHGen_dR_fixL;
  TH2F* h2_rHGen_dR_vsEnergia;
  TH2F* h2_rHGen_dR_vsPt;
  TH2F* h2_rHGen_dRcosh_vsEta;
  TH2F* h2_rHGen_dRcoshP_vsEta;

  TProfile2D* tp2_Energy_vs_Eta_dR;
  TProfile2D* tp2_Energy_vs_Eta_dDist;


  TH2F* h2_dPhivsdEta_rhGen;
  TH2F* h2_dPhivsdEta_rhGen_fixL;
  TH2F* h2_dPhivsdEta_rhAxis;
  TH2F* h2_dPhivsdEta_GenAxis;

  //eta bins 1.7-1.9   1.9-2.1 - 2.1-2.3   2.3-2.5 
  TH1F* h_rhGen_Radius_Eta[8];
  TH1F* h_rhGen_RadiusDdeta_Eta[8];
  TH1F* h_rhGen_RadiusPdeta_Eta[8];
  TH1F* hAverageTime_Eta_dRadius[8][3];
  TH1F* hAverageTime_Eta_dRadius_fixL[8][2];

  TH1F* hTotHits_Eta_dRadius[8][3];
  TH1F* hTotHitsWithTime_Eta_dRadius[8][3];
  TH1F* hFractionHitsWithTime_Eta_dRadius[8][3];
  TH1F* hFractionEvents_HitsWithTime_Eta_dRadius[8][3];
  TProfile* tpFractionHitsWithTime_Eta_dRadius_vsPt[8][3];
  TH1F* hTotHits_Eta_dRadius_fixL[8][2];
  TH1F* hTotHitsWithTime_Eta_dRadius_fixL[8][2];
  TH1F* hFractionHitsWithTime_Eta_dRadius_fixL[8][2];


  int totEvtsEtaRadius[8][3];
  int totEvtsEtaRadius_withTime[8][3];


  float radiusEtaRad[8][3];
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
  //  _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("hgcalLayerClusters"));
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

  //  LayerOccupancy = fs->make<TH1F>("LayerOccupancy", "", 60, 0., 60.);

  //  h_numberOfMC = fs->make<TH1F>("h_numberOfMC", "", 60, 0., 60.);
  h_energyFractionInMC = fs->make<TH1F>("h_energyFractionInMC", "", 1000, 0., 2.);
  //  h2_energyFractionInMC_vs_VtxZ = fs->make<TH2F>("h2_energyFractionInMC_vs_VtxZ", "", 1000, -1., 1., 100, 0., 1.);
  h2_dPhiMCtoGen_vsPt = fs->make<TH2F>("h2_dPhiMCtoGen_vsPt", "", 300, 0., 300., 600, -0.3, 0.3);
  h_dPhiMCtoGen = fs->make<TH1F>("h_dPhiMCtoGen", "", 600, -3., 3.);

  h_dR_SeedCl2_multiCl = fs->make<TH1F>("h_dR_SeedCl2_multiCl", "", 1000, 0., 3.);
  h2_dR_SeedCl2_multiCl_vsLayerS2d = fs->make<TH2F>("h2_dR_SeedCl2_multiCl_vsLayerS2d", "", 100, 0., 100., 1000, 0., 5.);
  h2_dR_SeedCl2_multiCl_vsEta = fs->make<TH2F>("h2_dR_SeedCl2_multiCl_vsEta", "", 500, 1.5, 3.5, 1000, 0., 5.);
  tp_dR_SeedCl2_multiCl_vsEta = fs->make<TProfile>("tp_dR_SeedCl2_multiCl_vsEta", "", 500, 1.5, 3.5);
  h_n2dClPerMC = fs->make<TH1F>("h_n2dClPerMC", "", 50, 0., 50.);
  h_nMCPerEvt = fs->make<TH1F>("h_nMCPerEvt", "", 100, 0., 100.);

  h2_cellRadius_vsThick = fs->make<TH2F>("h2_cellRadius_vsThick", "", 500., 0., 500., 100., 0., 10.);

  h_EfracMCl_wrt_Best = fs->make<TH1F>("h_EfracMCl_wrt_Best", "", 500, 0., 1.1);

  for (int ij=0; ij<6; ++ij){
    h_hitTimeinRing_wrtGen[ij] = fs->make<TH1F>(Form("h_hitTimeinRing_wrtGen_%d", ij), "", 70., -0.2, 0.5);
    h_hitTimeinRing_wrtMCl[ij] = fs->make<TH1F>(Form("h_hitTimeinRing_wrtMCl_%d", ij), "", 70., -0.2, 0.5);
    h_hitTimeinRing_wrtGen_fixL[ij] = fs->make<TH1F>(Form("h_hitTimeinRing_wrtGen_fixL_%d", ij), "", 70., -0.2, 0.5);
    h_hitTimeinRing_wrtMCl_fixL[ij] = fs->make<TH1F>(Form("h_hitTimeinRing_wrtMCl_fixL_%d", ij), "", 70., -0.2, 0.5);
    h_allhitTimeinRing_wrtGen[ij] = fs->make<TH1F>(Form("h_allhitTimeinRing_wrtGen_%d", ij), "", 70., -0.2, 0.5);
    h_allhitTimeinRing_wrtGen_fixL[ij] = fs->make<TH1F>(Form("h_allhitTimeinRing_wrtGen_fixL_%d", ij), "", 70., -0.2, 0.5);
    //    h_allhitFractionRHWithTimeinRing_wrtGen[ij] = fs->make<TH1F>(Form("h_allhitFractionRHWithTimeinRing_wrtGen_%d", ij), "", 70., -0.2, 0.5);
    h_allhitAverageTimeinRing_wrtGen[ij] = fs->make<TH1F>(Form("h_allhitAverageTimeinRing_wrtGen_fixL_%d", ij), "", 70., -0.2, 0.5);
  }

  h2_hitTimeinRing_wrtGen_vsL = fs->make<TH2F>("h2_hitTimeinRing_wrtGen_vsL", "", 60., 0, 60., 70., -0.2, 0.5);
  h2_hitTimeinRing_wrtGen_vsL_fixL = fs->make<TH2F>("h2_hitTimeinRing_wrtGen_vsL_fixL", "", 60., 0, 60., 70., -0.2, 0.5);

  h_allRH_TimesOfRadius = fs->make<TH1F>("h_allRH_TimesOfRadius", "", 100., 0., 50.);
  h_allRH_TimesOfRadius_fixL = fs->make<TH1F>("h_allRH_TimesOfRadius_fixL", "", 100., 0., 50.);
  //  std::cout << " >>> costruito " << std::endl;

  h2_allRH_TimesVsDRaxis = fs->make<TH2F>("h2_allRH_TimesVsDRaxis", "", 500, 0., 10., 70., -0.2, 0.5);;
  tp_allRH_TimesVsDRaxis = fs->make<TProfile>("tp_allRH_TimesVsDRaxis", "", 500, 0., 10.);
  h2_allRH_TimesVsDEPaxis = fs->make<TH2F>("h2_allRH_TimesVsDEPaxis", "", 500, 0., 10., 70., -0.2, 0.5);;
  tp_allRH_TimesVsDEPaxis = fs->make<TProfile>("tp_allRH_TimesVsDEPaxis", "", 500, 0., 10.);
  h2_allRH_TimesVsDRaxis_fixL = fs->make<TH2F>("h2_allRH_TimesVsDRaxis_fixL", "", 500, 0., 10., 70., -0.2, 0.5);;
  tp_allRH_TimesVsDRaxis_fixL = fs->make<TProfile>("tp_allRH_TimesVsDRaxis_fixL", "", 500, 0., 10.);
  h2_allRH_TimesVsRadiusaxis = fs->make<TH2F>("h2_allRH_TimesVsRadiusaxis", "", 500, 0., 10., 70., -0.2, 0.5);;
  tp_allRH_TimesVsRadiusaxis = fs->make<TProfile>("tp_allRH_TimesVsRadiusaxis", "", 500, 0., 10.);
  h2_allRH_TimesVsRadiusaxis_fixL = fs->make<TH2F>("h2_allRH_TimesVsRadiusaxis_fixL", "", 500, 0., 10., 70., -0.2, 0.5);;
  tp_allRH_TimesVsRadiusaxis_fixL = fs->make<TProfile>("tp_allRH_TimesVsRadiusaxis_fixL", "", 500, 0., 10.);

  dPhivsDEta_Time = fs->make<TProfile2D>("dPhivsDEta_Time", "", 500, -0.2, 0.2, 500, -0.2, 0.2);
  dRvsDZ_Time = fs->make<TProfile2D>("dRvsDZ_Time", "", 500, 0., 5., 500., 0., 100.);



  h_2DClSumE = fs->make<TH1F>("h_2DClSumE", "", 500, 0., 2.);
  tp_2DClSumE_vsDR_wrtGen = fs->make<TProfile>("tp_2DClSumE_vsDR_wrtGen", "", 5, 0., 5);
  tp_2DClSumE_vsDR_wrtMCl = fs->make<TProfile>("tp_2DClSumE_vsDR_wrtMCl", "", 5, 0., 5);
  tp_2DClSumE_vsRadius_wrtGen = fs->make<TProfile>("tp_2DClSumE_vsRadius_wrtGen", "", 5, 0., 5);
  tp_2DClSumE_vsRadius_wrtMCl = fs->make<TProfile>("tp_2DClSumE_vsRadius_wrtMCl", "", 5, 0., 5);


  h_mtcRHSumE = fs->make<TH1F>("h_mtcRHSumE", "", 500, 0., 2.);
  tp_mtcRHSumE_vsDR_wrtGen = fs->make<TProfile>("tp_mtcRHSumE_vsDR_wrtGen", "", 5, 0., 5);
  //  le* tp_mtcRHSumE_vsDR_wrtMCl
  tp_mtcRHSumE_vsRadius_wrtGen = fs->make<TProfile>("tp_mtcRHSumE_vsRadius_wrtGen", "", 5, 0., 5);
  //  le* tp_mtcRHSumE_vsRadius_wrtMCl;


  h_dEta_pos = fs->make<TH1F>("h_dEta_pos", "", 500, -0.1, 0.1);
  h_dPhi_pos = fs->make<TH1F>("h_dPhi_pos", "", 500, -0.1, 0.1);
  h_dEta_neg = fs->make<TH1F>("h_dEta_neg", "", 500, -0.1, 0.1);
  h_dPhi_neg = fs->make<TH1F>("h_dPhi_neg", "", 500, -0.1, 0.1);


  h_rHGen_dR = fs->make<TH1F>("h_rHGen_dR", "", 1000., 0., 5.);
  h_rHGen_dDist = fs->make<TH1F>("h_rHGen_dDist", "", 1000., 0., 10.);
  h2_rHGen_dR_vsEta = fs->make<TH2F>("h2_rHGen_dR_vsEta", "", 500., 1.5, 3., 1000., 0., 5.);
  h2_rHGen_dDist_vsEta = fs->make<TH2F>("h2_rHGen_dDist_vsEta", "", 500., 1.5, 3., 1000., 0., 10.);
  tp_rHGen_dDist_vsEta = fs->make<TProfile>("tp_rHGen_dDist_vsEta", "", 500., 1.5, 3.);
  h_rHGen_dR_fixL = fs->make<TH1F>("h_rHGen_dR_fixL", "", 1000., 0., 5.);
  h2_Energia_vsEta = fs->make<TH2F>("h2_Energia_vsEta", "", 500., 1.5, 3., 1000, 0., 1000);
  h2_Pt_vsEta = fs->make<TH2F>("h2_Pt_vsEta", "", 500., 1.5, 3., 1000, 0., 200.);
  h2_rHGen_dR_vsEnergia = fs->make<TH2F>("h2_rHGen_dR_vsEnergia", "", 1000., 0., 1000., 1000., 0., 5.);
  h2_rHGen_dR_vsPt = fs->make<TH2F>("h2_rHGen_dR_vsPt", "", 1000., 0., 200., 1000., 0., 5.);
  h2_rHGen_dRcosh_vsEta = fs->make<TH2F>("h2_rHGen_dRcosh_vsEta", "", 500., 1.5, 3., 1000., 0., 5.);
  h2_rHGen_dRcoshP_vsEta = fs->make<TH2F>("h2_rHGen_dRcoshP_vsEta", "", 500., 1.5, 3., 1000., 0., 5.);

  tp2_Energy_vs_Eta_dR = fs->make<TProfile2D>("tp2_Energy_vs_Eta_dR", "", 500., 1.5, 3., 1000, 0., 1000);
  tp2_Energy_vs_Eta_dDist = fs->make<TProfile2D>("tp2_Energy_vs_Eta_dDist", "", 500., 1.5, 3., 1000, 0., 1000);

  h2_dPhivsdEta_rhGen = fs->make<TH2F>("h2_dPhivsdEta_rhGen", "", 500, -0.2, 0.2, 500, -0.2, 0.2);
  h2_dPhivsdEta_rhGen_fixL = fs->make<TH2F>("h2_dPhivsdEta_rhGen_fixL", "", 500, -0.2, 0.2, 500, -0.2, 0.2);
  h2_dPhivsdEta_rhAxis = fs->make<TH2F>("h2_dPhivsdEta_rhAxis", "", 1000, -2., 2., 1000, -2., 2.);
  h2_dPhivsdEta_GenAxis = fs->make<TH2F>("h2_dPhivsdEta_GenAxis", "", 1000, -2., 2., 1000, -2., 2.);

  for(int ieta=0; ieta<8; ++ieta){
    
    std::cout << " ieta from = " << (1.5+ieta*0.2) << " to " << 1.5+0.2+ieta*0.2 << std::endl;
    h_rhGen_Radius_Eta[ieta] = fs->make<TH1F>(Form("h_rhGen_Radius_Eta_%.1f-%.1f", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2), "", 1000, 0., 10.); 
    h_rhGen_RadiusDdeta_Eta[ieta] = fs->make<TH1F>(Form("h_rhGen_RadiusDdeta_Eta_%.1f-%.1f", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2), "", 1000, 0., 10.); 
    h_rhGen_RadiusPdeta_Eta[ieta] = fs->make<TH1F>(Form("h_rhGen_RadiusPdeta_Eta_%.1f-%.1f", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2), "", 1000, 0., 10.); 

    for(int iRad=0; iRad<3; ++iRad){
      hAverageTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hAverageTime_Eta%.1f-%.1f_dRadius%d", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 500, -0.2, 0.5);
      hFractionHitsWithTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionHits_Eta%.1f-%.1f_dRadius%d", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 1000, 0, 1.);
      hFractionEvents_HitsWithTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_Hits_Eta%.1f-%.1f_dRadius%d", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 1000, 0, 1.);
      tpFractionHitsWithTime_Eta_dRadius_vsPt[ieta][iRad] = fs->make<TProfile>(Form("tpFractionHits_Eta%.1f-%.1f_dRadius%d_vsPt", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 1000, 0, 1000.);
      hTotHits_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hTotHits_Eta_%.1f-%.1f_dRadius%d", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 1000, 0, 100.);
      hTotHitsWithTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hTotHitsWithTime_Eta_%.1f-%.1f_dRadius%d", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 1000, 0, 100.);

      if(iRad > 1) continue;
      hAverageTime_Eta_dRadius_fixL[ieta][iRad] = fs->make<TH1F>(Form("hAverageTime_Eta%.1f-%.1f_dRadius%d_fixL", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 500, -0.2, 0.5);
      hFractionHitsWithTime_Eta_dRadius_fixL[ieta][iRad] = fs->make<TH1F>(Form("hFractionHits_Eta%.1f-%.1f_dRadius%d_fixL", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 1000, 0, 1.);
      hTotHits_Eta_dRadius_fixL[ieta][iRad] = fs->make<TH1F>(Form("hTotHits_Eta_%.1f-%.1f_dRadius%d_fixL", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 1000, 0, 100.);
      hTotHitsWithTime_Eta_dRadius_fixL[ieta][iRad] = fs->make<TH1F>(Form("hTotHitsWithTime_Eta_%.1f-%.1f_dRadius%d_fixL", (1.5+ieta*0.2), 1.5+0.2+ieta*0.2, iRad), "", 1000, 0, 100.);

    }
  }


  for(int yui=0; yui<8; ++yui){
    for(int wer=0; wer<3; ++wer)
      {
	totEvtsEtaRadius[yui][wer] = 0;
	totEvtsEtaRadius_withTime[yui][wer] = 0;
      }
  }


  for(int ij=0; ij<8; ++ij){
    radiusEtaRad[ij][0] = 2.;
    radiusEtaRad[ij][1] = 5.;
    radiusEtaRad[ij][2] = 10.;
  }

  /*
  radiusEtaRad[0][0] = 1.55;
  radiusEtaRad[0][1] = 2.15;
  radiusEtaRad[0][2] = 2.55;
  radiusEtaRad[1][0] = 1.4;
  radiusEtaRad[1][1] = 1.8;
  radiusEtaRad[1][2] = 2.25;
  radiusEtaRad[2][0] = 1.1;
  radiusEtaRad[2][1] = 1.6;
  radiusEtaRad[2][2] = 1.95;
  radiusEtaRad[3][0] = 1.15;
  radiusEtaRad[3][1] = 1.35;
  radiusEtaRad[3][2] = 1.7;
  radiusEtaRad[4][0] = 1.;
  radiusEtaRad[4][1] = 1.35;
  radiusEtaRad[4][2] = 1.7;
  radiusEtaRad[5][0] = 0.85;
  radiusEtaRad[5][1] = 1.1;
  radiusEtaRad[5][2] = 1.4;
  radiusEtaRad[6][0] = 0.65;
  radiusEtaRad[6][1] = 0.85;
  radiusEtaRad[6][2] = 1.05;
  radiusEtaRad[7][0] = 0.55;
  radiusEtaRad[7][1] = 0.7;
  radiusEtaRad[7][2] = 0.85;
  */
/*
  radiusEtaRad[0][0] = 1.7;
  radiusEtaRad[0][1] = 2.3;
  radiusEtaRad[0][2] = 3.;
  radiusEtaRad[1][0] = 1.6;
  radiusEtaRad[1][1] = 2.2;
  radiusEtaRad[1][2] = 2.9;
  radiusEtaRad[2][0] = 1.5;
  radiusEtaRad[2][1] = 2.1;
  radiusEtaRad[2][2] = 2.7;
  radiusEtaRad[3][0] = 1.3;
  radiusEtaRad[3][1] = 1.9;
  radiusEtaRad[3][2] = 2.4;
  radiusEtaRad[4][0] = 1.;
  radiusEtaRad[4][1] = 1.6;
  radiusEtaRad[4][2] = 2.1;
  radiusEtaRad[5][0] = 0.9;
  radiusEtaRad[5][1] = 1.4;
  radiusEtaRad[5][2] = 1.8;
  radiusEtaRad[6][0] = 0.85;
  radiusEtaRad[6][1] = 1.2;
  radiusEtaRad[6][2] = 1.6;
  radiusEtaRad[7][0] = 0.8;
  radiusEtaRad[7][1] = 1.1;
  radiusEtaRad[7][2] = 1.5;
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
  /*
  agpc->clear();
  arhc->clear();
  arhc_raw->clear();
  acdc->clear();
  amcc->clear();
  ascc->clear();
  apfcc->clear();
  acpc->clear();
  */

  // Energy_layer_calib.clear();
  // Energy_layer_calib_fraction.clear();
  // nHits_layer.clear();
  // nHitsWithTime_layer.clear();
  // for(unsigned int ij=0; ij<60; ++ij){
  //   Energy_layer_calib.push_back(0.);
  //   Energy_layer_calib_fraction.push_back(0.);
  //   nHits_layer.push_back(0.);
  //   nHitsWithTime_layer.push_back(0.);
  // }

  recHitTools.getEventSetup(iSetup);
  //  utilsMet = UsefulClasses();

  //  int npart = 0;
  /*
  int nhit  = 0;
  int nhit_raw = 0;
  int nsimclus = 0;
  int ncalopart = 0;
  */

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

  // Handle<reco::CaloClusterCollection> clusterHandle;
  // iEvent.getByToken(_clusters,clusterHandle);
  // const reco::CaloClusterCollection &clusters = *clusterHandle;

  // edm::PtrVector<reco::CaloCluster> clusterPtrs;
  // for( unsigned i = 0; i < clusterHandle->size(); ++i ){
  //   edm::Ptr<reco::CaloCluster> ptr(clusterHandle,i);
  //   clusterPtrs.push_back(ptr);
  // }

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

  // loop over caloParticles
  for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){

    float etaGen = it_caloPart->eta();
    float phiGen = it_caloPart->phi();

    UsefulClasses utilsMet = UsefulClasses(etaGen, phiGen); 
      

    int nmclusCount = 0;
    unsigned int MC_best_index = 0;
    float MC_best_energy = 0;
    //    std::cout << " >>> now MCL " << std::endl;
    for(unsigned int i = 0; i < multiClusters.size(); i++){
      if(multiClusters[i].eta() * it_caloPart->eta() < 0) continue;

      if(multiClusters[i].energy() > MC_best_energy){
	MC_best_index = i;
	MC_best_energy = multiClusters[i].energy();
      }
    }

    float sum2DClusterEnergy = 0;
    float sum2DClusterEnergy_vsDR_wrtGen[5];
    float sum2DClusterEnergy_vsDR_wrtMCl[5];
    float sum2DClusterEnergy_vsRadius_wrtGen[5];
    float sum2DClusterEnergy_vsRadius_wrtMCl[5];
    for(int ij=0; ij<5; ++ij){
      sum2DClusterEnergy_vsDR_wrtGen[ij] = 0;
      sum2DClusterEnergy_vsRadius_wrtGen[ij] = 0;
      sum2DClusterEnergy_vsDR_wrtMCl[ij] = 0;
      sum2DClusterEnergy_vsRadius_wrtMCl[ij] = 0;
    }
    //      std::cout << " >>> MC energy = " << multiClusters[i].energy() << " i  = " << i << std::endl;
    for(unsigned int i = 0; i < multiClusters.size(); i++){
      if(multiClusters[i].eta() * it_caloPart->eta() < 0) continue;

      float etaMC = multiClusters[i].eta();
      float phiMC = multiClusters[i].phi();


      h_EfracMCl_wrt_Best->Fill(multiClusters[i].energy() / multiClusters[MC_best_index].energy());

      h_energyFractionInMC->Fill(multiClusters[i].energy()/it_caloPart->energy()); 

      if(multiClusters[i].energy() >= 0.10 * multiClusters[MC_best_index].energy()) ++nmclusCount;

      if(MC_best_index == i){
	h_n2dClPerMC->Fill(multiClusters[i].size());
	float MC_pT = multiClusters[i].energy()/cosh(multiClusters[i].eta());
	h2_dPhiMCtoGen_vsPt->Fill(MC_pT, reco::deltaPhi(multiClusters[i].phi(), it_caloPart->phi()));
	h_dPhiMCtoGen->Fill(reco::deltaPhi(multiClusters[i].phi(), it_caloPart->phi()));
	
	int cl2dSeed = 0;
	//	int cl2dSeedRa = 0;
	unsigned int layerSeedRH = 0;
	//	std::cout << " >>> nuovo ciclo " << std::endl;
	//	std::cout << " >>> now 2D cluster " << std::endl;
	for(reco::HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
	    it!=multiClusters[i].end(); it++){

	  // std::cout << " >>> now get the seed 2dcluster (*it)->energy() = " << (*it)->energy() << "  of seed = " << (*(multiClusters[i].begin()+cl2dSeedRa))->energy() << std::endl;
	  // std::cout << " >>> cl2dSeed = " << cl2dSeed << std::endl;

	  // if((*it)->energy() > (*(it+cl2dSeed))->energy()){
	  //   cl2dSeed = it - multiClusters[i].begin();
	  //   std::cout << " clem = " << (*it)->energy() << std::endl;
	  // }
	  if((*it)->energy() > (*(multiClusters[i].begin()+cl2dSeed))->energy()){
	    cl2dSeed = it - multiClusters[i].begin();
	    //	    std::cout << " Ra = " << (*it)->energy() << std::endl;
	  }
	  
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
	    if(hit->time() > 0){
	      float rhX = recHitTools.getPosition(hf[j].first).x();
	      float rhY = recHitTools.getPosition(hf[j].first).y();
	      float rhL = recHitTools.getLayerWithOffset(hf[j].first);

	      float dRwrtGen = utilsMet.dsGenRecHit(etaGen, phiGen, rhL, rhX, rhY);
	      float dRwrtMCl = utilsMet.dsGenRecHit(etaMC, phiMC, rhL, rhX, rhY);

	      h2_cellRadius_vsThick->Fill(recHitTools.getSiThickness(hf[j].first), recHitTools.getRadiusToSide(hf[j].first));

	      //	      int timesOfRadius = int(dRwrtGen / recHitTools.getRadiusToSide(hf[j].first) );
	      int timesOfRadius = int(dRwrtGen);

	      if(timesOfRadius < 3) h_hitTimeinRing_wrtGen[timesOfRadius]->Fill(hit->time()-1);
	      else h_hitTimeinRing_wrtGen[3]->Fill(hit->time()-1);
	      if(timesOfRadius < 2) h2_hitTimeinRing_wrtGen_vsL->Fill(recHitTools.getLayerWithOffset(hf[j].first), hit->time()-1);

	      //	      timesOfRadius = int(dRwrtMCl / recHitTools.getRadiusToSide(hf[j].first));
	      timesOfRadius = int(dRwrtMCl);
	      if(timesOfRadius < 3) h_hitTimeinRing_wrtMCl[timesOfRadius]->Fill(hit->time()-1);
	      else h_hitTimeinRing_wrtMCl[3]->Fill(hit->time()-1);
	    }
	  }
	  //	  std::cout << " >>> 2nd rechit loop seedL = " << layerSeedRH << std::endl;
	  // 2nd loop to stick on the seed layer
	
	  for(unsigned int j = 0; j < hf.size(); j++){
            const HGCRecHit *hit = hitmap[hf[j].first];
	    if(recHitTools.getLayerWithOffset(hf[j].first) == layerSeedRH){
	      if(hit->time() > 0){
		float rhX = recHitTools.getPosition(hf[j].first).x();
		float rhY = recHitTools.getPosition(hf[j].first).y();
		float rhL = recHitTools.getLayerWithOffset(hf[j].first);
		
		float dRwrtGen = utilsMet.dsGenRecHit(etaGen, phiGen, rhL, rhX, rhY);
		float dRwrtMCl = utilsMet.dsGenRecHit(etaMC, phiMC, rhL, rhX, rhY);

		//		int timesOfRadius = int(dRwrtGen / recHitTools.getRadiusToSide(hf[j].first));
		int timesOfRadius = int(dRwrtGen);
		if(timesOfRadius < 5) h_hitTimeinRing_wrtGen_fixL[timesOfRadius]->Fill(hit->time()-1);
		else h_hitTimeinRing_wrtGen_fixL[5]->Fill(hit->time()-1);
		if(timesOfRadius < 3) h2_hitTimeinRing_wrtGen_vsL_fixL->Fill(recHitTools.getLayerWithOffset(hf[j].first), hit->time()-1);

		//		timesOfRadius = int(dRwrtMCl / recHitTools.getRadiusToSide(hf[j].first));
		timesOfRadius = int(dRwrtMCl);
		if(timesOfRadius < 5) h_hitTimeinRing_wrtMCl_fixL[timesOfRadius]->Fill(hit->time()-1);
		else h_hitTimeinRing_wrtMCl_fixL[5]->Fill(hit->time()-1);
	      }
	    }	  
	  }
	

	  sum2DClusterEnergy += (*it)->energy();

	  float eta2DCl = (*it)->eta();
	  float phi2DCl = (*it)->phi();

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
	} // loop over 2D cluster

	// FIXME LAYER dai recHit
	//    math::XYZPoint mclPosition = clusterTools->getMultiClusterPosition(multiClusters[i]);
	// std::cout << " mclPosition.eta = " << mclPosition.eta() << " mclPosition.phi = " << mclPosition.phi() 
	// 	      << " mcl.eta() = " << multiClusters[i].eta() << " phi = " << multiClusters[i].phi() << std::endl;

	float eta2dS = (*(multiClusters[i].begin() + cl2dSeed))->eta();
	float phi2dS = (*(multiClusters[i].begin() + cl2dSeed))->phi();

	h2_dR_SeedCl2_multiCl_vsEta->Fill(std::abs(etaMC), reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
	tp_dR_SeedCl2_multiCl_vsEta->Fill(std::abs(etaMC), reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
	h2_dR_SeedCl2_multiCl_vsLayerS2d->Fill(layerSeedRH, reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));
	h_dR_SeedCl2_multiCl->Fill(reco::deltaR(etaMC, phiMC, eta2dS, phi2dS));

	//	std::cout << " reco::deltaR = " << reco::deltaR(etaMC, phiMC, eta2dS, phi2dS) << " ok = " << sqrt(pow(etaMC-eta2dS, 2 ) + pow(reco::deltaPhi(phiMC,phi2dS), 2) )<< std::endl;
      } // MC best
    } // loop MC

    h_2DClSumE->Fill(sum2DClusterEnergy/it_caloPart->energy());
    for(int kk=0; kk<5; ++kk){
      tp_2DClSumE_vsDR_wrtGen->Fill(kk, sum2DClusterEnergy_vsDR_wrtGen[kk] / it_caloPart->energy());
      tp_2DClSumE_vsDR_wrtMCl->Fill(kk, sum2DClusterEnergy_vsDR_wrtMCl[kk] / it_caloPart->energy());
      tp_2DClSumE_vsRadius_wrtGen->Fill(kk, sum2DClusterEnergy_vsRadius_wrtGen[kk] / it_caloPart->energy());
      tp_2DClSumE_vsRadius_wrtMCl->Fill(kk, sum2DClusterEnergy_vsRadius_wrtMCl[kk] / it_caloPart->energy());
    }

    h_nMCPerEvt->Fill(nmclusCount);

    //start fixme => need to include in pt-gun matched

    /*
    float allRHSeedEnergy = 0;
    //    int allRHSeedIndex = 0;
    unsigned int allRHSeedLayer = 0;
    for(std::map<DetId,const HGCRecHit*>::const_iterator jk = hitmap.begin(); jk != hitmap.end(); ++jk){

      const HGCRecHit *hit = jk->second;
      const HGCalDetId detid = jk->first;
      
      if(hit->energy() > allRHSeedEnergy){
	allRHSeedEnergy = hit->energy();
	allRHSeedLayer = recHitTools.getLayerWithOffset(detid);
      }
      // if(hit->time()-1 > 0){
      // 	float hitEta = recHitTools.getEta(detid);
      // 	float hitPhi = recHitTools.getPhi(detid);
      // 	float dRwrtGen = reco::deltaR(hitEta, hitPhi, it_caloPart->eta(), it_caloPart->phi());

      // 	int timesOfRadius = int(dRwrtGen / recHitTools.getRadiusToSide(detid) );
      // 	if(timesOfRadius < 10) h_allhitTimeinRing_wrtGen[timesOfRadius]->Fill(hit->time()-1);
      // 	else h_allhitTimeinRing_wrtGen[10]->Fill(hit->time()-1);

      // 	h_allRH_TimesOfRadius->Fill(dRwrtGen / recHitTools.getRadiusToSide(detid));

      // 	h2_allRH_TimesVsDRaxis->Fill(dRwrtGen, hit->time()-1);
      // 	tp_allRH_TimesVsDRaxis->Fill(dRwrtGen, hit->time()-1);

      // }

    }
    for(std::map<DetId,const HGCRecHit*>::const_iterator jk = hitmap.begin(); jk != hitmap.end(); ++jk){

	const HGCRecHit *hit = jk->second;
	const HGCalDetId detid = jk->first;

	recHitTools.getLayerWithOffset(detid);
	if(recHitTools.getLayerWithOffset(detid) == allRHSeedLayer){

	if(hit->time()-1 > 0){
	  float hitEta = recHitTools.getEta(detid);
	  float hitPhi = recHitTools.getPhi(detid);
	  float dRwrtGen = reco::deltaR(hitEta, hitPhi, it_caloPart->eta(), it_caloPart->phi());

	  int timesOfRadius = int(dRwrtGen / recHitTools.getRadiusToSide(detid) );
	  if(timesOfRadius < 10) h_allhitTimeinRing_wrtGen_fixL[timesOfRadius]->Fill(hit->time()-1);
	  else h_allhitTimeinRing_wrtGen_fixL[10]->Fill(hit->time()-1);
	}
      }
    }
    */
    //end fixme 
  


    //shower axis by recHits
    float axisX = 0;
    float axisY = 0;
    float axisZ = 0;
    float sumEnergyToNorm = 0;
    GlobalPoint showerAxis;

    float sumMtcRHEnergy = 0;
    float sumMtcRHEnergy_vsDR_wrtGen[5];
    //    float sumMtcRHEnergy_vsDR_wrtMCl[5];
    float sumMtcRHEnergy_vsRadius_wrtGen[5];
    //    float sumMtcRHEnergy_vsRadius_wrtMCl[5];
    for(int ij=0; ij<5; ++ij){
      sumMtcRHEnergy_vsDR_wrtGen[ij] = 0;
      sumMtcRHEnergy_vsRadius_wrtGen[ij] = 0;
      //      sumMtcRHEnergy_vsDR_wrtMCl[ij] = 0;
      //      sumMtcRHEnergy_vsRadius_wrtMCl[ij] = 0;
    }
    //loop on rechit - matched to gen
    const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
    for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin(); it_sc != simClusterRefVector.end(); ++it_sc) {
      const SimCluster simCluster = (*(*it_sc));
      const std::vector<std::pair<uint32_t,float> > hits_and_fractions = simCluster.hits_and_fractions();

      
      for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
	//	unsigned int hitlayer = recHitTools.getLayerWithOffset(it_haf->first);
	DetId hitid = (it_haf->first);

	//	bool found = false;
	float rhEnergy = 0;
	//	float rhTime = -1;

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

      unsigned int allRHSeedLayer = 0; 
      float allRHSeedEnergy = 0;   

      float timePerEtaRadius[8][3];
      int totRHPerEtaRadius[8][3];
      int totRHPerEtaRadius_allTime[8][3];
      float timePerEtaRadius_fixL[8][2];
      int totRHPerEtaRadius_fixL[8][2];
      int totRHPerEtaRadius_allTime_fixL[8][2];
      for(int iet=0; iet<8; ++iet){
	for(int irad=0; irad<3; ++irad){
	  timePerEtaRadius[iet][irad] = 0.;
	  totRHPerEtaRadius[iet][irad] = 0;
	  totRHPerEtaRadius_allTime[iet][irad] = 0;
	  if(irad > 1) continue;
	  timePerEtaRadius_fixL[iet][irad] = 0.;
	  totRHPerEtaRadius_fixL[iet][irad] = 0;
	  totRHPerEtaRadius_allTime_fixL[iet][irad] = 0;
	}
      }
      //loop over hits
      for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
	DetId hitid = (it_haf->first);

	bool found = false;
	float rhEnergy = 0;
	float rhTime = -1;

	std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmap.find(hitid);
	if(trovatore == hitmap.end()){
	  continue;
	}
	else{
	  const HGCRecHit *hit = hitmap[hitid];
	  rhEnergy = hit->energy();
	  rhTime = hit->time() - 1;
	  found = true;
	  // if(recHitTools.getEta(hitid)*it_caloPart->eta() < 0){
	  //   std::cout << " >>> MAJOR PROBLEM!!! " << std::endl;

	    if(it_caloPart->eta() > 0){
	      h_dEta_pos->Fill(recHitTools.getEta(hitid) - it_caloPart->eta());
	      h_dPhi_pos->Fill(reco::deltaPhi(recHitTools.getPhi(hitid), it_caloPart->phi()) );
	    }
	    if(it_caloPart->eta() < 0){
	      h_dEta_neg->Fill(recHitTools.getEta(hitid) - it_caloPart->eta());
	      h_dPhi_neg->Fill(reco::deltaPhi(recHitTools.getPhi(hitid), it_caloPart->phi()) );
	    }

	  //   return;
	  // }

	  if(hit->energy()*it_haf->second > allRHSeedEnergy){
	    allRHSeedEnergy = hit->energy()*it_haf->second;
	    allRHSeedLayer = recHitTools.getLayerWithOffset(hitid);
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

	  h_rHGen_dR->Fill(Radius_2DClGen);
	  h_rHGen_dDist->Fill(dR_2DClGen);
	  h2_rHGen_dR_vsEta->Fill(std::abs(etaGen), Radius_2DClGen);
	  h2_rHGen_dDist_vsEta->Fill(std::abs(etaGen), dR_2DClGen);
	  tp_rHGen_dDist_vsEta->Fill(std::abs(etaGen), dR_2DClGen);
	  h2_Energia_vsEta->Fill(std::abs(etaGen), it_caloPart->energy());
	  h2_Pt_vsEta->Fill(std::abs(etaGen), it_caloPart->pt());
	  h2_rHGen_dR_vsEnergia->Fill(it_caloPart->energy(), Radius_2DClGen);
	  h2_rHGen_dR_vsPt->Fill(it_caloPart->pt(), Radius_2DClGen);

	  tp2_Energy_vs_Eta_dR->Fill(std::abs(etaGen), it_caloPart->energy(), Radius_2DClGen);
	  tp2_Energy_vs_Eta_dDist->Fill(std::abs(etaGen), it_caloPart->energy(), dR_2DClGen);

	  h2_rHGen_dRcosh_vsEta->Fill(std::abs(etaGen), Radius_2DClGen/cosh(etaGen));
	  h2_rHGen_dRcoshP_vsEta->Fill(std::abs(etaGen), Radius_2DClGen*cosh(etaGen));

	  //FIXME
	  int etaBin = int((std::abs(etaGen) - 1.5) / 0.2);
	  int iRadBin = -1;
	  for(int ir=0; ir<3; ++ir){
	    if(Radius_2DClGen < radiusEtaRad[etaBin][ir]) {
	      iRadBin = ir;
	      break;
	    }
	  }

	  if(iRadBin != -1) {
	    for(int ir=iRadBin; ir<=2; ++ir) totRHPerEtaRadius_allTime[etaBin][ir] += 1;
	  }

	  if(rhTime > -1.){
	    int timesOfRadius = int(Radius_2DClGen );

	    if(timesOfRadius < 5) h_allhitTimeinRing_wrtGen[timesOfRadius]->Fill(rhTime);
	    else h_allhitTimeinRing_wrtGen[5]->Fill(rhTime);

	    h_rhGen_Radius_Eta[etaBin]->Fill(Radius_2DClGen);
	    h_rhGen_RadiusDdeta_Eta[etaBin]->Fill(Radius_2DClGen / std::abs(rhEta - etaGen));
	    h_rhGen_RadiusPdeta_Eta[etaBin]->Fill(Radius_2DClGen * std::abs(rhEta - etaGen));

	    if(iRadBin != -1) {
	      for(int ir=iRadBin; ir<=2; ++ir){
		timePerEtaRadius[etaBin][ir] += rhTime; 
		totRHPerEtaRadius[etaBin][ir] += 1;
	      }
	    }

	    dPhivsDEta_Time->Fill(sqrt(pow(rhX - axX, 2)), sqrt(pow(rhY - axY, 2)), rhTime);
	    dRvsDZ_Time->Fill(sqrt(pow(rhX - axX, 2) + pow(rhY - axY, 2)),  std::abs(rhZ - axZ), rhTime);

	    h_allRH_TimesOfRadius->Fill(Radius_2DClGen);
	    h2_allRH_TimesVsDRaxis->Fill(dR_2DClGen, rhTime);
	    tp_allRH_TimesVsDRaxis->Fill(dR_2DClGen, rhTime);
	    h2_allRH_TimesVsDEPaxis->Fill(sqrt(pow(rhX - axX, 2) + pow(rhY - axY, 2)), rhTime);
	    tp_allRH_TimesVsDEPaxis->Fill(sqrt(pow(rhX - axX, 2) + pow(rhY - axY, 2)), rhTime);
	    h2_allRH_TimesVsRadiusaxis->Fill(Radius_2DClGen, rhTime);
	    tp_allRH_TimesVsRadiusaxis->Fill(Radius_2DClGen, rhTime);
	  }
	  
	}
      }// first loop over rechits


      //second to select seed layer
      for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
        DetId hitid = (it_haf->first);

        bool found = false;
        float rhTime = -1;

	std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmap.find(hitid);
        if(trovatore == hitmap.end()){
          continue;
        }
        else{
	  const HGCRecHit *hit = hitmap[hitid];
          rhTime = hit->time() - 1;
          found = true;
        }
        if(found && recHitTools.getLayerWithOffset(hitid) == allRHSeedLayer){
	  float rhEta = recHitTools.getEta(hitid);
	  float rhPhi = recHitTools.getPhi(hitid);
	  h2_dPhivsdEta_rhGen_fixL->Fill(reco::deltaPhi(rhPhi, it_caloPart->phi()), rhEta - it_caloPart->eta());

	  float dRwrtGen1 = reco::deltaR(etaGen, rhEta, phiGen, rhPhi);

	  float rhX = recHitTools.getPosition(hitid).x();
	  float rhY = recHitTools.getPosition(hitid).y();
	  float rhL = recHitTools.getLayerWithOffset(hitid);
	  
	  float dRwrtGen = utilsMet.dsGenRecHit(etaGen, phiGen, rhL, rhX, rhY);
	  
	  h_rHGen_dR_fixL->Fill(dRwrtGen);

	  int etaBin = int((std::abs(etaGen) - 1.5) / 0.2);
	  int iRadBin = -1;
	  for(int ir=0; ir<2; ++ir){
            if(dRwrtGen < radiusEtaRad[etaBin][ir]) {
              iRadBin = ir;
              break;
            }
          }

          if(iRadBin != -1){
	    for(int ir=iRadBin; ir<=1; ++ir) totRHPerEtaRadius_allTime_fixL[etaBin][ir] += 1;
	  }


	  if(rhTime > -1.){
	    //	    int timesOfdR = int(dRwrtGen1 / recHitTools.getRadiusToSide(hitid) );
	    //  int timesOfdR = int(dRwrtGen1);

	    h2_allRH_TimesVsDRaxis_fixL->Fill(dRwrtGen1, rhTime);
	    tp_allRH_TimesVsDRaxis_fixL->Fill(dRwrtGen1, rhTime);

	    int timesOfRadius = int(dRwrtGen);

	    if(iRadBin != -1) {
	      for(int ir=0; ir<=iRadBin; ++ir){
		timePerEtaRadius_fixL[etaBin][ir] += rhTime;
		totRHPerEtaRadius_fixL[etaBin][ir] += 1;
	      }
	    }


	    h_allRH_TimesOfRadius_fixL->Fill(dRwrtGen);
	    h2_allRH_TimesVsRadiusaxis_fixL->Fill(dRwrtGen, rhTime);
            tp_allRH_TimesVsRadiusaxis_fixL->Fill(dRwrtGen, rhTime);

	    if(timesOfRadius < 5) h_allhitTimeinRing_wrtGen_fixL[timesOfRadius]->Fill(rhTime);
	    else h_allhitTimeinRing_wrtGen_fixL[5]->Fill(rhTime);
	  }
	}
      }


      for(int iet=0; iet<8; ++iet){
	for(int irad=0; irad<3; ++irad){
	  totEvtsEtaRadius[iet][irad] += 1;

	  if(totRHPerEtaRadius[iet][irad] < 3) continue;

	  totEvtsEtaRadius_withTime[iet][irad] += 1;
	  // std::cout << " totRHPerEtaRadius[iet][irad] = " << totRHPerEtaRadius[iet][irad] << std::endl;
	  // std::cout << " totRHPerEtaRadius_allTime[iet][irad] = " << totRHPerEtaRadius_allTime[iet][irad] << std::endl;
	  if(totRHPerEtaRadius[iet][irad] != 0)
	    hAverageTime_Eta_dRadius[iet][irad]->Fill(1.*timePerEtaRadius[iet][irad]/totRHPerEtaRadius[iet][irad]);
	  if(totRHPerEtaRadius_allTime[iet][irad] != 0){
	    //	    std::cout << " frazione = " << 1.*totRHPerEtaRadius[iet][irad]/totRHPerEtaRadius_allTime[iet][irad] << std::endl;
	    hFractionHitsWithTime_Eta_dRadius[iet][irad]->Fill(1.*totRHPerEtaRadius[iet][irad]/totRHPerEtaRadius_allTime[iet][irad]);	  
	    tpFractionHitsWithTime_Eta_dRadius_vsPt[iet][irad]->Fill(it_caloPart->energy(), 1.*totRHPerEtaRadius[iet][irad]/totRHPerEtaRadius_allTime[iet][irad]);
	    hTotHits_Eta_dRadius[iet][irad]->Fill(totRHPerEtaRadius_allTime[iet][irad]);
	    hTotHitsWithTime_Eta_dRadius[iet][irad]->Fill(totRHPerEtaRadius[iet][irad]);
	  }
	  if(irad > 1) continue;
	  if(totRHPerEtaRadius_fixL[iet][irad] != 0)
	    hAverageTime_Eta_dRadius_fixL[iet][irad]->Fill(1.*timePerEtaRadius_fixL[iet][irad]/totRHPerEtaRadius_fixL[iet][irad]);
	  if(totRHPerEtaRadius_allTime_fixL[iet][irad] != 0){
	    hFractionHitsWithTime_Eta_dRadius_fixL[iet][irad]->Fill(1.*totRHPerEtaRadius_fixL[iet][irad]/totRHPerEtaRadius_allTime_fixL[iet][irad]);	  
	    hTotHits_Eta_dRadius_fixL[iet][irad]->Fill(totRHPerEtaRadius_allTime_fixL[iet][irad]);
	    hTotHitsWithTime_Eta_dRadius_fixL[iet][irad]->Fill(totRHPerEtaRadius_fixL[iet][irad]);
	  }
	}
      }

    }
    h_mtcRHSumE -> Fill(sumMtcRHEnergy/it_caloPart->energy());
    for(int ij=0; ij<5; ++ij){
      tp_mtcRHSumE_vsDR_wrtGen->Fill(ij, sumMtcRHEnergy_vsDR_wrtGen[ij]/it_caloPart->energy());
      //tp_mtcRHSumE_vsDR_wrtMCl->Fill(ij, sumMtcRHEnergy_vsDR_wrtGen[ij]/it_caloPart->energy());
      tp_mtcRHSumE_vsRadius_wrtGen->Fill(ij, sumMtcRHEnergy_vsRadius_wrtGen[ij]/it_caloPart->energy());
      //tp_mtcRHSumE_vsRadius_wrtMCl;
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

  for(int iet=0; iet<8; ++iet){
    for(int irad=0; irad<3; ++irad){
    if(totEvtsEtaRadius[iet][irad] != 0)  hFractionEvents_HitsWithTime_Eta_dRadius[iet][irad]->Fill(1.*totEvtsEtaRadius_withTime[iet][irad]/totEvtsEtaRadius[iet][irad]);
    }
  }
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
