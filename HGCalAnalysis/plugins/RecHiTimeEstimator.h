#ifndef RecHiTimeEstimator_h
#define RecHiTimeEstimator_h

// user include files
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TRandom3.h"
#include "TProfile.h"
using namespace std;




//legend
/*
CS   = cell standard (1cm2)
CH   = cell half (0.5cm2) => effect ~half timing resolution
F20  = Floor 20ps time at high S/N
F30  = Floor 30ps time at high S/N
LB   = life at beginning for charge collection efficiency
LE   = life end for charge collection efficiency (70% 50% 50% for 300, 200, 100)
*/


class RecHiTimeEstimator
{
  
public:
  explicit RecHiTimeEstimator(const edm::ParameterSet& ps);
  ~RecHiTimeEstimator();

  double getTimeHit(int thick, double SoverN);
  double getTimeHitFixThr();
  void correctTime(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits);
  void correctTimeFixThr(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits);

  void setOptions(int cellType, float floor, int liveAge);

  /* void timeCSF20LB(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCSF30LB(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCHF20LB(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCHF30LB(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCSF20LE(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCSF30LE(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCHF20LE(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCHF30LE(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */

  void setEventSetup(const edm::EventSetup& es);

private:

  hgcal::RecHitTools recHitTools;

  TF1* timeResolution;

  float paramA[3];
  float parErrA[3];
  float paramC[3];
  float parErrC[3];

  float floorValue;

  float chargeCollEff[3];

  float cellSize[3];

  std::vector<float> scaleCorrection;
  std::vector<float> weights;

  float keV2GeV;
  float keV2MeV;

  double keV2fC[2];
  double keV2MIP;

  double noisefC[2];
  double noiseMIP;

  //for cell type
  double fCPerMIP[3];
};

#endif
