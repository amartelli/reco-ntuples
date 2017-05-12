// user include files
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include "TRandom3.h"

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
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "RecoNtuples/HGCalAnalysis/plugins/RecHiTimeEstimator.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

RecHiTimeEstimator::RecHiTimeEstimator(const edm::ParameterSet& ps){

  //  std::cout << " >>> constructor " << std::endl;
  //section type
  keV2fC[0] =  ps.getParameter<double>("HGCEE_keV2fC");
  keV2fC[1] =  ps.getParameter<double>("HGCHEF_keV2fC");
  keV2MIP = ps.getParameter<double>("HGCHB_keV2MIP");

  noisefC[0] = (ps.getParameter<std::vector<double> >("HGCEE_noisefC")).at(0);
  noisefC[1] = (ps.getParameter<std::vector<double> >("HGCEF_noisefC")).at(1);
  noiseMIP = ps.getParameter<double>("HGCBH_noiseMIP");

  //cell type
  fCPerMIP[0] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(0);
  fCPerMIP[1] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(1);
  fCPerMIP[2] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(2);

  const auto& rcorr = ps.getParameter<std::vector<double> >("thicknessCorrection");
  scaleCorrection.clear();
  for( auto corr : rcorr ) {
    scaleCorrection.push_back(1.0/corr);
  }

  const auto& dweights = ps.getParameter<std::vector<double> >("dEdXweights");
  for( auto weight : dweights ) {
    weights.push_back(weight);
  }

  keV2GeV = 1e-6;
  keV2MeV = 1e-3;

  //100um
  paramA[0] = 1.00;
  parErrA[0] = 0.01;
  paramC[0] = 0.009;
  parErrC[0] = 0.001;

  //200um
  paramA[1] = 1.06;
  parErrA[1] = 0.02;
  paramC[1] = 0.008;
  parErrC[1] = 0.001;

  //300um
  paramA[2] = 1.11;
  parErrA[2] = 0.02;
  paramC[2] = 0.010;
  parErrC[2] = 0.001;

  timeResolution = new TF1("timeSi100", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
}


void RecHiTimeEstimator::setEventSetup(const edm::EventSetup& es){
  recHitTools.getEventSetup(es);
}


double RecHiTimeEstimator::getTimeHit(int thick, double SoverN){
  timeResolution->SetParameters(paramA[thick], paramC[thick]);

  //resolution from TB results with floor of 20ps at high S/N
  double sigma = 0.2;
  if(SoverN > 1) sigma = timeResolution->Eval(SoverN);
  if(SoverN > 20) sigma = 0.02;

  TRandom3* rand = new TRandom3(); 
  double smearing = rand->Gaus(0., sigma);
  return smearing;
}



double RecHiTimeEstimator::getTimeHitFixThr(){
  //flat resolution at 50ps
  double sigma = 0.05;

  TRandom3* rand = new TRandom3(); 
  double smearing = rand->Gaus(0., sigma);
  return smearing;
}


void RecHiTimeEstimator::correctTime(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits){

  for(HGCRecHitCollection::const_iterator it_hit = rechits.begin(); it_hit < rechits.end(); ++it_hit) {
    const DetId detid = it_hit->detid();

    int thick = (detid.det() != DetId::Forward) ? -1 : recHitTools.getSiThickness(detid) / 100. - 1.;

    int sectionType = -1;
    if(detid.subdetId() == HGCEE) sectionType = 0;
    else if(detid.subdetId() == HGCHEF) sectionType = 1;
    else if(detid.subdetId() == HGCHEB) sectionType = 2;

    HGCRecHit myrechit(*it_hit);
    float energy = it_hit->energy();
    float time = it_hit->time();

    if(sectionType == -1 || thick == -1){
      myrechit.setTime(-1.);
      Newrechits->push_back(myrechit);
      continue;
    }

    unsigned int layer = recHitTools.getLayerWithOffset(detid);

    double sigmaNoiseMIP = noiseMIP;
    if(sectionType != 2) sigmaNoiseMIP = noisefC[sectionType]/fCPerMIP[thick];

    int energyMIP = 0.;
    if(sectionType == 2) energyMIP = energy/keV2GeV * keV2MIP;
    else if(sectionType == 0 || sectionType == 1) energyMIP = energy/scaleCorrection.at(thick)/keV2GeV / (weights.at(layer)/keV2MeV);

    if(energyMIP > 3.){
      float SoverN = energyMIP / sigmaNoiseMIP;
      double smearedTime = getTimeHit(thick, SoverN);
      myrechit.setTime((time - 1.) * (1+smearedTime) + 1.);
    }
    else myrechit.setTime(-1.); 

    Newrechits->push_back(myrechit);
  }
  return;
}


void RecHiTimeEstimator::correctTimeFixThr(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits){

  for(HGCRecHitCollection::const_iterator it_hit = rechits.begin(); it_hit < rechits.end(); ++it_hit) {
    const DetId detid = it_hit->detid();
    int thick = (detid.det() != DetId::Forward) ? -1 : recHitTools.getSiThickness(detid) / 100. - 1.;

    int sectionType = -1;
    if(detid.subdetId() == HGCEE) sectionType = 0;
    else if(detid.subdetId() == HGCHEF) sectionType = 1;
    else if(detid.subdetId() == HGCHEB) sectionType = 2;

    HGCRecHit myrechit(*it_hit);
    float energy = it_hit->energy();
    float time = it_hit->time();

    if(sectionType == -1 || thick == -1){
      myrechit.setTime(-1.);
      Newrechits->push_back(myrechit);

      continue;
    }

    unsigned int layer = recHitTools.getLayerWithOffset(detid);

    int energyMIP = 0.;
    if(sectionType == 2) energyMIP = energy/keV2GeV * keV2MIP;
    else if(sectionType == 1 || sectionType == 0) energyMIP = energy/scaleCorrection.at(thick)/keV2GeV / (weights.at(layer)/keV2MeV);

    float energyCharge = 0.;
    // from SimCalorimetry/HGCalSimProducers/src/HGCHEbackDigitizer.cc L 58
    if(sectionType == 2) energyCharge = energyMIP * 1.; 
    else if(sectionType == 1 || sectionType == 0) energyCharge = energyMIP * fCPerMIP[thick];


    if(energyCharge > 60){
      double smearedTime = getTimeHitFixThr();
      myrechit.setTime((time-1) * (1+smearedTime) + 1.);
    }
    else myrechit.setTime(-1.); 

    Newrechits->push_back(myrechit);
  }
  return;
}


