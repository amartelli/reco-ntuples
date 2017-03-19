// #include <iostream>
// #include <cmath>
// #include <string>
// #include <vector>
// #include "TH1.h"
// #include "TCanvas.h"

#include "RecoNtuples/HGCalAnalysis/interface/UsefulClasses.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// using std::string;
// using std::vector;
// using std::cout;
// using std::endl;
// using std::map;


UsefulClasses::UsefulClasses(float genEta, float genPhi):_genEta(genEta), _genPhi(genPhi)
{}

UsefulClasses::~UsefulClasses(){
}
/*
 float deltaPhi( float phi1, float phi2) {   
  float dPhi = phi1 - phi2;                                                                                                                                       
  float pi = 3.14159265;                                                                                                                               
  if     ( dPhi <=-pi) dPhi += 2.0*pi;                                                                                                                
  else if( dPhi >  pi) dPhi -= 2.0*pi; 
  return dPhi;  
}
*/
 float UsefulClasses::deltaEta(float eta1, float eta2){
    float dEta = eta1 - eta2;
    return dEta;
}

 float UsefulClasses::deltaR(float eta1, float eta2, float phi1, float phi2) {
    float dEta = deltaEta(eta1, eta2);
    float dPhi = reco::deltaPhi(phi1, phi2);
    return sqrt( dEta*dEta + dPhi*dPhi );
}

 float UsefulClasses::deltaX(float x1, float x2, float y1, float y2) {
    float dX = x1 - x2;
    float dY = y1 - y2;
    return sqrt( dX*dX + dY*dY );
}

// need correct sign for z
 float UsefulClasses::etaPhiZtoX( float eta, float phi, float z )
{
  float t = exp( -1 * eta );
  float x = z * 2 * t * cos(phi) / ( 1 - t*t );
  return x;
}

// need correct sign for z
 float UsefulClasses::etaPhiZtoY( float eta, float phi, float z )
{
  float t = exp( -1 * eta );
  float y = z * 2 * t * sin(phi) / ( 1 - t*t );
  return y;
}

// go from HGC layer to z in cm
 float UsefulClasses::layerToZ( int layer, float eta ) 
{
    map<int, float> lToZ;
    lToZ[0]  = 320.75;
    lToZ[1]  = 321.50;
    lToZ[2]  = 322.73;
    lToZ[3]  = 323.48;
    lToZ[4]  = 324.71;
    lToZ[5]  = 325.46;
    lToZ[6]  = 326.69;
    lToZ[7]  = 327.44;
    lToZ[8]  = 328.67;
    lToZ[9]  = 329.42; //first set
    lToZ[10] = 330.73;
    lToZ[11] = 331.60;
    lToZ[12] = 332.91;
    lToZ[13] = 333.78;
    lToZ[14] = 335.09;
    lToZ[15] = 335.96;
    lToZ[16] = 337.27;
    lToZ[17] = 338.14;
    lToZ[18] = 339.45;
    lToZ[19] = 340.32; //second set
    lToZ[20] = 341.77;
    lToZ[21] = 342.84;
    lToZ[22] = 344.29;
    lToZ[23] = 345.36;
    lToZ[24] = 346.81;
    lToZ[25] = 347.88;
    lToZ[26] = 349.33;
    lToZ[27] = 350.40; //third set
    lToZ[28] = 356.33;
    lToZ[29] = 361.01;
    lToZ[30] = 365.69;
    lToZ[31] = 370.37;
    lToZ[32] = 375.05;
    lToZ[33] = 379.73;
    lToZ[34] = 384.41;
    lToZ[35] = 389.09;
    lToZ[36] = 393.77;
    lToZ[37] = 398.45;
    lToZ[38] = 403.13;
    lToZ[39] = 407.81; //fourth set
    
    float z = lToZ[ layer ];
    if( eta < 0 ) z *= -1.;
    return z;
}


int UsefulClasses::Ztolayer(float z, float eta ) 
{
    map<int, float> lToZ;
    lToZ[0]  = 320.75;
    lToZ[1]  = 321.50;
    lToZ[2]  = 322.73;
    lToZ[3]  = 323.48;
    lToZ[4]  = 324.71;
    lToZ[5]  = 325.46;
    lToZ[6]  = 326.69;
    lToZ[7]  = 327.44;
    lToZ[8]  = 328.67;
    lToZ[9]  = 329.42; //first set
    lToZ[10] = 330.73;
    lToZ[11] = 331.60;
    lToZ[12] = 332.91;
    lToZ[13] = 333.78;
    lToZ[14] = 335.09;
    lToZ[15] = 335.96;
    lToZ[16] = 337.27;
    lToZ[17] = 338.14;
    lToZ[18] = 339.45;
    lToZ[19] = 340.32; //second set
    lToZ[20] = 341.77;
    lToZ[21] = 342.84;
    lToZ[22] = 344.29;
    lToZ[23] = 345.36;
    lToZ[24] = 346.81;
    lToZ[25] = 347.88;
    lToZ[26] = 349.33;
    lToZ[27] = 350.40; //third set
    lToZ[28] = 356.33;
    lToZ[29] = 361.01;
    lToZ[30] = 365.69;
    lToZ[31] = 370.37;
    lToZ[32] = 375.05;
    lToZ[33] = 379.73;
    lToZ[34] = 384.41;
    lToZ[35] = 389.09;
    lToZ[36] = 393.77;
    lToZ[37] = 398.45;
    lToZ[38] = 403.13;
    lToZ[39] = 407.81; //fourth set
    
    // float zedValue = (z > 0 ? z : -1.*z);
    // int layerOK = find(lToZ.begin(), lToZ.end(), zedValue)->first;
 
   unsigned int nstep = 0;
    std::map<int,float>::iterator it2 = lToZ.begin();
    ++it2;

    int layer = 0;
    for(std::map<int,float>::iterator it=lToZ.begin(); it!=lToZ.end(); ++it){
      if(nstep >= lToZ.size()-1) continue;
      //      std::cout << " it->second = " << it->second << " (it2)->second = " << (it2)->second << std::endl; 
      if(eta < 0 && z * -1. >= it->second && z * -1. < (it2)->second) layer = it->first;
      else{
	if(z >= it->second && z < (it2)->second) layer = it->first;
      }
      ++nstep;
      ++it2;
    }

    //    std::cout << " func Z = " << z << " layer = " << layer << std::endl;
    //    if(layerOK != layer) std::cout << " BIG PROOOOOBLEM " << std::endl;
    return layer;
}

 float UsefulClasses::dsGenRecHit(float genEta, float genPhi, int recHitLayer, float recHitX, float recHitY )
{

    float genZ = layerToZ( recHitLayer, genEta );
    float genX = etaPhiZtoX( genEta, genPhi, genZ );
    float genY = etaPhiZtoY( genEta, genPhi, genZ );
    float ds   = deltaX( genX, recHitX , genY, recHitY );

    return ds; 
}

 float UsefulClasses::dsGenRecoObj(float genEta, float genPhi, float recoObjZ, float recoObjX, float recoObjY )
{

    float genZ = recoObjZ;
    float genX = etaPhiZtoX( genEta, genPhi, genZ );
    float genY = etaPhiZtoY( genEta, genPhi, genZ );
    float ds   = deltaX( genX, recoObjX , genY, recoObjY );

    return ds; 
}


 vector<TH1*> UsefulClasses::initHists( string name )
{
  vector<TH1*> v;
  string tempName  = "h" + name + "Eta";
  string tempTitle = name + " eta";
  v.push_back(     new TH1F( tempName.c_str(), tempTitle.c_str(), 50, -5., 5. ) );
  tempName     = "h" + name + "Phi";
  tempTitle    = name + " phi";
  v.push_back( new TH1F( tempName.c_str(), tempTitle.c_str(), 50, -3.14159, 3.14159 ) );
  tempName     = "h" + name + "Energy";
  tempTitle    = name + " energy";
  v.push_back( new TH1F( tempName.c_str(), tempTitle.c_str(), 50, 0., 500. ) );
  tempName     = "h" + name + "Pt";
  tempTitle    = name + " pt";
  v.push_back( new TH1F( tempName.c_str(), tempTitle.c_str(), 50, 0., 50. ) );
  tempName     = "hMatched" + name + "dR";
  tempTitle    = "deltaR between gen particle and " + name;
  v.push_back( new TH1F( tempName.c_str(), tempTitle.c_str(), 50, 0., 1. ) );
  tempName     = "hMatched" + name + "EFrac";
  tempTitle    = "fraction of gen energy within " + name + "s";
  v.push_back( new TH1F( tempName.c_str(), tempTitle.c_str(), 50, 0., 1.5 ) );
  tempName     = "hMatched" + name + "PtFrac";
  tempTitle    = "fraction of gen pt within " + name + "s";
  v.push_back( new TH1F( tempName.c_str(), tempTitle.c_str(), 50, 0., 1.5 ) );
  tempName     = "hNum" + name + "NearGen";
  tempTitle    = "number of " + name + "with dR < 0.1 to gen";
  v.push_back( new TH1F( tempName.c_str(), tempTitle.c_str(), 50, 0, 50 ) );
  return v;
}


 void UsefulClasses::writeHists( vector<TH1*> hists )
{
  for( uint i=0; i<hists.size(); i++ ) hists[i]->Write();
}

 void UsefulClasses::drawAndSaveHists( vector<TH1*> hists, TCanvas *c, bool doLogs)
{
  for( uint i=0; i<hists.size(); i++ )  
  {
    c->SetLogy( 0 );
    hists[i]->Draw("hist");
    hists[i]->GetXaxis()->SetTitleSize(0.05);
    hists[i]->GetXaxis()->SetTitleOffset(0.8);
    string histName = string( hists[i]->GetName() );
    string saveName = "HGCPlots/" + histName + ".pdf";
    c->Print( saveName.c_str() );
    saveName = "HGCPlots/" + histName + ".png";
    c->Print( saveName.c_str() );
    c->SetLogy( 1 );
    if( doLogs ) {
      saveName = "HGCPlots/" + histName + "Log.pdf";
      c->Print( saveName.c_str() );
      saveName = "HGCPlots/" + histName + "Log.png";
      c->Print( saveName.c_str() );
    }
  }
}


//DEFINE_FWK_MODULE(UsefulClasses);
