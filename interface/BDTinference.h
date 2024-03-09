#ifndef bdtinference_h
#define bdtinference_h

// Standard libraries
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <xgboost/c_api.h> 

#include "DataFormats/Math/interface/deltaR.h"

float sigmoid(float x){
  return (1./(1 + std::exp(-1.*x)));
}

struct jet_legacy_t {
    float pt = -999.;
    float eta = -999.;
    float phi = -999.;
    float mass = -999.;
    float chEmEF = -999.;
    float chHEF = -999.;
    float neEmEF = -999.;
    float neHEF = -999.;
    float muEF = -999.;
    float muonSubtrFactor = -999.;
    float chFPV0EF = -999.;
    float nMuons = -999.;
    float nElectrons = -999.;
    float nConstituents = -999.;
    float btagDeepB = -999.;
    float btagDeepC = -999.;
    float qgl = -999.;
    float puIdDisc = -999.;
    float muonIdx1 = -999.;
    float muonIdx2 = -999.;
    //jet_legacy_t():
        //pt(-999), eta(-999), phi(-999), mass(-999), chEmEF(-999),
        //chHEF(-999), neEmEF(-999), neHEF(-999), muEF(-999),
        //muonSubtrFactor(-999), chFPV0EF(-999), nMuons(-999), nElectrons(-999),
        //nConstituents(-999), btagDeepB(-999), btagDeepC(-999), qgl(-999),
        //puIdDisc(-999), muonIdx1(-999), muonIdx2(-999) {}
};

struct jet_ul_t {
    float pt = -999.;
    float eta = -999.;
    float phi = -999.;
    float mass = -999.;
    float chEmEF = -999.;
    float chHEF = -999.;
    float neEmEF = -999.;
    float neHEF = -999.;
    float muEF = -999.;
    float muonSubtrFactor = -999.;
    float chFPV0EF = -999.;
    float nMuons = -999.;
    float nElectrons = -999.;
    float nConstituents = -999.;
    float btagDeepB = -999.;
    // float btagDeepC = -999.;
    float qgl = -999.;
    float puIdDisc = -999.;
    float muonIdx1 = -999.;
    float muonIdx2 = -999.;
};

struct jet_scouting_t {
    float pt = -999.;
    float eta = -999.;
    float phi = -999.;
    float mass = -999.;
};

bool jet_legacy_sort (const jet_legacy_t& jA, const jet_legacy_t& jB)
{
  return (jA.pt > jB.pt);
}

bool jet_ul_sort (const jet_ul_t& jA, const jet_ul_t& jB)
{
  return (jA.pt > jB.pt);
}

bool jet_scouting_sort (const jet_scouting_t& jA, const jet_scouting_t& jB)
{
  return (jA.pt > jB.pt);
}

struct muon_legacy_t {
    float eta = -999.;
    float phi = -999.;
    float pt = -999.;
    float ptErr = -999.;
    float dxy = -999.;
    float dxyErr = -999.;
    float dz = -999.;
    float dzErr = -999.;
    float ip3d = -999.;
    float sip3d = -999.;
    float charge = -999.;
    float tightId = -999.;
    float softMva = -999.;
    float pfRelIso03_all = -999.;
    float miniPFRelIso_all = -999.;
    float jetIdx = -999.;
    //muon_legacy_t():
        //eta(-999), phi(-999), pt(-999), ptErr(-999), dxy(-999), dxyErr(-999), dz(-999), dzErr(-999),
        //ip3d(-999), sip3d(-999), charge(-999), tightId(-999), softMva(-999),
        //pfRelIso03_all(-999), miniPFRelIso_all(-999), jetIdx(-999) {}
};

struct muon_scouting_t {
    float eta = -999.;
    float phi = -999.;
    float pt = -999.;
    float normalizedChi2 = -999.;
    float ecalIso = -999.;
    float hcalIso = -999.;
    float trackIso = -999.;
    float dxy = -999.;
    float dxyErr = -999.;
    float dz = -999.;
    float dzErr = -999.;
    float charge = -999.;
    float isTracker = -999.;
    float isGlobal = -999.;
    float isPFmatched = -999.;
    float isStandalone = -999.;
    float nStations = -999.;
    float nValidPixelHits = -999.;
    float nValidStripHits = -999.;
    float nTrackerLayersWithMeasurement = -999.;
    float nPixelLayersWithMeasurement = -999.;
};

bool muon_legacy_sort (const muon_legacy_t& mA, const muon_legacy_t& mB)
{
  return (mA.pt > mB.pt);
}

bool muon_scouting_sort (const muon_scouting_t& mA, const muon_scouting_t& mB)
{
  return (mA.pt > mB.pt);
}

struct muonsv_legacy_t {
    float chi2 = -999.;
    float pAngle = -999.;
    float dlen = -999.;
    float dlenSig = -999.;
    float dxy = -999.;
    float dxySig = -999.;
    float mu1pt = -999.;
    float mu1eta = -999.;
    float mu1phi = -999.;
    float mu2pt = -999.;
    float mu2eta = -999.;
    float mu2phi = -999.;
    float x = -999.;
    float y = -999.;
    float z = -999.;
    float deltaR = -999.;
    //muonsv_legacy_t():
        //chi2(-999), pAngle(-999), dlen(-999), dlenSig(-999), dxy(-999), dxySig(-999),
        //mu1pt(-999), mu1eta(-999), mu1phi(-999), mu2pt(-999), mu2eta(-999), mu2phi(-999),
        //x(-999), y(-999), z(-999), deltaR(-999) {}
};

struct muonsv_scouting_t {
    float chi2 = -999.;
    float pAngle = -999.;
    float dlen = -999.;
    float dlenSig = -999.;
    float dxy = -999.;
    float dxySig = -999.;
    float mu1pt = -999.;
    float mu1eta = -999.;
    float mu1phi = -999.;
    float mu2pt = -999.;
    float mu2eta = -999.;
    float mu2phi = -999.;
    float x = -999.;
    float y = -999.;
    float z = -999.;
    float deltaR = -999.;
    float deltaEta = -999.;
    float deltaPhi = -999.;
    float pt = -999.;
};

bool muonsv_legacy_sort (const muonsv_legacy_t& svA, const muonsv_legacy_t& svB)
{
  return (svA.dlen > svB.dlen);
}

bool muonsv_ul_sort (const muonsv_legacy_t& svA, const muonsv_legacy_t& svB)
{
  return (svA.chi2 < svB.chi2) && (svA.chi2 != -999.);
}

bool muonsv_scouting_sort (const muonsv_scouting_t& svA, const muonsv_scouting_t& svB)
{
  return (svA.chi2 < svB.chi2) && (svA.chi2 != -999.);
}

struct sv_legacy_t {
    float pt = -999.;
    float eta = -999.;
    float phi = -999.;
    float mass = -999.;
    float x = -999.;
    float y = -999.;
    float z = -999.;
    float dxy = -999.;
    float dxySig = -999.;
    float dlen = -999.;
    float dlenSig = -999.;
    float pAngle = -999.;
    float chi2 = -999.;
    float ndof = -999.;
    //sv_legacy_t():
        //pt(-999), eta(-999), phi(-999), mass(-999),
        //x(-999), y(-999), z(-999), dxy(-999), dxySig(-999), dlen(-999), dlenSig(-999),
        //pAngle(-999), chi2(-999), ndof(-999) {}
};

struct sv_scouting_t {
    float mass = -999.;
    float x = -999.;
    float y = -999.;
    float z = -999.;
    float dxy = -999.;
    float dxySig = -999.;
    float dlen = -999.;
    float dlenSig = -999.;
    float pAngle = -999;
    float chi2 = -999.;
    float ndof = -999.;
    float ntracks = -999.;
};



bool sv_legacy_sort (const sv_legacy_t& svA, const sv_legacy_t& svB)
{
  return (svA.dlen > svB.dlen);
}

bool sv_ul_sort (const sv_legacy_t& svA, const sv_legacy_t& svB)
{
  return (svA.chi2 < svB.chi2) && (svA.chi2 != -999.);
}

bool sv_scouting_sort (const sv_scouting_t& svA, const sv_scouting_t& svB)
{
  return (svA.dlen > svB.dlen);
}


class BDTinference {
  public:
  BDTinference ();
    BDTinference (std::string filename);
    ~BDTinference ();
    std::vector<float> get_bdt_outputs(std::vector<float> inputs);

  private:
    BoosterHandle booster_;

};

#endif
