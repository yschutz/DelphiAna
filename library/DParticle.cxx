//////////////////////////////////////////////////////////
// DParticle class: a delphi particle
// Y. Schutz 13 December 2018
//////////////////////////////////////////////////////////
#include "DParticle.h"
#include "DData.h"
#include <iostream>

//==========================================================================
DParticle::DParticle() : fDecayPartnerIndex(-1), fPdgCode(-1)
{
 // empty defaults ctor
  SetName(""); 
  SetTitle(""); 
}
//==========================================================================
DParticle::DParticle(Int_t pdgCode, const TLorentzVector &mom4, Int_t decay, Double_t vr, Double_t vz) :
    fMomentum4(mom4), fPdgCode(pdgCode), fVtxRphi(vr), fVtxZ(vz)
{
  // ctor
  if (decay != 0) fDecayPartnerIndex = decay -1; //carefull: particles index starts at 1 in Fortran and 0 in C++
  else fDecayPartnerIndex = decay; 
  SetNameTitle(); 
}

//==========================================================================
DParticle::DParticle(const DParticle& part) :
  TNamed(part),
  fDecayPartnerIndex(part.fDecayPartnerIndex),
  fMomentum4(part.fMomentum4),
  fPdgCode(part.fPdgCode),
  fVtxRphi(part.fVtxRphi),
  fVtxZ(part.fVtxZ)
{
  // Copy constructor
}

//==========================================================================
DParticle& DParticle::operator=(const DParticle& part)
{
  // Assignment operator
  if(this!=&part) {
    TNamed::operator=(part); 
    fDecayPartnerIndex = part.fDecayPartnerIndex;
    fMomentum4 = part.fMomentum4;
    fPdgCode = part.fPdgCode;
    fVtxRphi= part.fVtxRphi;
    fVtxZ = part.fVtxZ;
  }
  return *this;
}

//==========================================================================
void DParticle::Clear(Option_t *option)
{
  SetName(""); 
  SetTitle(""); 
  // clears all data 
  fPdgCode = 0; 
  fMomentum4.SetXYZM(0., 0., 0., -1.); 
}  
   
//==========================================================================
Double_t DParticle::Eta() const
{
  // pseudorapidity

  Double_t pmom = Momentum(); 
  Double_t pz   = Momentum3().Z();  
  if (pmom != TMath::Abs(pz)) 
    return 0.5*TMath::Log((pmom + pz) / (pmom - pz));
  else 
    return 1.e30;
}

//==========================================================================
void DParticle::Print(Int_t pindex) const
{
  // print data of particle 
  printf("part # = %3i %10s:%20s pdg = %6i  mass = %10.8f \n", pindex, GetName(), GetTitle(), fPdgCode, Mass());
}

//==========================================================================
void DParticle::Set(Int_t pdg, Double_t px, Double_t py, Double_t pz, Double_t mass, Int_t decay, Double_t vr, Double_t vz)
{
  // set the attributes of a Dparticle 
  if (decay != 0) fDecayPartnerIndex = decay - 1; //carefull: particles index starts at 1 in Fortran and 0 in C++
  else fDecayPartnerIndex = decay; 
  fPdgCode           = pdg;   
  fMomentum4.SetXYZM(px, py, pz, mass);   
  fVtxRphi = vr; 
  fVtxZ = vz; 
  SetNameTitle(); 
}

//==========================================================================
void DParticle::SetNameTitle()
{
  // set the name and the title according to fPdgCode
  TString name(""); 
  TString title("");  
      
  if (fPdgCode == 0) {
    name = TString("unidentified"); 
    title = TString("no identification available"); 
  } else if (fPdgCode == -11) {
    name = TString("e-"); 
    title = TString("electron"); 
  } else if (fPdgCode == 11) {
    name = TString("e+"); 
    title = TString("positron"); 
  } else if (fPdgCode == 13) {
    name = TString("mu+"); 
    title = TString("anti muon"); 
  } else if (fPdgCode == -13) {
    name = TString("mu-"); 
    title = TString("muon"); 
  } else if (fPdgCode == 22) {
    name = TString("gamma"); 
    title = TString("photon"); 
  } else if (fPdgCode == 11) {
    name = TString("pi0"); 
    title = TString("neutral pion"); 
  } else if (fPdgCode == 130) {
    name = TString("K0L"); 
    title = TString("K0 long"); 
  } else if (fPdgCode == 211) {
    name = TString("pi+"); 
    title = TString("anti charged pion"); 
  } else if (fPdgCode == -211) {
    name = TString("pi-"); 
    title = TString("charged pion"); 
  } else if (fPdgCode == 310) {
    name = TString("K0S"); 
    title = TString("K0 short"); 
  } else if (fPdgCode == 311) {
    name = TString("K0"); 
    title = TString("Neutral K"); 
  } else if (fPdgCode == 321) {
    name = TString("K+"); 
    title = TString("anti charged kaon"); 
  } else if (fPdgCode == -321) {
    name = TString("K-"); 
    title = TString("charged kaon"); 
  } else if (fPdgCode == 2121) {
    name = TString("n"); 
    title = TString("neutron"); 
  } else if (fPdgCode == 2212) {
    name = TString("p"); 
    title = TString("proton"); 
  } else if (fPdgCode == -2212) {
    name = TString("anti p"); 
    title = TString("anti proton"); 
  } else if (fPdgCode == 3122) {
    name = TString("Lambda"); 
    title = TString("lambda baryon"); 
  } else {
    name = TString("unknownino"); 
    title = TString("a new particle ?"); 
  }
  SetName(name); 
  SetTitle(title); 
}
