//////////////////////////////////////////////////////////
// DParticle class: a delphi particle
// Y. Schutz 13 December 2018
//////////////////////////////////////////////////////////

#ifndef DParticle_h
#define DParticle_h

#include <TNamed.h>
#include <TLorentzVector.h> 
#include <TVector3.h>

class DParticle : public TNamed {
  
  public:
   DParticle(); 
   DParticle(Int_t pdgCode, const TLorentzVector &mom4, Int_t decay, Double_t vr, Double_t vz);
   virtual ~DParticle() {;}
   DParticle(const DParticle& part); 
   DParticle& operator=(const DParticle& part);

   virtual void Clear(Option_t *option ="");
   Double_t Energy()  const   { return fMomentum4.E(); } 
   Double_t Eta() const; 
   Double_t Mass() const      { return fMomentum4.M(); }
   TVector3 Momentum3() const { return fMomentum4.Vect(); } 
   Double_t Momentum()  const { return fMomentum4.P(); }
   Double_t Pt() const        { return fMomentum4.Pt(); }
   Int_t    PdgCode() const   { return fPdgCode; }
   Double_t Phi()const        { return TMath::RadToDeg() * fMomentum4.Phi(); } 
   Double_t PhiRad()const     { return fMomentum4.Phi(); } 
   void     Print(Int_t pindex) const;  
   void     Set(Int_t pdg, Double_t px, Double_t py, Double_t pz, Double_t mass, Int_t decay, Double_t vr, Double_t vz);
   Double_t Theta()const      { return TMath::RadToDeg() * fMomentum4.Theta(); } 
   Double_t ThetaRad()const      { return fMomentum4.Theta(); } 
   Double_t VtxR() const      { return fVtxRphi; }
   Double_t VtxZ() const      { return fVtxZ; }

   private:
   void     SetNameTitle();  

   Int_t                     fDecayPartnerIndex; //! index of the decay partner of this particle (-1 means it is a prompt paticle)
   TLorentzVector            fMomentum4;         //! particle 4-momentum
   Int_t                     fPdgCode;           //! particle pdg code 
   Double_t                  fVtxRphi;           //! RPhi position of the origine vertex
   Double_t                  fVtxZ;              //! Z position of the origine vertex

  ClassDef(DParticle,1);
};
#endif
