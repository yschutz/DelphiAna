//////////////////////////////////////////////////////////
// DEvent class: a delphi event
// Y. Schutz 20 December 2018
//////////////////////////////////////////////////////////
#ifndef DEvent_h
#define DEvent_h

#include <TNamed.h>

#include <map>

// class TClonesArray; 
class TFile; 
class TTree; 
class DEvent : public TNamed {
  
  public:
   DEvent(); 
   virtual ~DEvent();
   
   Double_t         Aplanarity() const            { return fSphericity[2] * 3. / 2.;}                                               
   Int_t            ChargedMultiplicity() const   { return fChargedMul; }
   Double_t         ChargedEnergy() const; 
   Long64_t         EventNumber() const           { return fEventNumber; }
   Int_t            Multiplicity() const          { return fMul; }
   Double_t         Oblateness() const            { return fThrust[1] - fThrust[2]; }
   TClonesArray *   Particles() const             { return fParticles; }
   void             SetECM(Double_t e)            { fECM = e; }
   void             SetEventNumber(Long64_t evt)  { fEventNumber = evt; } 
   void             SetMul(Int_t mul, Int_t cmul) { fMul = mul; fChargedMul = cmul; }
   void             SetSpericity(Double_t s1, Double_t s2, Double_t s3) 
                                                  { fSphericity[0] = s1; fSphericity[1] = s2; fSphericity[2] = s3;}
   void             SetThrust(Double_t t1, Double_t t2, Double_t t3) 
                                                  { fThrust[0] = t1; fThrust[1] = t2; fThrust[2] = t3;}
   Double_t         Sphericity() const            { return (fSphericity[1] + fSphericity[2]) * 3. / 2.;}                                               
   Double_t         Thrust() const                { return fThrust[0]; }                                                
   Double_t         ThrustMajor() const           { return fThrust[1]; }  

   private:
   Int_t                     fChargedMul;             //! charges particles multiplicity
   Double_t                  fECM;                    //! Center-of-mass energy of the collision) 
   Long64_t                  fEventNumber;            //! The event number in the TTree 
   Int_t                     fMul;                    //! total multiplicity
   TClonesArray              *fParticles;             //! list of particles 
   Double_t                  fSphericity[3];          //! componemts of sphericity
   Double_t                  fThrust[3];              //! components of thrust 


  ClassDef(DEvent,1);
};
#endif