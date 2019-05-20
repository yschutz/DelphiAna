//////////////////////////////////////////////////////////
// DEvent class: a delphi event
// Y. Schutz 20 December 2018
//////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TLorentzVector.h>    
#include <TString.h>

#include "DEvent.h"
#include "DParticle.h"

//==========================================================================
DEvent::DEvent()
{
    // empty default ctor
    
    fParticles = new TClonesArray("DParticle", 100); 
    SetName(""); 
    SetTitle(""); 
}

//==========================================================================
DEvent::~DEvent()
{
    // dtor
    fParticles->Delete(); 
    delete fParticles; 
}

//==========================================================================
Double_t DEvent::ChargedEnergy() const
{
   // calculates the charged energy assuming all charged particles
   TLorentzVector mom4;
   Double_t ene = 0.0;
   Int_t index = 0; 
   TIter next(fParticles); 
   while (DParticle *part = (DParticle*)next()) {
      if (index > ChargedMultiplicity()) break; 
      ene += part->Energy(); 
      index++; 
   }
   return ene / fECM;
} 

//==========================================================================
void DEvent::Clear(Option_t *option)
{
  // clear event
  SetName("");
  SetTitle("");
  fBeamSpotTrans = -1e9;
  fBeamSpotLong = -1e9;
  fChargedMul = -1;
  fECM = -1e9;
  fEventNumber = -1;
  fMul = -1;
  fParticles->Clear("C");
  fSphericity[0] = fSphericity[1] = fSphericity[2] = -1e9;
  fThrust[0] = fThrust[1] = fThrust[2] = -1e9;
}
