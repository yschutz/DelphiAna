#ifndef dpool_h
#define dpool_h
#include <vector>
#include "TH2D.h"
#include "TClonesArray.h"
#include "DParticle.h"

class DPool{
 public:
  DPool(Int_t minMult, Int_t fullMult):fIsReady(0),fIsFull(0),fMult(0),fMinMult(minMult),fFullMult(fullMult),fTracks(),fCurrentEvent(0){}
  ~DPool(){}

  void UpdatePool(TClonesArray* array){
    
    // create TClonesArrays if pool is not full yet
    if (!fIsFull) {
      fTracks.push_back(new TClonesArray("DParticle",10));
      fCurrentEvent = fTracks.size()-1;
    }
    
    // clear current event
    fTracks[fCurrentEvent]->Clear("C"); 
    
    // fill TClonesArray for current event
    Int_t nTracks=0;
    for (Int_t i=0;i<array->GetEntriesFast();i++) {
      DParticle* tr = (DParticle*) array->UncheckedAt(i);
      new ((*fTracks[fCurrentEvent])[nTracks++]) DParticle(*tr);
    }
    
    if (fIsFull) {
      // if pool is full, just increment event counter 
      fCurrentEvent++; if (fCurrentEvent>=fTracks.size()) fCurrentEvent=0;
    } else {
      // if pool is not full, update multiplicity counters
      fMult = fTracks.size();//+= array->GetEntriesFast();
      if (fMult>fMinMult) fIsReady=1;
      if (fMult>fFullMult) fIsFull=1;
    }
  }
  Bool_t        IsReady()           { return fIsReady;       }
  Int_t         GetCurrentNEvents() { return fTracks.size(); }
  TClonesArray* GetEvent(Int_t i)   { return fTracks[i];     }
 private:
  Bool_t       fIsReady;
  Bool_t       fIsFull;
  Int_t        fMult;
  Int_t        fMinMult;
  Int_t        fFullMult;
  std::vector<TClonesArray*> fTracks;
  UInt_t       fCurrentEvent;
};

class DPoolManager{
 public:
  DPoolManager(const char *poolName, Int_t minMult, Int_t fullMult, Int_t nCent, Double_t* vCent, Int_t nZvtx, Double_t* vZvtx){
    fHist = new TH2D(Form("hCentZvtxHist_%s",poolName),"",nCent,vCent,nZvtx,vZvtx);
    fPools = new DPool**[nCent];
    for (Int_t ic=0;ic<nCent;ic++) {
      fPools[ic] = new DPool*[nZvtx];
      for (Int_t iz=0;iz<nZvtx;iz++) {
        fPools[ic][iz] = new DPool(minMult,fullMult);
      }
    }
  }
  ~DPoolManager(){}
  DPool* GetEventPool(Double_t cent, Double_t zvtx){
    Int_t ic = fHist->GetXaxis()->FindBin(cent)-1;
    Int_t iz = fHist->GetYaxis()->FindBin(zvtx)-1;
    if (ic<0 || ic>fHist->GetNbinsX()) return 0;
    if (iz<0 || iz>fHist->GetNbinsY()) return 0;
    return fPools[ic][iz];
  }
  void PrintInfo(){
    printf("====\n");
    for (Int_t ic=0;ic<fHist->GetNbinsX();ic++) {
      for (Int_t iz=0;iz<fHist->GetNbinsY();iz++) {
        printf("%i  ",fPools[ic][iz]->GetCurrentNEvents());
        Int_t nTracks = 0;
        for (Int_t i=0;i<fPools[ic][iz]->GetCurrentNEvents();i++) nTracks+=fPools[ic][iz]->GetEvent(i)->GetEntriesFast();
//        printf("%i  ",nTracks);
      }
      printf("\n");
    }
  }
 private:
  DPool*** fPools;
  TH2D* fHist;
};

#endif
