// ROOT includes
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TList.h>
#include <Riostream.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TProfile.h>

const char *names[3] = {"High Mult", "Intermediate Mult", "Low Mult"};

TH1D *ProjectDEta(TH2D *histo, Double_t mindeta, Double_t maxdeta);
void ScaleBins(TH2D *histo);

void analyzeCorr(const char *filename = "../data/histos.root",
		 const char *var = "Eta",
		 Int_t firstZbin = 0, Int_t lastZbin = 0,
		 Double_t minDeta = 1.2, Double_t maxDEta = 4.0,
		 Int_t sumMethod = 2,
		 Int_t extMethod = 4,
		 Int_t phiRebin = 1,
		 Int_t etaRebin = 1,
		 Int_t firstTrigBin = 0, //0,
		 Int_t firstAssocBin = 0,//0,
		 Double_t asScale = -1.0
		 )
{
  const Int_t nBinsMult = 3;
  const Int_t nBinsPt = 5;
  const Int_t nBinsZvtx = 1;
  
  TFile *f = TFile::Open(filename);

  TList *l = (TList *)f->Get("histos");
  TList *lMix = (TList *)f->Get("histos");
  if (!lMix) lMix = l;

  TH1D *histTrig[nBinsMult];
  TH1D *histTrigBin[nBinsMult][nBinsZvtx];
  for(Int_t i = 0; i < nBinsMult; ++i) {
    histTrig[i] = 0x0;
    for(Int_t j = 0; j < nBinsZvtx; ++j) {
      histTrigBin[i][j] = 0x0;
    }
  }

  {
    for(Int_t iCent = 0; iCent < nBinsMult; ++iCent) {
      for(Int_t j = firstZbin; j <= lastZbin; ++j) {
	      for(Int_t iPtBin = firstTrigBin; iPtBin < nBinsPt; ++iPtBin) {
          std::cout << "********** " << Form("fHistTrig_Mult%02d_Z%02d_PtBin%02d",iCent,j,iPtBin) << std::endl;
	        TH1D *histTrigTrk = (TH1D*)l->FindObject(Form("fHistTrig_Mult%02d_Z%02d_PtBin%02d",iCent,j,iPtBin));
	        histTrigTrk->Sumw2();
	        if (!histTrigBin[iCent][j])
	          histTrigBin[iCent][j] = (TH1D*)histTrigTrk->Clone(Form("histTrigBin_%d_%d",iCent,j));
	        else
	          histTrigBin[iCent][j]->Add(histTrigTrk);
	      }
      }
    }
  }

  TH2D *hDphiEtaSum[nBinsMult][nBinsZvtx];
  TH2D *hDphiEtaMixSum[nBinsMult][nBinsZvtx];
  for(Int_t i = 0; i < nBinsMult; ++i) {
    for(Int_t j = 0; j < nBinsZvtx; ++j) {
      hDphiEtaSum[i][j] = 0x0;
      hDphiEtaMixSum[i][j] = 0x0;
    }
  }

  {
    for(Int_t iCent = 0; iCent < nBinsMult; ++iCent) {
      for(Int_t j = firstZbin; j <= lastZbin; ++j) {
	      for(Int_t iBin = 0; iBin < nBinsPt; ++iBin) {
	        for(Int_t jBin = 0; jBin < nBinsPt; ++jBin) {
	          if (iBin >= firstTrigBin && jBin >= firstAssocBin) {
              TH2D *hDPhiEta = (TH2D*)l->FindObject(Form("fHistDPhi%s_Mult%02d_Z%02d_PtBin%02d_%02d",var,iCent,j,iBin,jBin));
	            hDPhiEta->Sumw2();
	            if (iBin >= TMath::Max(firstTrigBin,firstAssocBin) && jBin >= TMath::Max(firstTrigBin,firstAssocBin))
		            ScaleBins(hDPhiEta);
	            if (!hDphiEtaSum[iCent][j])
		            hDphiEtaSum[iCent][j] = (TH2D*)hDPhiEta->Clone(Form("hDphiEtaSum_%02d_%02d",iCent,j));
	            else
		            hDphiEtaSum[iCent][j]->Add(hDPhiEta);
	          }
	          if ((iBin >= firstTrigBin && jBin >= firstAssocBin) || (jBin >= firstTrigBin && iBin >= firstAssocBin)) {
	            TH2D *hDPhiEtaMix = (TH2D*)lMix->FindObject(Form("fHistDPhi%sMix_Mult%02d_Z%02d_PtBin%02d_%02d",var,iCent,j,iBin,jBin));
	            hDPhiEtaMix->Sumw2();
      	      if (!hDphiEtaMixSum[iCent][j])
		            hDphiEtaMixSum[iCent][j] = (TH2D*)hDPhiEtaMix->Clone(Form("hDphiEtaMixSum_%02d_%02d",iCent,j));
	            else
		            hDphiEtaMixSum[iCent][j]->Add(hDPhiEtaMix);
	          }
	        }
	      }
      }
    }
  }


  {
    for(Int_t iCent = 0; iCent < nBinsMult; ++iCent) {
      for(Int_t iz = firstZbin; iz <= lastZbin; ++iz) {
	      Double_t sum = 0;
	      if (sumMethod == 0) {
	  // Take sum at maximum acceptance in deta
	        for(Int_t k = 0; k <= (hDphiEtaMixSum[iCent][iz]->GetNbinsY()+1); ++k) {
	          Double_t sumPart = 0;
	          for(Int_t j = 1; j <= hDphiEtaMixSum[iCent][iz]->GetNbinsX(); ++j) {
		          sumPart += hDphiEtaMixSum[iCent][iz]->GetBinContent(j,k);
	          }
	          if (sumPart > sum) sum = sumPart;
	        }
	      }
	      if (sumMethod == 1) {
	  // Take sum as integral in deta
	        for(Int_t k = 0; k <= (hDphiEtaMixSum[iCent][iz]->GetNbinsY()+1); ++k) {
	          for(Int_t j = 1; j <= hDphiEtaMixSum[iCent][iz]->GetNbinsX(); ++j) {
		          sum += hDphiEtaMixSum[iCent][iz]->GetBinContent(j,k);
	          }
	        }
	      }
	      if (sumMethod == 2) {
	  // Take sum close to (0,0)
	        Int_t bin1x = hDphiEtaMixSum[iCent][iz]->GetXaxis()->FindBin(-1e-6);
	        Int_t bin2x = hDphiEtaMixSum[iCent][iz]->GetXaxis()->FindBin(+1e-6);
	        Int_t bin1y = hDphiEtaMixSum[iCent][iz]->GetYaxis()->FindBin(-1e-6);
	        Int_t bin2y = hDphiEtaMixSum[iCent][iz]->GetYaxis()->FindBin(+1e-6);
	        for(Int_t k = bin1y; k <= bin2y; ++k) {
	          for(Int_t j = bin1x; j <= bin2x; ++j) {
		          sum += hDphiEtaMixSum[iCent][iz]->GetBinContent(j,k);
	          }
	        }
	      }
	      printf("Sum = %f (%d %d)\n",sum,iCent,iz);
	      if (sum == 0)
	        exit(0);
	      else {
	  //	if (sum > 0) {
	        for(Int_t k = 0; k <= (hDphiEtaMixSum[iCent][iz]->GetNbinsY()+1); ++k) {
	          for(Int_t j = 1; j <= hDphiEtaMixSum[iCent][iz]->GetNbinsX(); ++j) {
	            hDphiEtaMixSum[iCent][iz]->SetBinContent(j,k,
						  histTrigBin[iCent][iz]->GetBinContent(1)/sum*hDphiEtaMixSum[iCent][iz]->GetBinContent(j,k));
	            hDphiEtaMixSum[iCent][iz]->SetBinError(j,k,
						  histTrigBin[iCent][iz]->GetBinContent(1)/sum*hDphiEtaMixSum[iCent][iz]->GetBinError(j,k));
	          }
	        }
	      }
       	hDphiEtaSum[iCent][iz]->Divide(hDphiEtaMixSum[iCent][iz]);
      }
    }
  }

  TH2D *hDphiEtaInt[nBinsMult];
  {
    for(Int_t iCent = 0; iCent < nBinsMult; ++iCent) {
      hDphiEtaInt[iCent] = (TH2D*)hDphiEtaSum[iCent][firstZbin]->Clone(Form("hDphiEtaInt_%02d",iCent));
      hDphiEtaInt[iCent]->Reset();
      Int_t nZeros = 0;
      for(Int_t j = 1; j <= hDphiEtaInt[iCent]->GetNbinsX(); ++j) {
	      for(Int_t k = 0; k <= (hDphiEtaInt[iCent]->GetNbinsY()+1); ++k) {
	        Double_t sum = 0;
	        Double_t sumErr = 0;
	        for(Int_t iz = firstZbin; iz <= lastZbin; ++iz) {
	          if (hDphiEtaSum[iCent][iz]->GetBinContent(j,k) > 0) {
	            sum += (hDphiEtaSum[iCent][iz]->GetBinContent(j,k)/hDphiEtaSum[iCent][iz]->GetBinError(j,k)/hDphiEtaSum[iCent][iz]->GetBinError(j,k));
	            sumErr += (1./hDphiEtaSum[iCent][iz]->GetBinError(j,k)/hDphiEtaSum[iCent][iz]->GetBinError(j,k));
	          }
	          else {
	            nZeros++;
	          }
	        }
	        if (sumErr > 0) {
	          hDphiEtaInt[iCent]->SetBinContent(j,k,sum/sumErr);
	          hDphiEtaInt[iCent]->SetBinError(j,k,TMath::Sqrt(2./sumErr));
	        }
	      }
	//	printf("Zeros %d %d -> %d\n",iCent,i,nZeros);
      }
    }
  }

  TH2D *hSubtracted2D = (TH2D*)hDphiEtaInt[0]->Clone("hSubtracted2D");
  hSubtracted2D->Add(hDphiEtaInt[nBinsMult-1],asScale);

  if (1) {
    TCanvas *c1 = new TCanvas("c1");
    c1->Divide(3,2);
    for(Int_t iCent = 0; iCent < nBinsMult; ++iCent) {
      c1->cd(iCent+1);
      if (strcmp(var,"Eta") == 0) hDphiEtaInt[iCent]->GetYaxis()->SetRangeUser(-maxDEta,maxDEta);
      hDphiEtaInt[iCent]->SetStats(0);
      hDphiEtaInt[iCent]->GetXaxis()->SetTitle("#Delta#varphi (rad)");
      if (strcmp(var,"Eta") == 0) hDphiEtaInt[iCent]->GetYaxis()->SetTitle("#Delta#eta");
      if (strcmp(var,"Theta") == 0) hDphiEtaInt[iCent]->GetYaxis()->SetTitle("#Delta#theta (rad)");
      hDphiEtaInt[iCent]->GetZaxis()->SetTitle("Assoc yield per trigger");
      hDphiEtaInt[iCent]->SetTitle(names[iCent]);
      hDphiEtaInt[iCent]->GetXaxis()->SetTitleOffset(1.5);
      hDphiEtaInt[iCent]->GetYaxis()->SetTitleOffset(1.5);
      hDphiEtaInt[iCent]->GetZaxis()->SetTitleOffset(1.5);
      hDphiEtaInt[iCent]->Draw("surf5");
    }
    c1->cd(4);
    if (strcmp(var,"Eta") == 0) hSubtracted2D->GetYaxis()->SetRangeUser(-maxDEta,maxDEta);
    //    if (strcmp(var,"Eta") == 0)
    hSubtracted2D->GetZaxis()->SetRangeUser(-0.02,1.05*hSubtracted2D->GetMaximum());
    hSubtracted2D->SetStats(0);
    hSubtracted2D->GetXaxis()->SetTitle("#Delta#varphi (rad)");
    if (strcmp(var,"Eta") == 0) hSubtracted2D->GetYaxis()->SetTitle("#Delta#eta");
    if (strcmp(var,"Theta") == 0) hSubtracted2D->GetYaxis()->SetTitle("#Delta#theta (rad)");
    hSubtracted2D->GetZaxis()->SetTitle("Assoc yield per trigger");
    hSubtracted2D->SetTitle("High Mult - Low Mult");
    hSubtracted2D->GetXaxis()->SetTitleOffset(1.5);
    hSubtracted2D->GetYaxis()->SetTitleOffset(1.5);
    hSubtracted2D->GetZaxis()->SetTitleOffset(1.5);
    hSubtracted2D->Draw("surf5");
    //    if (strcmp(var,"Eta") == 0) {
      c1->cd(5);
      TH2D *hSubtracted2Dzoom = (TH2D*)hSubtracted2D->Clone(Form("%s_zoomed",hSubtracted2D->GetName()));
      hSubtracted2Dzoom->SetTitle("High Mult - Low Mult (Zoomed)");
      hSubtracted2Dzoom->GetZaxis()->SetRangeUser(-0.02,0.15*hSubtracted2D->GetMaximum());
      hSubtracted2Dzoom->Draw("surf5");
      //    }
    
  }

  return;
  
  TCanvas *cc = new TCanvas("cc");
  TH1D *hSubtracted = ProjectDEta(hSubtracted2D,minDeta,maxDEta);
  hSubtracted->SetName("hSubtracted");

  TH1D *hDphi[nBinsMult];
  TCanvas *cTT = new TCanvas("cTT");
  cTT->cd(1);
  for(Int_t iCent = 0; iCent < nBinsMult; ++iCent) {
    hDphi[iCent] = ProjectDEta(hDphiEtaInt[iCent],minDeta,maxDEta);
    hDphi[iCent]->SetName(Form("hDPhi_%d",iCent));
  }
  delete cTT;

  if (1) {
    TCanvas *c2 = new TCanvas("c2");
    c2->Divide(3,2);
    for(Int_t iCent = 0; iCent < nBinsMult; ++iCent) {
      c2->cd(iCent+1);
      hDphi[iCent]->Draw();
    }
    c2->cd(4);
    hSubtracted->Draw();
  }

  {
    TCanvas *cTemp = new TCanvas;
    //    TF1 *funcB = new TF1("funcB","[0]+2.*[0]*[1]*(TMath::Cos(x)-1.)");
    TF1 *funcB = new TF1("funcB","[0]+[0]*[1]*TMath::Exp(-(TMath::Abs(x)-TMath::Pi())*(TMath::Abs(x)-TMath::Pi())/[2]/[2])");
    //    funcB->SetParLimits(1,-10.,10.);
    funcB->SetParameter(2,1.);
    hDphi[nBinsMult-1]->Fit(funcB,"I");
    TF1 *funcB2 = new TF1("funcB2","[0]",-0.5,0.5);
    funcB2->SetLineColor(kGreen);
    hDphi[nBinsMult-1]->Fit(funcB2,"+R");
    Double_t b = funcB2->GetParameter(0);
    Double_t sigma = funcB->GetParameter(2);
    TF1 *func = 0x0;
    if (extMethod == 1) func = new TF1("func","[0]+2.*([0]+[3])*[2]*TMath::Cos(2.*x)+2.*([0]+[3])*[1]*TMath::Cos(x)");
    if (extMethod == 2) func = new TF1("func","[0]+2.*([0]+[3])*[2]*TMath::Cos(2.*x)+2.*([0]+[3])*[1]*TMath::Cos(x)+2.*([0]+[3])*[4]*TMath::Cos(3.*x)");
    if (extMethod == 3) func = new TF1("func","[0]+2.*([0]+[3])*[2]*TMath::Cos(2.*x)+([0]+[3])*[1]*TMath::Exp(-(TMath::Abs(x)-TMath::Pi())*(TMath::Abs(x)-TMath::Pi())/[4]/[4])");
    if (extMethod == 4) func = new TF1("func","[0]+2.*([0]+[3])*[2]*TMath::Cos(2.*x)+([0]+[3])*[1]*TMath::Exp(-(TMath::Abs(x)-TMath::Pi())*(TMath::Abs(x)-TMath::Pi())/[4]/[4])+2.*([0]+[3])*[5]*TMath::Cos(3.*x)");
    func->FixParameter(3,-asScale*b);
    if ((extMethod == 3) || (extMethod == 4)) func->FixParameter(4,sigma);
    hSubtracted->Fit(func,"I");
    hSubtracted->Fit(func,"I");
    printf("Results => b0 = %f  a0 = %f  a2 = %f   v2 = %f +- %f\n", b,func->GetParameter(0),func->GetParameter(2)*(b+func->GetParameter(0)),
	                                                                   TMath::Sqrt(func->GetParameter(2)),0.5*func->GetParError(2)/TMath::Sqrt(func->GetParameter(2)));
    delete cTemp;
    TFile fOut("resultsCorr.root","UPDATE");
    hSubtracted->SetName(Form("%s_%d_%d_%.1f_%.1f_%d_%d_%d_%d_%d_%d_%.1f",
			      hSubtracted->GetName(),
			      firstZbin,lastZbin,
			      minDeta,maxDEta,
			      sumMethod,extMethod,
			      phiRebin,etaRebin,
			      firstTrigBin,firstAssocBin,
			      asScale));
    hSubtracted->Write();
    fOut.Close();
   }
}

TH1D *ProjectDEta(TH2D *histo, Double_t mindeta, Double_t maxdeta)
{
  TH1D *hOut = (TH1D*)histo->ProjectionX(Form("%s_px",histo->GetName()));
  hOut->Reset();

  for(Int_t i = 1; i <= histo->GetNbinsX(); ++i) {
    TH1D *hDEta = (TH1D*)histo->ProjectionY(Form("%s_%d",histo->GetName(),i),i,i);
    if (1) {
      for(Int_t j = 0; j <= (hDEta->GetNbinsX()/2); ++j) {
	      hDEta->SetBinContent(j,
			  hDEta->GetBinContent(j)+
			  hDEta->GetBinContent(hDEta->GetNbinsX()+1-j));
	      hDEta->SetBinError(j,
			   TMath::Sqrt(hDEta->GetBinError(j)*hDEta->GetBinError(j)+
				       hDEta->GetBinError(hDEta->GetNbinsX()+1-j)*hDEta->GetBinError(hDEta->GetNbinsX()+1-j)));
      }
      TF1 *func1 = new TF1(Form("func1_%s",hDEta->GetName()),"pol1");
      hDEta->Fit(func1,"Q","",-maxdeta+1e-6,-mindeta-1e-6);
      hOut->SetBinContent(i,func1->Integral(-maxdeta,-mindeta));
      hOut->SetBinError(i,func1->IntegralError(-maxdeta,-mindeta));
    }
    else {
      TF1 *func1 = new TF1(Form("func1_%s",hDEta->GetName()),"pol1");
      hDEta->Fit(func1,"Q","",mindeta+1e-6,maxdeta-1e-6);
      Double_t int1 = func1->Integral(mindeta,maxdeta);
      Double_t int1err = func1->IntegralError(mindeta,maxdeta);
      TF1 *func2 = new TF1(Form("func2_%s",hDEta->GetName()),"pol1");
      hDEta->Fit(func2,"Q+","",-maxdeta+1e-6,-mindeta-1e-6);
      Double_t int2 = func2->Integral(-maxdeta,-mindeta);
      Double_t int2err = func2->IntegralError(-maxdeta,-mindeta);
      Double_t errSum = 1./int1err/int1err + 1./int2err/int2err;
      hOut->SetBinContent(i,(int1/int1err/int1err+int2/int2err/int2err)/errSum);
      hOut->SetBinError(i,TMath::Sqrt(1./errSum));
    }
  }
  return hOut;
}

void ScaleBins(TH2D *histo)
{
  for(Int_t i = 1; i <= histo->GetNbinsX(); ++i) {
    for(Int_t j = 0; j <= (histo->GetNbinsY()+1); ++j) {
      // histo->SetBinContent(i,j,
      // 			   histo->GetBinContent(i,j)/2.);
      // histo->SetBinError(i,j,
      // 			 histo->GetBinError(i,j)/TMath::Sqrt(2.));
      histo->SetBinError(i,j,
			 histo->GetBinError(i,j)*TMath::Sqrt(2.));
    }
  }
}
