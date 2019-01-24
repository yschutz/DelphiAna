#include <iostream>
#include <TNamed.h>
#include <TString.h>
#include <TSystem.h>
#include "library/DData.h"
void DAna()
{
	const Int_t kVerbose = 0;
	DData *data = new DData("DXXGeV", "e+e- XX GeV", "../data/");
	Float_t ecm = data->GetEcm(); 
	char cecm[9] = "";
	sprintf(cecm, "%.3f", ecm);
	TString name(data->GetName());
	name = name.ReplaceAll("XX", cecm); 
	TString title(data->GetTitle()); 
	title = title.ReplaceAll("XX", cecm);
	std::cout << "************ " << title << "************ " << std::endl; 
	data->SetName(name); 
	data->SetTitle(title); 
	data->SetVerbosity(0);  

	std::cout << "DAna:......." << std::endl; 

	// * Loops over all events 
	// 	  ** Loops over all particles      (DData::kAll)
    //                  charged particles  (DData::kCharged)
    //                  hadrons            (DData::kHadrons)
    //                  charged hadrons    (DData::kChargedHadrons)
	//       *** task Control: DData::kControlHisto --> compares histo produced by Delphi software and Root (DData::kCharged)
	//                         DData::kControlDiff  --> makes relative difference between histo produced by Delphi software and Root (DData::kCharged)
	//       *** task Correlation: implementation in progress

	TString opt1(gSystem->Getenv("OPT1"));  
	TString opt2(gSystem->Getenv("OPT2")); 
	TString param(gSystem->Getenv("PARAM")); 
	DData::Eopt oo1; 
	DData::Epopt oo2 = DData::kAll; 
	DData::Eparam par; 
	if (opt1.Contains("CH"))
		oo1 = DData::kControlHisto; 
	else if (opt1.Contains("CD"))
		oo1 = DData::kControlDiff; 
	else if (opt1.Contains("CO"))
		oo1 = DData::kCorrelation; 
	else if (opt1.Contains("SH"))
		oo1 = DData::kSingleHisto; 
	else {
		printf("!! Unknown option %s\n", opt1.Data());
		exit(1);
	} 
	if (opt2.Contains("A"))
		oo2 = DData::kAll; 
	else if (opt2.Contains("CH"))
		oo2 = DData::kChargedHadrons; 
	else if (opt2.Contains("C"))
		oo2 = DData::kCharged; 
	else if (opt2.Contains("H"))
		oo2 = DData::kHadrons; 
	
	if (param.Contains("1"))
		par = DData::kThrust; 
	else if (param.Contains("2"))
		par = DData::kMultiplicity; 

	if (oo1 == DData::kCorrelation) {
		TString smulL(gSystem->Getenv("LOWM"));
		TString smulH(gSystem->Getenv("HIGHM"));
		data->SetMulBin(smulL.Atoi(), smulH.Atoi());
		TString mix(gSystem->Getenv("MIX"));
		data->SetMix(mix.Atoi()); 
	}
	data->Run(oo1, oo2, par);
	data->WriteOutput(); 
}