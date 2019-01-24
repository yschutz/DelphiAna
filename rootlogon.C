{
    TString incpath(gSystem->GetWorkingDirectory()); 
    gSystem->AddIncludePath(incpath.Prepend("-I")); 
    incpath.Append("/library");
    gSystem->AddIncludePath(incpath); 
}