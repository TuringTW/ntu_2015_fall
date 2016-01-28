void test(){
  TFile *_file0 = TFile::Open("muon.root");
  ntuple->Draw("nobs");



}
