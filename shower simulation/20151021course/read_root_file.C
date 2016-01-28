void read_root_file(){
	TFile *f = TFile::Open("ntuple.root");
	ntuple->Draw("x");
}