void profile_muon_press(){
	TCanvas * c1 = new TCanvas("c1", "c1", 1200, 500);	
	TFile *f1 = TFile::Open("muon.root");
	// TNtuple *ntuple = (TNtuple*)gDirectory->Get("ntuple");
	// TNtuple *ntuple = (TNtuple*)f1->Get("ntuple");
	
	c1 ->Divide(2, 0);
	c1 ->cd(1);
	ntuple ->SetMarkerSize(25);
	ntuple ->Draw("nobs:press*1000");

	c1 ->cd(2);
	TProfile *p1 = new TProfile("p1", "p1", 20, 0.9, 1.1, 0, 2000);
	ntuple ->Draw("nobs:press*1000>>p1");
	p1->Fit("pol1");

}