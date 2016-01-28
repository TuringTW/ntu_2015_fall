void randomN2(){
	// TH2F* hist = new TH2F("test", "test", 1000, -5, 5, 1000, -5, 5);
	TNtuple* ntuple = new TNtuple("test", "ntuple", "x:y");
	 TCanvas* c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
	TRandom *r = new TRandom();
	r->SetSeed(23);
	for (int i = 0; i < 400000; ++i){

		float x = r->Gaus(0, 1);
		float y = r->Gaus(0, 1);
		
		ntuple->Fill(x, y);
	}
	c1->Divide(2, 1);
	c1->cd(1);
	ntuple->Draw("x:y");
	// hist->Draw();
	TH1F* hist = new TH1F("hist", "hist", 1000, -5, 5);
	c1->cd(2);
	ntuple->Draw("x+y>>hist");
	hist->Draw("");
	hist->Fit("gaus");
}