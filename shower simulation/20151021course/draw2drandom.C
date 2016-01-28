int draw2drandom(){
	TCanvas* c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
	TH2F* hist = new TH2F("test", "test", 100, -10, 10, 100, -10, 10);
	TRandom *r = new TRandom();
	r->SetSeed(212312413);
	for (int i = 0; i < 100000; ++i){

		// float x, y;
		float x = r->Gaus(0, 2);
		float y = r->Gaus(0, 2);
		// r->Rannor(x, y);
		// float y = pow(3*x, 1.0/3);
		
		hist->Fill(x, y);
	}
	c1->Divide(3, 3);
	c1->cd(1);
	hist->Draw("cont");
	c1->cd(2);

	hist->Draw("box");
	// hist->Fit("");
}