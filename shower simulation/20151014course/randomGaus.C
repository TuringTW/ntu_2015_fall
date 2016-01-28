void randomGaus(){
	TH1F* hist = new TH1F("test", "test", 10000, -10, 10);
	TRandom *r = new TRandom();
	r->SetSeed(212312413);
	for (int i = 0; i < 100000; ++i){

		float x = r->Gaus(0, 3);
		// float y = pow(3*x, 1.0/3);
		
		hist->Fill(x);
	}

	hist->Draw();
	hist->Fit("gaus", "", "", -10, 10);
}