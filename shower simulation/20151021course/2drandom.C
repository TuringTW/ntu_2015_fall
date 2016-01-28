int 2drandom(){

	TH2F* hist = new TH2F("test", "test", 100, -10, 10, 100, -10, 10);
	TRandom *r = new TRandom();
	r->SetSeed(212312413);
	for (int i = 0; i < 100000; ++i){

		float x, y;
		r->Rannor(x, y);
		// float y = pow(3*x, 1.0/3);
		
		hist->Fill(x, y);
	}

	hist->Draw("lego");
	hist->Fit("gaus", "", "", -10, 10);
}