void randomN1(){
	TH2F* hist = new TH2F("test", "test", 1000, -5, 5, 1000, -5, 5);
	TRandom *r = new TRandom();
	r->SetSeed(212312413);
	for (int i = 0; i < 100000; ++i){

		float x = r->Gaus(0, 1);
		float y = r->Gaus(0, 1);
		
		hist->Fill(x, y);
	}

	hist->Draw();
}