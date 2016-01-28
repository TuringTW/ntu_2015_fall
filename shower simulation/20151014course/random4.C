void random4(){
	TH1F* hist = new TH1F("test", "test", 10000, -10, 10);
	TH1F* hist2 = new TH1F("test", "test", 10000, -10, 10);
	TRandom *r = new TRandom();
	r->SetSeed(212312413);
	for (int i = 0; i < 100000000; ++i){

		float x = r->Rndm()*20-10;
		float g = exp(-x*x/2.);
		// float s = sin(x);
		float y = r->Rndm();
		if (y < g){
		    hist->Fill(x);
		}
		// if(y < s){
		// 	hist2->Fill(x);
		// }
	}
	hist->Draw();
	// hist->Fit("fit");
	// hist2->Draw();
}