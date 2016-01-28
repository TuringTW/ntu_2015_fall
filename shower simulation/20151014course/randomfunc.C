
void randomfunc(){
	TH1F* hist = new TH1F("test", "test", 100, 0, 1);
	TRandom *r = new TRandom();
	r->SetSeed(212312413);
	for (int i = 0; i < 1000000000; ++i)
	{

		float x = r->Rndm();
	    hist->Fill(x);
	}
	hist->Draw();
}