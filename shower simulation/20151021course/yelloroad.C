int yelloroad(){
	TH1F *yr = new TH1F("yellowroad", "yellowroad", 100, 0, 50); 
	TRandom *r = new TRandom();
	r->SetSeed();
	for (int i = 0; i < 10000; ++i)
	{

		int x = 0;
		double n = r->Rndm();
		while(n>0.5){
			x++;
			n = r->Rndm();
		}
		yr->Fill(x);

	}
	yr->Draw();
	yr->Fit("expo");
}