void pion_lifetime(){
	double energy = 1; //gev
	double M_pion = 0.13977; //gev
	double Pi_lifetime = 0.26033e-6; //s
	double c = 3e8;//m/s
	TCanvas *c1 = new TCanvas("test", "test");
	TRandom *r = new TRandom();
	TH1F *h1 = new TH1F("h1", "h1", 1000, 0, 5e-6);
	TH1F *h2 = new TH1F("h2", "h2", 1000, 0, 10);
	
	for (int i = 0; i < 10000; ++i)
	{
		//life time
		double t = r->Exp(Pi_lifetime);
		h1->Fill(t);
		//decay length

		double gamma = energy/M_pion;
		double length = c*t*gamma/1e3;
		h2->Fill(length);

	}
	c1->Divide(2,1);
	c1->cd(1);
	h1->Draw();
	c1->cd(2);
	h2->Draw();

}