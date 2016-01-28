TRandom *r;
void psplit(double all, double &a, double &b){
	double frac = r->Rndm();
	a = all*frac;
	b = all*(1-frac);
}
void pizza_spit(){
	r = new TRandom();
	TH1F *h1 = new TH1F("h1", "h1", 100, 0, 1);

	double pizza_mass = 1; //kg

	for (int i = 0; i < 10.q0; ++i)
	{
		double pis_A;
		double pis_B;
		psplit(pizza_mass, pis_A, pis_B);
		h1->Fill(pis_B);
	}

	h1->Draw();
}