TRandom *r;
void pizza_tax(){
	r = new TRandom();
	TNtuple *ntuple = new TNtuple("ntuple", "ntuple", "A:n");
	TCanvas *c1 = new TCanvas();
	TH1F *h1 = new TH1F("h1", "h1", 100, 0, 1);
	TH1F *npis = new TH1F("npis", "npis", 100, 0, 10);

	double tax = 0.05;//kg
	
	for (int i = 0; i < 1000; ++i)
	{
		double pizza_mass = 1; //kg
		int n_pis = 0;
		while(pizza_mass > 0){
			double frac = r->Rndm();
			double A = pizza_mass*frac;
			pizza_mass = pizza_mass*(1-frac)-tax;
			n_pis++;
			h1->Fill(A);		
				
		}
		npis->Fill(n_pis);	

	}
	c1->Divide(2,1);
	
	c1->cd(1);
	h1->Draw();

	c1->cd(2);
	npis->Draw();
}



/*
double ntuple

*/