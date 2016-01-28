TRandom *r;
double M_pion = 0.13977;
int n_branches = 4;

void Nsplit(double all, double branch[], int branches){
	double frac[branches];
	double frac_sum = 0;
	for (int i = 0; i < branches; ++i)
	{
		frac[i] = r->Rndm();
		frac_sum += frac[i];
	}

	for (int i = 0; i < branches; ++i)
	{
		branch[i] = all*frac[i]/frac_sum;
	}
}
void twosplit(double all, double &a, double &b){
	double frac = r->Rndm();
	a = all*frac;
	b = all*(1-frac);
}
void proton_pion(){
	r = new TRandom();
	TCanvas *c1 = new TCanvas();
	TH1F *Epion = new TH1F("Epion", "Epion", 100, 0, 140);
	TH1F *Npion = new TH1F("Npion", "Npion", 100, 0, 50);

	
	for (int i = 0; i < 100000; ++i)
	{
		double E_init_proton = 1000; // GeV
		double E_proton;
		double E_init_pion;

		twosplit(E_init_proton, E_proton, E_init_pion);
		
		double pion_branch[n_branches];
		
		Nsplit(E_init_pion, pion_branch, n_branches);
		int n_pion = 0;
		for (int i = 0; i < n_branches; ++i)
		{
			double energy_in = pion_branch[i];
			while(energy_in > 0){
				double E_a, E_b;
				twosplit(energy_in, E_a, E_b);
				double pion_E = E_a + M_pion;
				Epion->Fill(pion_E);

				energy_in = E_b - M_pion;
				if (energy_in>0)
				{
					n_pion++;
				}
				// double frac = r->Rndm();
				// double A = E_init_pion*frac;
				// pizza_mass = pizza_mass*(1-frac)-tax;
				// n_pis++;
				// h1->Fill(A);		
					
			}

		}
		Npion->Fill(n_pion);
			
	}
	c1->Divide(2,1);
	
	c1->cd(1);
	Epion->Draw();

	c1->cd(2);
	Npion->Draw();
	Npion->Fit("gaus");
}



/*
double ntuple

*/