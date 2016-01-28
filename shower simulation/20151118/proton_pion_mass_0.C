TRandom *r;
double M_pion = 1;
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
void proton_pion_mass_0(){
	r = new TRandom();
	TCanvas *c1 = new TCanvas();
	TH1F *Epion = new TH1F("Epion", "Epion", 100, 0, 140);
	shift = 2700;
	TH1F *Npion_0 = new TH1F("Npion", "Npion", 100, shift, shift + 500);
	shift = 0;
	TH1F *Npion_1 = new TH1F("Npion", "Npion", 100, shift + 0, shift + 50);
	shift = 0;
	TH1F *Npion_2 = new TH1F("Npion", "Npion", 100, shift + 0, shift + 50);

	M_pion = 0;	
	for (int i = 0; i < 40000; ++i)
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
					// cout << n_pion <<endl;
				}
			}
		}
		Npion_0->Fill(n_pion);
	}
	M_pion = 0.5;	
	for (int i = 0; i < 40000; ++i)
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
			}
		}
		Npion_1->Fill(n_pion);
	}
	M_pion = 1;	
	for (int i = 0; i < 40000; ++i)
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
			}
		}
		Npion_2->Fill(n_pion);
	}

	c1->Divide(3,1);
	
	c1->cd(1);
	Npion_0->Draw();
	Npion_0->Fit("gaus");

	c1->cd(2);
	Npion_1->Draw();
	Npion_1->Fit("gaus");

	c1->cd(3);
	Npion_2->Draw();
	Npion_2->Fit("gaus");

}



/*
double ntuple

*/