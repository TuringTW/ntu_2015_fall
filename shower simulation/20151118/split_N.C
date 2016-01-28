TRandom *r;
const int branches = 4;
void Nsplit(double all, double branch[]){
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
void split_N(){
	r = new TRandom();
	TH1F *h1 = new TH1F("h1", "h1", 100, 0, 1);

	double pizza_mass = 1; //kg
	for (int i = 0; i < 4000; ++i)
	{
		D
		double branch[branches];
		Nsplit(pizza_mass, branch);
		// for (int i = 0; i < branches; ++i)
		// {
		// 	h1->Fill(branch[i]);
		// }
		h1->Fill(branch[1]);
	}

	h1->Draw();
}