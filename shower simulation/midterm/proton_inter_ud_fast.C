TRandom *r;
void proton_inter_ud_fast(){
	TH1F* random = new TH1F("random", "random", 100, 0, 10);
	TF1* fun = new TF1("myfun", "exp(-x)", 0, 10);
	r = new TRandom();
	r->SetSeed(403);

	double lamda_proton = 87; //cm
	double phi = 70*3.14/180;
	for (int i = 0; i < 10000; ++i)
	{
		double x = r->Exp(1)*lamda_proton;
		double x_cm = -8*1e5*log(-(x/(1.02*1e-3*1e5*8)-1));
		double x_km = x_cm/1e5;
		random->Fill(x_km);
		printf("%f\n", x_km);
	}
	random->Draw();
	// random->Fit("expo");
}
