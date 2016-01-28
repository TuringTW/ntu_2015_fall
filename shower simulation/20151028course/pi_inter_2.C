double myfun(double *x, double *par){
  return exp(-x[0]*par[0]/par[1]);
}
void pi_inter_2(){
	TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
	TH1F* random = new TH1F("random", "random", 200, 0, 5000);
	// TF1* fun = new TF1("myfun", myfun, 0, 5000, 2);
	TF1* fun = new TF1("myfun", "exp(-x)", 0, 10);

	double rho_air = 1.2754e-03; 
	double lamda_pi = 130; //cm

	// fun->SetParameters(rho_air, lamda_pi);

	c1->Divide(1, 2);
	c1->cd(1);
	fun->Draw();

	for (int i = 0; i < 10000; ++i)
	{
		double x = fun->GetRandom()*lamda_pi;
		double x_cm = x/rho_air;
		double x_m = x_cm/100;
		random->Fill(x_m);
	}
	c1->cd(2);
	random->Draw();
	random->Fit("expo");
}
