double M_e = 0.105658;
double M_mu = 511e-6;
double dedx_muon(double *x, double *par){
	double E = x[0];
	int material_code = par[0];
	// cout << material_code << endl;

	double K = 0.307e-3;
	double I = 85.7e-9;
	double ZA = 0.49919;
	double delta = 0;
	if (material_code == 1)
	{
		I = 136.3e-9;
		ZA = 0.50;
	}

	double z = 1;
	double p = sqrt(E*E - M_mu*M_mu);
	double gamma = sqrt(1.+pow(p/M_mu, 2));
	double beta = sqrt(1.-1./(gamma*gamma));
	double g_b = p/M_mu;
	double tmax = 2*M_e*p*p/(M_mu*M_mu + M_e*M_e + 2*M_e*E);
	double dedx = K*z*z*ZA/(beta*beta)*(1/2.*log(2*M_e*g_b*g_b*tmax/(I*I)) - beta*beta -delta/2.);
	cout << dedx << endl;
	return dedx;
}

void draw_dedx(){
	TF1 *f1 = new TF1("f1", dedx_muon, M_mu, 1000, 1);
	TF1 *f2 = new TF1("f2", dedx_muon, M_mu, 1000, 1);

	f1->SetParameter(0, 0);
	f1->Draw();
	f2->SetParameter(0, 1);
	f2->SetLineColor(4);
	f2->Draw("same");


}