double M_e = 0.105658;
double M_mu = 511e-6;

TF1 *fn_muon_dxde;


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
	return dedx;
}
double dxde_muon(double *x, double *par){
	double dedx = dedx_muon(x, par);
	double dxde = 1./dedx;
	return dxde;
}
double get_muon_range(double E, int material_code){
	fn_muon_dxde->SetParameters(0, material_code);
	if (E<0.2){return 0;} // ignore low energy
	else{
		double R = fn_muon_dxde->Integral(0.2, E);
		return R;
	}
}
void draw_depth(){
	TRandom * r = new TRandom();	

	TH1F * muon_depth = new TH1F("muon_depth", "muon_depth", 100, 0, 500);
	fn_muon_dxde = new TF1("f1", dxde_muon, M_mu, 1000, 1);

	for (int i = 0; i < 1000; ++i)
	{
		double E_mu = 30; // gev
		int material_code = 1; //rock
		double range_muon = get_muon_range(E_mu, material_code);

		double X = r ->Exp(range_muon);
		double density_rock = 2.65; //g/cm3
		double penetration_depth = X/density_rock/100;

		muon_depth ->Fill(penetration_depth);

	}

	muon_depth ->Draw();
}