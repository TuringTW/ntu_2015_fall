
double pi = TMath::Pi();
double N_Av = 6.02E23;
double rho_earth = 5.52; //g/cm^3
double R_earth = 6.4E8 ;//cm
double R_Det = 1.0E5; //cm
double rho_Det = 0.9 ;//g/cm^3
double threshold_e = 1e17;

double Xsec_neutrino(double energy, int current){
	double Xsec_cc = 5.53E-36*pow(energy, 0.363);
	double Xsec_nc = 2.51E-36*pow(energy, 0.363);

	if (current==0)return Xsec_cc;
	if(current == 1)return Xsec_nc;
	if(current == 2)return Xsec_nc+Xsec_cc;
	else return 0;
}
double cord_length(double R, double theta){
	double kapa = 2*theta-pi;
	double c = 2*R*sin(kapa/2.);
	return c;
}

void neutrinoTest(){
	TRandom3 *r = new TRandom3();
	TNtuple *ntuple = new TNtuple("ntuple", "ntuple", "e:theta");
	TF1 *fun = new TF1("test", "exp(x)", 0, 20);
	
	

	for (int i = 0; i < 600000; ++i)
	{
		double theta = acos(r->Rndm()-1);
		double energy = 1.0E18;//GeV
		double theta_deg = theta*180./pi;
		double cord_length_cm = cord_length(R_earth, theta);
		double cord_length_x = cord_length_cm*rho_earth;

		double sigma = Xsec_neutrino(energy/1.0e9, 1); //only NC // in GeV
		double lambda = 1./(N_Av*sigma);
		double propagation_length = r->Exp(lambda);

		while(cord_length_x > propagation_length){
			energy = energy/2.;  //interaction inside the earth 
			double sigma = Xsec_neutrino(energy/1.0e9, 1); //only NC // in GeV
			double lambda = 1./(N_Av*sigma);
			double propagation_l_one_time = r->Exp(lambda);
			propagation_length = propagation_length + propagation_l_one_time;
			
		}

		if (energy>=threshold_e)
		{
			ntuple->Fill(energy, theta_deg);		
		}
		
		
	}
	TH1F *hist = new TH1F("hist", "hist", 180, 0, 180);
	ntuple->Draw("theta>>hist");
	hist->SetXTitle("theta(Degree)");
	hist->SetYTitle("Entry");

}

