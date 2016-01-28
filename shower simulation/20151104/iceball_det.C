double pi = TMath::Pi();
double N_Av = 6.02E23;
double rho_earth = 5.52; //g/cm^3
double R_earth = 6.4E8 ;//cm
double R_Det = 1.0E5; //cm
double rho_Det = 0.9 ;//g/cm^3

double Xsec_neutrino(double energy, int current){
	double Xsec_cc = 5.53E-36*pow(energy, 0.363);
	double Xsec_nc = 2.51E-36*pow(energy, 0.363);

	if (current==0)return Xsec_cc;
	else if(current == 1)return Xsec_nc;
	else if(current == 2)return Xsec_nc+Xsec_cc;
	// else return 0;
}
double cord_length(double R, double theta){
	double kapa = 2*theta - pi;
	double c = 2*R*sin(kapa/2.);
	return c;
}
double Probability_inter(double x, double lambda){
	return (1-exp(-x/lambda));
}
void iceball_det(){

	TRandom3 *r = new TRandom3();
	TNtuple *ntuple = new TNtuple("ntuple", "ntuple", "e:theta:weight");
	TF1 *fun = new TF1("test", "exp(x)", 0, 20);
	for (int i = 0; i < 4; ++i)
	{
		
		double energy = 1.0E07*pow(10, i);
		double sigma = Xsec_neutrino(energy, 2);
		double lambda = 1./(N_Av*sigma);

		for (int i = 0; i < 400000; ++i)
		{
			double theta = acos(r->Rndm()-1);
			
			double theta_deg = theta*180./pi;
			double cord_length_cm = cord_length(R_earth, theta);
			double cord_length_x = cord_length_cm*rho_earth;

			double probability_passing_earth = 1 - Probability_inter(cord_length_x, lambda);
			double x_det = 2*R_Det*rho_Det;
			double probability_detection = Probability_inter(x_det, lambda);

			double Probability_evt = probability_passing_earth*probability_detection;

			ntuple->Fill(energy, theta_deg, Probability_evt);
		}
		
		
	}
	TH1F *hist1= new TH1F("hist1", "hist1", 100, 0, 180);
	TH1F *hist2 = new TH1F("hist2", "hist2", 100, 0, 180);
	TH1F *hist3 = new TH1F("hist3", "hist3", 100, 0, 180);
	TH1F *hist4 = new TH1F("hist4", "hist4", 100, 0, 180);
	
	ntuple->Draw("theta>>hist1", "(e==1.0E07)*weight");
	ntuple->Draw("theta>>hist2", "(e==1.0E08)*weight");
	ntuple->Draw("theta>>hist3", "(e==1.0E09)*weight");
	ntuple->Draw("theta>>hist4", "(e==1.0E10)*weight");
	
	hist1->SetLineColor(1);
	hist1->Draw("h");
	hist2->SetLineColor(2);
	hist2->Draw("h same");
	hist3->SetLineColor(3);
	
	hist3->Draw("h same");
	hist4->SetLineColor(4);
	hist4->Draw("h same");

}

