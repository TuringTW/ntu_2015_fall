TRandom *r;
double lamda_proton = 87; //g/cm^2
double phi = TMath::Pi()/4; //45degree
double step = 0.001; //km


double density(double h){
	double rho = 1.02*exp(-h/8)*1e-3;
	return rho;
}
double pro_inter(double density, double lamda, double step){
	double pro = (1./lamda)*step*density*1e5;
	return pro;
}
double get_inter_height(){
	double height = 0;
		while(r->Rndm() > pro_inter(density(height), lamda_proton, step/cos(phi))){
			height = height + step;
		}
	return height;
}
void proton_inter(){
	TH1F* h1 = new TH1F("h1", "h1", 2000, 0, 10);
	r = new TRandom();
	r->SetSeed(4343);
	for (int i = 0; i < 40000; ++i){
		double height = get_inter_height();
		h1->Fill(height);
	}
	h1->Draw();
	h1->Fit("gaus");
	h1->SetXTitle("Height(Km)");
	h1->SetYTitle("Entry");
}
