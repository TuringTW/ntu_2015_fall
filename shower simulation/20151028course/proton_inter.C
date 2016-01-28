TRandom *r;
double lamda_proton = 87;
double phi = 0.5;
double step = 1/cos(phi); //cm

float density(float h){
	float rho = 1.02*exp(-h/8)*1e-3;
	return rho;
}
double pro_inter(float density, float lamda, float step){
	double pro = 1./lamda*step*density*1e5;
	return pro;
}
float get_inter_height(){
	float height = 100;
		while(r->Rndm() > pro_inter(density(height), lamda_proton, step)){
			height = height - 1;
		}
	return height;
}
void proton_inter(){
	TH1F* h1 = new TH1F("h1", "h1", 150, 0, 150);
	r = new TRandom();
	r->SetSeed(212312413);
	for (int i = 0; i < 10000; ++i){
		double height = get_inter_height();
		h1->Fill(height);
	}
	h1->Draw();
}
