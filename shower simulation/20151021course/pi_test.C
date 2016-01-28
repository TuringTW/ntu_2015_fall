void pi_test(){

	TCanvas* c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);

	TH1F* h1 = new TH1F("h1", "h1", 200, 0, 5000);
	TNtuple* ntuple1 = new TNtuple("ntuple1", "ntuple1", "x");
	TRandom *r = new TRandom();

	double rho_air = 1.2754e-03;
	double lamda_pi = 130;
	double step = 100; //cm
	double prob_interaction = 1./lamda_pi*step*rho_air;
	r->SetSeed(212312413);

	for (int i = 0; i < 100000; ++i){

		double travel_l = 0;

		while(r->Rndm() > prob_interaction){
			travel_l +=(step)*0.01;
		}
		
		ntuple1->Fill(travel_l);
	}
	c1->Divide(2, 1);
	c1->cd(1);
       	ntuple1->Draw("x>>h1");

	h1->Draw();


}
