double M_pion = 0.13977; //gev
double M_muon = 0.105658;
double M_neu = 0;

void plot_particle(TLorentzVector particle, double vx, double vy, double vz, int color){
	TParticle *particle_Neu = new TParticle();
	particle_Neu->SetMomentum(particle);
	particle_Neu->SetProductionVertex(vx, vy, vz, 0);
	particle_Neu->SetLineColor(color); //blue 
	particle_Neu->Draw();
}
void drawing_pion_decay(){

	double Pi_lifetime = 0.26033e-6; //s
	double c = 3e8;//m/s
	TCanvas *c1 = new TCanvas("test", "test");
	TRandom *r = new TRandom();
	

	TView *view = TView::CreateView(1, 0, 0);
	view->ShowAxis();
	view->SetRange(-5, -5, 0, 5, 5, 10);
	
	

	for (int i = 0; i < 100; ++i)
	{
		double px = 0.3;
		double py = 0.3;
		double pz = -1; // GeV
		double energy = sqrt(px*px+py*py+pz*pz+M_pion*M_pion);
		
		//life time
		double t = r->Exp(Pi_lifetime);
		//decay length
		double gamma = energy/M_pion;
		double length = c*t*gamma/1e3;
		// decay position
		double vx = 0;
		double vy = 0;
		double vz = 10 - length;

		TLorentzVector pion(px, py, pz, energy);

		double decay_particle_mass[2] = {M_muon, M_neu};
		TGenPhaseSpace event;
		event.SetDecay(pion, 2, decay_particle_mass);
		event.Generate();

		TLorentzVector Muon = *(event.GetDecay(0));
		TLorentzVector Neu = *(event.GetDecay(1));
		
		plot_particle(Muon, vx, vy, vz, 2);
		plot_particle(Neu, vx, vy, vz, 4);


	}


}