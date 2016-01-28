TRandom3 *r;
int pdg_proton_p = 2112;
int pdg_proton_n = -2112;
int pdg_pion_p = 211;
int pdg_pion_0 = 111;
int pdg_pion_n = -211;
int pdg_muon_p = -13;
int pdg_muon_n = 13;
int pdg_neutrino = 14;

double tau_pion = 0.26033e-6; //s

const int n_branches = 4;

double M_proton = 0.938;
double M_pi = 0.13977;
double M_pi0 = 0.13498;
double M_mu = 0.105658;
double M_neu = 0;
double M_e = 511e-06; //GeV

double rho0 = 1.021/1000; //g/cm3
double A=8; //g/cm3   ==>air density

const double pi = TMath::Pi();
TNtuple *ntuple1 = new TNtuple("ntuple1", "ntuple1", "ea:eb");

struct particle
{
	int id;
	int flag;
	double e;
	double px;
	double py;
	double pz;
	double vertex_x;
	double vertex_y;
	double vertex_z;
	double fx;
	double fy;
	double fz;
};
const int bank_size = 1e5;
particle Ptcl_bank[bank_size];
int n_particle = 0;

void Ptcl_Register(int id, int flag, double px, double py, double pz, double x, double y, double z){
	double mass = 0;
	
	if (id == pdg_pion_0) mass = M_pi0;
	if (id == pdg_pion_p) mass = M_pi;
	if (id == pdg_proton_p ) mass = M_proton;


	Ptcl_bank[n_particle].id = id;
	Ptcl_bank[n_particle].flag = flag;
	Ptcl_bank[n_particle].e = sqrt(px*px+py*py+pz*pz+mass*mass);
	Ptcl_bank[n_particle].px = px;
	Ptcl_bank[n_particle].py = py;
	Ptcl_bank[n_particle].pz = pz;
	Ptcl_bank[n_particle].vertex_x = x;
	Ptcl_bank[n_particle].vertex_y = y;
	Ptcl_bank[n_particle].vertex_z = z;
	

	n_particle++;

}
double X2H(double x, double h0, double theta){
	double costheta = fabs(cos(theta));
	double h;
	if (costheta!=0)
	{
		h = -A*log(x*costheta/(rho0*1e5*A)+exp(-h0/A));
	}else{
		h = h0;
	}
	return h;
}

void Get_Int_Posi(int id, double theta, double phi,double x, double y, double z, double &vertex_x, double &vertex_y, double &vertex_z){
	double lambda = 0;
	if (id == pdg_proton_p){lambda = 87;} //cm
	else if (id == pdg_pion_p){lambda = 116;} //cm
	else{
		cout << id << endl;
	}
	double X = r->Exp(lambda);
	double h0 = z;
	vertex_z = X2H(X, h0, theta);
	double radius = vertex_z*tan(theta);
	vertex_x = radius*cos(phi);
	vertex_y = radius*sin(phi);
}

void Nsplit(double all, double branch[], int branches){
	double frac[branches];
	double frac_sum = 0;
	for (int i = 0; i < branches; ++i)
	{
		frac[i] = r->Rndm();
		frac_sum += frac[i];
	}

	for (int i = 0; i < branches; ++i)
	{
		branch[i] = all*frac[i]/frac_sum;
	}
}
void twosplit(double all, double &a, double &b){
	double frac = r->Rndm();
	a = all*frac;
	b = all*(1-frac);
}

void Hillas_Split(double id, double E, double theta, double phi, double vertex_x, double vertex_y, double vertex_z){

	double E_init_proton = E; // GeV
	double E_proton = 0;
	double E_init_pion = 0;

	double mass = 0;
	if (id == pdg_pion_p) mass = M_pi;
	else mass = M_proton;
	
	while(E_proton<mass){
		twosplit(E_init_proton, E_proton, E_init_pion);
	}
	
	double p = sqrt(E_proton*E_proton-mass*mass);
	double px = -p*sin(theta)*cos(phi);
	double py = -p*sin(theta)*sin(phi);
	double pz = -p*cos(theta);


	Ptcl_Register(id, 1, px, py, pz, vertex_x, vertex_y, vertex_z);

	double pion_branch[n_branches];
	
	Nsplit(E_init_pion, pion_branch, n_branches);
	for (int i = 0; i < n_branches; ++i)
	{
		double energy_in = pion_branch[i];
		while(energy_in > 0){
			double E_a, E_b;
			twosplit(energy_in, E_a, E_b);
			
			energy_in = E_b - M_pi;
			if (energy_in>0)
			{
				// ?
				double px = -E_a*sin(theta)*cos(phi);
				double py = -E_a*sin(theta)*sin(phi);
				double pz = -E_a*cos(theta);

				Ptcl_Register(pdg_pion_p, 1, px, py, pz, vertex_x, vertex_y, vertex_z);
			}
		}
	}

	
}
double get_pion_decay_length(double Pi_lifetime, double energy){
		double t = r->Exp(Pi_lifetime);
		double c = 3e8;
		double gamma = energy/M_pi;
		double length = c*t*gamma/1e3;
		return length;
}


void recursive_inter(){
	TCanvas *c1 = new TCanvas("test", "test", 600, 700);
	TView *view = TView::CreateView(1, 0, 0);
	view->ShowAxis();
	view->SetRange(-10, -10, 0, 10, 10, 20);

	TNtuple *ntuple = new TNtuple("ntuple", "ntuple", "id:flag:e:px:py:pz:x:y:z");
	r = new TRandom3();
	r ->SetSeed(12232);
	for (int ievertex_nt = 0; ievertex_nt < 1; ievertex_nt++)
	{
		n_particle = 0;

		double E = 10;
		double theta = 0;
		// double theta = acos(r->Rndm());
		double phi = r->Rndm()*pi*2;
		double x = 0, y = 0, z = 1.0e10; //infinity
		double vertex_z, vertex_x, vertex_y; //km

		Get_Int_Posi(pdg_proton_p, theta, phi, x, y, z, vertex_x, vertex_y, vertex_z);

		Hillas_Split(pdg_proton_p, E, theta, phi, vertex_x, vertex_y, vertex_z);

		double proton_threshold = 10;
		double pion_threshold = 10;
		int i_particle = 0;
		while(i_particle < n_particle){
			int id = Ptcl_bank[i_particle].id;
			int flag = Ptcl_bank[i_particle].flag;
			double e = Ptcl_bank[i_particle].e;
			double x = Ptcl_bank[i_particle].vertex_x;
			double y = Ptcl_bank[i_particle].vertex_y;
			double z = Ptcl_bank[i_particle].vertex_z;
			if (id == pdg_proton_p && flag == 1 && e > proton_threshold)
			{
				Get_Int_Posi(id, theta, phi, x, y, z, vertex_x, vertex_y, vertex_z);
				if (vertex_z > 0)
				{
					Hillas_Split(id, e, theta, phi, vertex_x, vertex_y, vertex_z);
					Ptcl_bank[i_particle].flag = 0;
				}
			}
			if (id == pdg_pion_p &&flag == 1)
			{

				Get_Int_Posi(id, theta, phi, x, y, z, vertex_x, vertex_y, vertex_z)	;			
				double l_to_interaction = sqrt((vertex_x-x)*(vertex_x-x)+(vertex_y-y)*(vertex_y-y)+(vertex_z-z)*(vertex_z-z));
				double l_to_decay = get_pion_decay_length(tau_pion, e);
				
				// cout <<"e=" << e << "d=" << l_to_decay << " i=" << l_to_interaction << endl;

				if (l_to_decay>l_to_interaction)
				{
					if (vertex_z > 0)
					{
						Hillas_Split(id, e, theta, phi, vertex_x, vertex_y, vertex_z);
						Ptcl_bank[i_particle].flag = 0;
					}
				}else{
					Ptcl_bank[i_particle].flag = 0;
					double px = Ptcl_bank[i_particle].px;
					double py = Ptcl_bank[i_particle].py;
					double pz = Ptcl_bank[i_particle].pz; // GeV
					TLorentzVector pion(px, py, pz, e);
					double decay_particle_mass[2] = {M_mu, M_neu};
					TGenPhaseSpace event;
					event.SetDecay(pion, 2, decay_particle_mass);
					event.Generate();
					cout << "dsfaadf" << endl;
					TLorentzVector Muon = *(event.GetDecay(0));
					TLorentzVector Neu = *(event.GetDecay(1));
					//Muon
					TParticle *particle_Muon = new TParticle();
					particle_Muon->SetMomentum(Muon);
					particle_Muon->SetProductionVertex(vertex_x, vertex_y, vertex_z, 0);
					particle_Muon->SetLineColor(2);
					particle_Muon->Draw();
					// Neu
					TParticle *particle_Neu = new TParticle();
					particle_Neu->SetMomentum(Neu);
					particle_Neu->SetProductionVertex(vertex_x, vertex_y, vertex_z, 0);
					particle_Neu->SetLineColor(4);
					particle_Neu->Draw();

					// Ptcl_Register(pdg_muon_p, 0, Muon.Px(), Muon.Py(), Muon.Pz(), vertex_x, vertex_y, vertex_z);
					// Ptcl_Register(pdg_neutrino, 0, Neu.Px(), Neu.Py(), Neu.Pz(), vertex_x, vertex_y, vertex_z);


				}
					
			}
			i_particle++;
		}
		for (int i_particle = 0; i_particle < n_particle; ++i_particle){
			int id = Ptcl_bank[i_particle].id;
			int flag = Ptcl_bank[i_particle].flag;
			double e = Ptcl_bank[i_particle].e;
			double px = Ptcl_bank[i_particle].px;
			double py = Ptcl_bank[i_particle].py;
			double pz = Ptcl_bank[i_particle].pz;
			double x = Ptcl_bank[i_particle].vertex_x;
			double y = Ptcl_bank[i_particle].vertex_y;
			double z = Ptcl_bank[i_particle].vertex_z;
			ntuple->Fill(id, flag, e, px, py, pz, x, y, z);

			// TLorentzVector particle_temp(px, py, pz, e);
			// TParticle *particle = new TParticle();
			// particle->SetMomentum(particle_temp);
			// particle->SetProductionVertex(vertex_x, vertex_y, vertex_z, 0);
			// if (id==pdg_proton_p){	particle->SetLineColor(2);}
			// if (id==pdg_pion_p){	particle->SetLineColor(4);}
			// if (id==pdg_neutrino){	particle->SetLineColor(2);particle->Draw();}
			// if (id==pdg_muon_p){	particle->SetLineColor(4);particle->Draw();}
			


		}
	}
	// ntuple->Draw("z");
} 
