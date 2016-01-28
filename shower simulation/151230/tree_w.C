const int MAX_N = 1000;

typedef struct event_t
{
	int evt;
	int n_mu;
	float theta;
	float phi;
	float e;
	float x[MAX_N];
	float y[MAX_N];
	float z[MAX_N];
	float t[MAX_N];
};

TRandom * r;

void tree_w(){
	r = new Random();
	TFile *f1 = new TFile("muon.root", "recreate");
	TTree tree_muon("muon_tree", "apple tree");
	event_t muon_event;
	tree_muon.Branch("evt", &muon_event.evt, "evt/I");
	tree_muon.Branch("n_mu", &muon_event.n_mu, "n_mu/I");
	tree_muon.Branch("e", &muon_event.e, "e/F");
	tree_muon.Branch("theta", &muon_event.theta, "theta/F");
	tree_muon.Branch("phi", &muon_event.phi, "phi/F");
	tree_muon.Branch("x", muon_event.x, "x[n_mu]/F");
	tree_muon.Branch("y", muon_event.y, "y[n_mu]/F");
	tree_muon.Branch("z", muon_event.z, "z[n_mu]/F");
	tree_muon.Branch("t", muon_event.t, "t[n_mu]/F");


	for (int ievt = 0; ievt < 100; ++ievt)
	{
		int n_mu = r ->Rndm()*100;
		int N_obs = 0;

		while(N_obs < n_mu){
			float x_muon_impact = r ->Rndm();
			float y_muon_impact = r ->Rndm();
			float z_muon_impact = r ->Rndm();
			float t_obs = r ->Rndm();

			muon_event.x[]
		}
	}
}