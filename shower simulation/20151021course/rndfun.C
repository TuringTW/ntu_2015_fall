int rndfun()
{
	TCanvas *c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
	TH1F* random = new TH1F("random", "random", 1000, -10, 10);
	TF1* fun = new TF1("fun", "exp(0.1*x)*sin(x)+3", -10, 10);
	c1->Divide(1, 2);
	c1->cd(1);
	fun->Draw();

	for (int i = 0; i < 10000; ++i)
	{
		double x = fun->GetRandom();
		random->Fill(x);
	}
	c1->cd(2);
	random->Draw();
}