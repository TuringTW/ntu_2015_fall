void random2(){
	float x=7;
	float pi = 3.14159266;
	TH1F* hist = new TH1F("test", "test", 100, 0, 1);

	for (int i = 0; i < 1000; ++i)
	{
		x = (pi*x) - int(pi*x);
		printf("x%d=%f\n", i+1, x);
	    hist->Fill(x);
	}
	hist->Draw();
}