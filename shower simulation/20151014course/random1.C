void random1(){
	int x=7;
	int m=100;
	int a=2;
	for (int i = 0; i < 100; ++i)
	{
		x = (a*x)%m;
		printf("x%d=%d\n", i+1, x);
	}
}