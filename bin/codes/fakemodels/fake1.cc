// this is included int fakemodel.cc

void GetY(int iY,string Yname,vector<double> &X,double &Y,double &SigmaY){
	int ipar,NPars=X.size();
	Crandy *randy=new Crandy(iY+time(NULL));
	double Lambda=2.5,arg;
	arg=0.0;
	for(ipar=0;ipar<NPars;ipar++){
		arg+=randy->ran()*X[ipar]/(2.0*PI);
	}
	Y=50.0*randy->ran()*sin(arg/Lambda)+50.0*randy->ran()*cos(arg/Lambda);
	arg=0.0;
	for(ipar=0;ipar<NPars;ipar++){
		arg+=randy->ran()*X[ipar]/(2.0*PI);
	}
	Y+=50.0*randy->ran()*exp(arg/(Lambda*NPars));
	SigmaY=0.0;
	if(iY==0)
		printf("-----------------------\n");
	printf("%s  Y=%g +/- %g\n",Yname.c_str(),Y,SigmaY);	
}
