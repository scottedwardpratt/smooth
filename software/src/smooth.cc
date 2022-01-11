#include "emulator.h"
using namespace std;
CRandy CSmooth::randy=NULL;

CSmooth::CSmooth(unsigned int MaxRank_set,unsigned int NPars_set,unsigned int NTrainingPoints_set){
	int j,rank,isame;
	NPars=NPars_Set;
	MaxRank=MaxRank_Set;
	NTrainingPoints=NTrainingPoints_Set;
	if(MaxRank>5){
		printf("MaxRank=%d is too big, being reset to 5\n",MaxRank);
		MaxRank=5;
	}
	vector<int> i,countsame;
	countsame.resize(MaxRank+1);
	i.resize(MaxRank+1);
	factorial.resize(MaxRank+1);
	factorial[0]=factorial[1]=1;
	for(j=2;j<=MaxRank;j++)
		factorial[j]=j*factorial[j-1];

	A1.resize(NPars);
	A2.resize(NPars);
	A3.resize(NPars);
	A4.resize(NPars);
	A5.resize(NPars);
	dupfactor1.resize(NPars);
	dupfactor2.resize(NPars);
	dupfactor3.resize(NPars);
	dupfactor4.resize(NPars);
	dupfactor5.resize(NPars);
	for(i[0]=0;i[0]<NPars;i[0]++){
		A2[i[0]].resize(i[0]+1);
		A3[i[0]].resize(i[0]+1);
		A4[i[0]].resize(i[0]+1);
		A5[i[0]].resize(i[0]+1);
		dupfactor2[i[0]].resize(i[0]+1);
		dupfactor3[i[0]].resize(i[0]+1);
		dupfactor4[i[0]].resize(i[0]+1);
		dupfactor5[i[0]].resize(i[0]+1);
		A1[i[0]]=0.0;
		rank=1;
		for(j=0;j<rank;j++)
			countsame[j]=0;
		isame=0;
		countsame[isame]=1;
		for(j=1;j<rank;j++){
			if(i[j]!=i[j-1])
				isame+=1;
			countsame[isame]+=1;
		}
		dupfactor1[i[0]]=factorial[rank];
		for(j=0;j<rank;j++)
			dupfactor1[i[0]]/=factorial[countsame[j]];
		
		
		//
		for(i[1]=0;i[1]<=i[0];i[1]++){
			A3[i[0]][i[1]].resize(i[1]+1);
			A4[i[0]][i[1]].resize(i[1]+1);
			A5[i[0]][i[1]].resize(i[1]+1);
			dupfactor3[i[0]][i[1]].resize(i[1]+1);
			dupfactor4[i[0]][i[1]].resize(i[1]+1);
			dupfactor5[i[0]][i[1]].resize(i[1]+1);
			A2[i[0]][i[1]]=0.0;
			rank=2;
			for(j=0;j<rank;j++)
				countsame[j]=0;
			isame=0;
			countsame[isame]=1;
			for(j=1;j<rank;j++){
				if(i[j]!=i[j-1])
					isame+=1;
				countsame[isame]+=1;
			}
			dupfactor2[i[0]][i[1]]=factorial[rank];
			for(j=0;j<rank;j++){
				dupfactor2[i[0]][i[1]]/=factorial[countsame[j]];
			}

			//
			for(i[2]=0;i[2]<=i[1];i[2]++){
				A4[i[0]][i[1]][i[2]].resize(i[2]+1);
				A5[i[0]][i[1]][i[2]].resize(i[2]+1);
				dupfactor4[i[0]][i[1]][i[2]].resize(i[2]+1);
				dupfactor5[i[0]][i[1]][i[2]].resize(i[2]+1);
				A3[i[0]][i[1]][i[2]]=0.0;
				rank=3;
				for(j=0;j<rank;j++)
					countsame[j]=0;
				isame=0;
				countsame[isame]=1;
				for(j=1;j<rank;j++){
					if(i[j]!=i[j-1])
						isame+=1;
					countsame[isame]+=1;
				}
				dupfactor3[i[0]][i[1]][i[2]]=factorial[rank];
				for(j=0;j<rank;j++)
					dupfactor3[i[0]][i[1]][i[2]]/=factorial[countsame[j]];
				
				//
				for(i[3]=0;i[3]<=i[2];i[3]++){
					A5[i[0]][i[1]][i[2]][i[3]].resize(i[3]+1);
					dupfactor5[i[0]][i[1]][i[2]][i[3]].resize(i[3]+1);
					A4[i[0]][i[1]][i[2]][i[3]]=0.0;
					rank=4;
					for(j=0;j<rank;j++)
						countsame[j]=0;
					isame=0;
					countsame[isame]=1;
					for(j=1;j<rank;j++){
						if(i[j]!=i[j-1])
							isame+=1;
						countsame[isame]+=1;
					}
					dupfactor4[i[0]][i[1]][i[2]][i[3]]=factorial[rank];
					for(j=0;j<rank;j++)
						dupfactor4[i[0]][i[1]][i[2]][i[3]]/=factorial[countsame[j]];
					//					
					for(i[4]=0;i[4]<=i[3];i[4]++){
						A5[i[0]][i[1]][i[2]][i[3]][i[4]]=0.0;
						rank=5;
						for(j=0;j<rank;j++)
							countsame[j]=0;
						isame=0;
						countsame[isame]=1;
						for(j=1;j<rank;j++){
							if(i[j]!=i[j-1])
								isame+=1;
							countsame[isame]+=1;
						}
						dupfactor5[i[0]][i[1]][i[2]][i[3]][i[4]]=factorial[rank];
						for(j=0;j<rank;j++)
							dupfactor5[i[0]][i[1]][i[2]][i[3]][i[4]]/=factorial[countsame[j]];
					}
				}
			}
		}
	}
}

double CSmooth::CalcY(vector<double> &theta){
	vector<int> i;
	double lambda2=LAMBDA*LAMBDA;
	double lambda3=lambda2*LAMBDA;
	double lambda4=lambda2*lambda2;
	double lambda5=lambda3*lambda2;
	lambda2=lambda2*factorial[2];
	lambda3=lambda3*factorial[3];
	lambda4=lambda4*factorial[4];
	lambda5=lambda5*factorial[5];
	i.resize(MaxRank+1);
	double F=A0;
	for(i[0]=0;i[0]<NPars;i[0]++){
		F+=dupfactor1[i[1]]*A1[i[0]]
			*theta[i[0]]/LAMBDA;
		for(i[1]=0;i[1]<=i[0];i[1]++){
			F+=dupfactor2[i[0]][i[1]]*A2[i[0]][i[1]]
				*theta[i[0]]*theta[i[1]]/lambda2;
			for(i[2]=0;i[2]<=i[1];i[2]++){
				F+=dupfactor3[i[0]][i[1]][i[2]]*A3[i[0]][i[1]][i[2]]
					*theta[i[0]]*theta[i[1]]*theta[i[2]]/lambda3;
				for(i[3]=0;i[3]<=i[2];i[3]++){
					F+=dupfactor4[i[0]][i[1]][i[2]][i[3]]*A4[i[0]][i[1]][i[2]][i[3]]
						*theta[i[0]]*theta[i[1]]*theta[i[2]]*theta[i[3]]/lambda4;
					for(i[4]=0;i[4]<=i[3];i[4]++){
						F+=dupfactor5[i[0]][i[1]][i[2]][i[3]][i[4]]*A5[i[0]][i[1]][i[2]][i[3]][i[4]]
							*theta[i[0]]*theta[i[1]]*theta[i[2]]*theta[i[3]]*theta[i[4]]/lambda5;
					}
				}
			}
		}
	}
	return F;
}

void CSmooth::GetPreFactors(unsigned int NTrainingPts,vector<double> &theta,vector<double> &Prefactor){
	if(Prefactor.size()!=NTrainingPts){
		Prefactor.resize(NTrainingPts)
	}
	vector<int> i;
	double lambda2=LAMBDA*LAMBDA;
	double lambda3=lambda2*LAMBDA;
	double lambda4=lambda2*lambda2;
	double lambda5=lambda3*lambda2;
	lambda2=lambda2*factorial[2];
	lambda3=lambda3*factorial[3];
	lambda4=lambda4*factorial[4];
	lambda5=lambda5*factorial[5];
	i.resize(MaxRank+1);
	unsigned int ip=0;
	Prefactor[ip]=1.0;
	ip+=1;
	for(i[0]=0;i[0]<NPars;i[0]++){
		if(ip<NTrainingPts)
			Prefactor[ip]=dupfactor1[i[1]]
				*theta[i[0]]/LAMBDA;
		ip+=1;
		if(ip>=NTrainingPts)
			break;
		for(i[1]=0;i[1]<=i[0];i[1]++){
			if(ip<NTrainingPts)
				Prefactor[ip]=dupfactor2[i[0]][i[1]]
					*theta[i[0]]*theta[i[1]]/lambda2;
			ip+=1;
			if(ip>=NTrainingPts)
				break;
			for(i[2]=0;i[2]<=i[1];i[2]++){
				if(ip<NTrainingPts)
					Prefactor[ip]=dupfactor3[i[0]][i[1]][i[2]]
						*theta[i[0]]*theta[i[1]]*theta[i[2]]/lambda3;
				ip+=1;
				if(ip>=NTrainingPts)
					break;
				for(i[3]=0;i[3]<=i[2];i[3]++){
					if(ip<NTrainingPts)
						Prefactor[ip]=dupfactor4[i[0]][i[1]][i[2]][i[3]]
							*theta[i[0]]*theta[i[1]]*theta[i[2]]*theta[i[3]]/lambda4;
					ip+=1;
					if(ip>=NTrainingPts)
						break;
					for(i[4]=0;i[4]<=i[3];i[4]++){
						if(ip<NTrainingPts)
							Prefactor[ip]=dupfactor5[i[0]][i[1]][i[2]][i[3]][i[4]]
							*theta[i[0]]*theta[i[1]]*theta[i[2]]*theta[i[3]]*theta[i[4]]/lambda5;
						ip+=1;
						if(ip>=NTrainingPts)
							break;
					}
				}
			}
		}
	}
}

void CSmooth::SetParsFromX(vector<double> x){
	unsigned int NTrainingPts=x.size();
	vector<int> i;
	i.resize(MaxRank+1);
	unsigned int ip=0;
	A0=x[ip];
	ip+=1;
	for(i[0]=0;i[0]<NPars;i[0]++){
		if(ip<NTrainingPts)
			A1[i[1]]=x[ip]
		ip+=1;
		if(ip>=NTrainingPts)
			break;
		for(i[1]=0;i[1]<=i[0];i[1]++){
			if(ip<NTrainingPts)
				A2[i[0]][i[1]]=x[ip];
			ip+=1;
			if(ip>=NTrainingPts)
				break;
			for(i[2]=0;i[2]<=i[1];i[2]++){
				if(ip<NTrainingPts)
					A3[i[0]][i[1]][i[2]]=x[ip];
				ip+=1;
				if(ip>=NTrainingPts)
					break;
				for(i[3]=0;i[3]<=i[2];i[3]++){
					if(ip<NTrainingPts)
						A4[i[0]][i[1]][i[2]][i[3]]=x[ip];
					ip+=1;
					if(ip>=NTrainingPts)
						break;
					for(i[4]=0;i[4]<=i[3];i[4]++){
						if(ip<NTrainingPts)
							A5[i[0]][i[1]][i[2]][i[3]][i[4]]=x[ip];
						ip+=1;
						if(ip>=NTrainingPts)
							break;
					}
				}
			}
		}
	}
}

void CSmooth::SetA_Zero(){
	vector<int> i,countsame;
	i.resize(MaxRank+1);
	A0=0.0;
	for(i[0]=0;i[0]<NPars;i[0]++){
		A1[i[0]]=0.0;
		for(i[1]=0;i[1]<=i[0];i[1]++){
			A2[i[0]][i[1]]=0.0;
			for(i[2]=0;i[2]<=i[1];i[2]++){
				A3[i[0]][i[1]][i[2]]=0.0;
				for(i[3]=0;i[3]<=i[2];i[3]++){
					A4[i[0]][i[1]][i[2]][i[3]]=0.0;
					for(i[4]=0;i[4]<=i[3];i[4]++){
						A5[i[0]][i[1]][i[2]][i[3]][i[4]]=0.0;
					}
				}
			}
		}
	}
}

void CSmooth::SetA_RanGauss(double Amag){
	vector<int> i,countsame;
	i.resize(MaxRank+1);
	A0=Amag*randy->ran_gauss();
	for(i[0]=0;i[0]<NPars;i[0]++){
		A1[i[0]]=Amag*randy->ran_gauss();
		for(i[1]=0;i[1]<=i[0];i[1]++){
			A2[i[0]][i[1]]=Amag*randy->ran_gauss();
			for(i[2]=0;i[2]<=i[1];i[2]++){
				A3[i[0]][i[1]][i[2]]=Amag*randy->ran_gauss();
				for(i[3]=0;i[3]<=i[2];i[3]++){
					A4[i[0]][i[1]][i[2]][i[3]]=Amag*randy->ran_gauss();
					for(i[4]=0;i[4]<=i[3];i[4]++){
						A5[i[0]][i[1]][i[2]][i[3]][i[4]]=Amag*randy->ran_gauss();
					}
				}
			}
		}
	}
}

void CSmooth::SetA_RanSech(double Amag){
	// proportional to 1/cosh(A/Amag)
	double r;
	vector<int> i,countsame;
	i.resize(MaxRank+1);
	r=randy->ran();
	A0=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
	for(i[0]=0;i[0]<NPars;i[0]++){
		r=randy->ran();
		A1[i[0]]=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
		for(i[1]=0;i[1]<=i[0];i[1]++){
			r=randy->ran();
			A2[i[0]][i[1]]=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
			for(i[2]=0;i[2]<=i[1];i[2]++){
				r=randy->ran();
				A3[i[0]][i[1]][i[2]]=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
				for(i[3]=0;i[3]<=i[2];i[3]++){
					r=randy->ran();
					A4[i[0]][i[1]][i[2]][i[3]]=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
					for(i[4]=0;i[4]<=i[3];i[4]++){
						r=randy->ran();
						A5[i[0]][i[1]][i[2]][i[3]][i[4]]=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
					}
				}
			}
		}
	}
}

void CSmooth::SetA_Constant(double Amag){
	vector<int> i,countsame;
	i.resize(MaxRank+1);
	A0=Amag;
	for(i[0]=0;i[0]<NPars;i[0]++){
		A1[i[0]]=Amag;
		for(i[1]=0;i[1]<=i[0];i[1]++){
			A2[i[0]][i[1]]=Amag;
			for(i[2]=0;i[2]<=i[1];i[2]++){
				A3[i[0]][i[1]][i[2]]=Amag;
				for(i[3]=0;i[3]<=i[2];i[3]++){
					A4[i[0]][i[1]][i[2]][i[3]]=Amag;
					for(i[4]=0;i[4]<=i[3];i[4]++){
						A5[i[0]][i[1]][i[2]][i[3]][i[4]]=Amag;
					}
				}
			}
		}
	}
}


//Note that if I have n1 of the same theta and n2 of the same theta,,  s.t. n1+n2+n3.. = Nrank, then the number of combinations is Nrank!/(n1!n2!n3!...). Because we are dividing Nrank!, our counting factor is 1/(n1!n2!n3!...)