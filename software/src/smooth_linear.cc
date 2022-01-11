#include "emulator.h"
#include "smooth.h"
using namespace std;

CRandy CSmooth_Linear::randy=NULL;

CSmooth_Linear::CSmooth_Linear(unsigned int MaxRank_set,unsigned int NPars_set,unsigned int NTrainingPoints_set){
	unsigned int ic,rank,j,isame;
	vector<unsigned int> countsame;
	vector<unsigned int> factorial;
	vector<unsigned int> i;
	NTrainingPoints=NTrainingPoints_set;
	NPars=NPars_set;
	MaxRank=MaxRank_set;
	factorial.resize(MaxRank+1);
	factorial[0]=factorial[1]=1;
	for(j=2;j<=MaxRank;j++)
		factorial[j]=j*factorial[j-1];
	if(MaxRank>5){
		printf("MaxRank=%d is too big, being reset to 5\n",MaxRank);
		MaxRank=5;
	}
	i.resize(MaxRank+1);
	countsame.resize(MaxRank);
	
	ic=0;
	
	A.push_back();
	dupfactor.push_back(1.0);
	IPar.push_back();
	ic+=1;
	
	rank=1;
	for(i[0]=0;i[0]<NPars;i[0]++){
		A.push_back();
		dupfactor.push_back();
		IPar.push_back();
		IPar[ic].resize(rank);
		for(ir=0;ir<rank;ir++){
			IPar[ic][ir]=i[ir];
		}
		for(j=0;j<rank;j++)
			countsame[j]=0;
		isame=0;
		countsame[isame]=1;
		for(j=1;j<rank;j++){
			if(i[j]!=i[j-1])
				isame+=1;
			countsame[isame]+=1;
		}
		dupfactor[ic]=factorial[rank];
		for(j=0;j<rank;j++)
			dupfactor[ic]/=factorial[countsame[j]];
		ic+=1;
	}
	
	rank=2;
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			A.push_back();
			dupfactor.push_back();
			IPar.push_back();
			IPar[ic].resize(rank);
			for(ir=0;ir<rank;ir++){
				IPar[ic][ir]=i[ir];
			}
			for(j=0;j<rank;j++)
				countsame[j]=0;
			isame=0;
			countsame[isame]=1;
			for(j=1;j<rank;j++){
				if(i[j]!=i[j-1])
					isame+=1;
				countsame[isame]+=1;
			}
			dupfactor[ic]=factorial[rank];
			for(j=0;j<rank;j++)
				dupfactor[ic]/=factorial[countsame[j]];
			ic+=1;
		}
	}
	
	rank=3;
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				A.push_back();
				dupfactor.push_back();
				IPar.push_back();
				IPar[ic].resize(rank);
				for(ir=0;ir<rank;ir++){
					IPar[ic][ir]=i[ir];
				}
				for(j=0;j<rank;j++)
					countsame[j]=0;
				isame=0;
				countsame[isame]=1;
				for(j=1;j<rank;j++){
					if(i[j]!=i[j-1])
						isame+=1;
					countsame[isame]+=1;
				}
				dupfactor[ic]=factorial[rank];
				for(j=0;j<rank;j++)
					dupfactor[ic]/=factorial[countsame[j]];
				ic+=1;
			}
		}
	}
	
	rank=4;
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					A.push_back();
					dupfactor.push_back();
					IPar.push_back();
					IPar[ic].resize(rank);
					for(ir=0;ir<rank;ir++){
						IPar[ic][ir]=i[ir];
					}
					for(j=0;j<rank;j++)
						countsame[j]=0;
					isame=0;
					countsame[isame]=1;
					for(j=1;j<rank;j++){
						if(i[j]!=i[j-1])
							isame+=1;
						countsame[isame]+=1;
					}
					dupfactor[ic]=factorial[rank];
					for(j=0;j<rank;j++)
						dupfactor[ic]/=factorial[countsame[j]];
					ic+=1;
				}
			}
		}
	}
	
	rank=5;
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					for(i[4]=0;i[4]<=i[3];i[4]++){
						A.push_back();
						dupfactor.push_back();
						IPar.push_back();
						IPar[ic].resize(rank);
						for(ir=0;ir<rank;ir++){
							IPar[ic][ir]=i[ir];
						}
						for(j=0;j<rank;j++)
							countsame[j]=0;
						isame=0;
						countsame[isame]=1;
						for(j=1;j<rank;j++){
							if(i[j]!=i[j-1])
								isame+=1;
							countsame[isame]+=1;
						}
						dupfactor[ic]=factorial[rank];
						for(j=0;j<rank;j++)
							dupfactor[ic]/=factorial[countsame[j]];
						ic+=1;
					}
				}
			}
		}
	}
	NCoefficients=ic;
	if(NCoefficients!=A.size()){
		printf("size mismatch\n");
		exit(1);
	}
}

double CSmooth_Linear::CalcY(vector<double> theta){
	int ic,ir;
	double answer=0.0,term;
	for(ic=0;ic<NCoefficients;ic++){
		term=A[ic]*dupfactor[ic];
		for(ir=0;ir<IPar[ic].size();ir++){
			term*=theta[IPar[ic][ir]]
		}
		answer+=term;
	}
	return answer;
}

void CSmooth_Linear::Export(CSmooth *smooth){
	unsigned int ic,ir,rank;
	double answer=0.0,term;
	smooth->A0=A[0];
	for(ic=0;ic<NCoefficients;ic++){
		rank=IPar[ic].size();
		if(rank==0)
			smooth->A0=A[ic];
		if(rank==1)
			smooth->A1[IPar[ic][0]]=A[ic];
		if(rank==2)
			smooth->A2[IPar[ic][0]][IPar[ic][2]]=A[ic];
		if(rank==3)
			smooth->A3[IPar[ic][0]][IPar[ic][2]][IPar[ic][3]]=A[ic];
		if(rank==4)
			smooth->A4[IPar[ic][0]][IPar[ic][2]][IPar[ic][3]][IPar[ic][4]]=A[ic];
		if(rank==5)
			smooth->A5[IPar[ic][0]][IPar[ic][2]][IPar[ic][3]][IPar[ic][4]][IPar[ic][5]]=A[ic];
	}
}

void CSmooth_Linear::Import(CSmooth *smooth){
	unsigned int ic,rank;
	ic=0;
	
	A[0]=smooth->A0;
	ic+=1;
	
	rank=1;
	for(i[0]=0;i[0]<NPars;i[0]++){
		A[ic]=smooth->A1[i[0]];
		ic+=1;
	}
	
	rank=2;
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			A[ic]=smooth->A2[i[0]][i[1]];
			ic+=1;
		}
	}
	
	rank=3;
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				A[ic]=smooth->A3[i[0]][i[1]][i[2]];
				ic+=1;
			}
		}
	}
	
	rank=4;
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					A[ic]=smooth->A4[i[0]][i[1]][i[2]][i[3]];
					ic+=1;
				}
			}
		}
	}
	
	rank=5;
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					for(i[4]=0;i[4]<=i[3];i[4]++){
						A[ic]=smooth->A5[i[0]][i[1]][i[2]][i[3]][i[4]];
						ic+=1;
					}
				}
			}
		}
	}
}

void CSmooth_Linear::SetA_Zero(){
	for(int ic=0;ic<NCoefficients;ic++){
		A[ic]=0;
	}
}

void CSmooth_Linear::SetA_RanGauss(double Amag){
	for(int ic=0;ic<NCoefficients;ic++)
		A[ic]=Amag*randy->ran_gauss();
}

void CSmooth_Linear::SetA_RanSech(double Amag){
	for(int ic=0;ic<NCoefficients;ic++)
		A[ic]=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
}

void CSmooth_Linear::SetA_Constant(double Amag){
	for(int ic=0;ic<NCoefficients;ic++)
		A[ic]=Amag;
}
