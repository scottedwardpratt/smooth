//#include "emulator.h"
#include "smooth.h"
#include "gslmatrix.h"
using namespace std;
vector<unsigned int> CSmooth::factorial={};

CSmooth::CSmooth(){
	//
}

CSmooth::CSmooth(CparameterMap *parmap){
	NPars=parmap->getI("Smooth_NPars",0);
	UseRFactor=parmap->getB("Smooth_UseRFactor",false);
	unsigned int ic,j,isame,ir;
	vector<unsigned int> countsame;
	vector<unsigned int> dummy;
	vector<unsigned int> i;
	MaxRank=5;
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
	rank.resize(ic+1);
	rank[ic]=0;
	dupfactor.push_back(1.0);
	IPar.resize(ic+1);
	IPar[0].resize(rank[ic]);
	
	ic+=1;
	for(i[0]=0;i[0]<NPars;i[0]++){
		dupfactor.push_back(0);
		IPar.resize(ic+1);
		rank.resize(ic+1);
		rank[ic]=1;
		IPar[ic].resize(rank[ic]);
		for(ir=0;ir<rank[ic];ir++){
			IPar[ic][ir]=i[ir];
		}
		for(j=0;j<rank[ic];j++)
			countsame[j]=0;
		isame=0;
		countsame[isame]=1;
		for(j=1;j<rank[ic];j++){
			if(i[j]!=i[j-1])
				isame+=1;
			countsame[isame]+=1;
		}
		dupfactor[ic]=factorial[rank[ic]];
		for(j=0;j<rank[ic];j++)
			dupfactor[ic]/=factorial[countsame[j]];
		ic+=1;
	}
	
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			dupfactor.push_back(0);
			IPar.resize(ic+1);
			rank.resize(ic+1);
			rank[ic]=2;
			IPar[ic].resize(rank[ic]);
			for(ir=0;ir<rank[ic];ir++){
				IPar[ic][ir]=i[ir];
			}
			for(j=0;j<rank[ic];j++)
				countsame[j]=0;
			isame=0;
			countsame[isame]=1;
			for(j=1;j<rank[ic];j++){
				if(i[j]!=i[j-1])
					isame+=1;
				countsame[isame]+=1;
			}
			dupfactor[ic]=factorial[rank[ic]];
			for(j=0;j<rank[ic];j++)
				dupfactor[ic]/=factorial[countsame[j]];
			ic+=1;
		}
	}
	
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				dupfactor.push_back(0);
				IPar.resize(ic+1);
				rank.resize(ic+1);
				rank[ic]=3;
				IPar[ic].resize(rank[ic]);
				for(ir=0;ir<rank[ic];ir++){
					IPar[ic][ir]=i[ir];
				}
				for(j=0;j<rank[ic];j++)
					countsame[j]=0;
				isame=0;
				countsame[isame]=1;
				for(j=1;j<rank[ic];j++){
					if(i[j]!=i[j-1])
						isame+=1;
					countsame[isame]+=1;
				}
				dupfactor[ic]=factorial[rank[ic]];
				for(j=0;j<rank[ic];j++)
					dupfactor[ic]/=factorial[countsame[j]];
				ic+=1;
			}
		}
	}
	
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					dupfactor.push_back(0);
					IPar.resize(ic+1);
					rank.resize(ic+1);
					rank[ic]=4;
					IPar[ic].resize(rank[ic]);
					for(ir=0;ir<rank[ic];ir++){
						IPar[ic][ir]=i[ir];
					}
					for(j=0;j<rank[ic];j++)
						countsame[j]=0;
					isame=0;
					countsame[isame]=1;
					for(j=1;j<rank[ic];j++){
						if(i[j]!=i[j-1])
							isame+=1;
						countsame[isame]+=1;
					}
					dupfactor[ic]=factorial[rank[ic]];
					for(j=0;j<rank[ic];j++)
						dupfactor[ic]/=factorial[countsame[j]];
					ic+=1;
				}
			}
		}
	}
	
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					for(i[4]=0;i[4]<=i[3];i[4]++){
						dupfactor.push_back(0);
						IPar.resize(ic+1);
						rank.resize(ic+1);
						rank[ic]=5;
						IPar[ic].resize(rank[ic]);
						for(ir=0;ir<rank[ic];ir++){
							IPar[ic][ir]=i[ir];
						}
						for(j=0;j<rank[ic];j++)
							countsame[j]=0;
						isame=0;
						countsame[isame]=1;
						for(j=1;j<rank[ic];j++){
							if(i[j]!=i[j-1])
								isame+=1;
							countsame[isame]+=1;
						}
						dupfactor[ic]=factorial[rank[ic]];
						for(j=0;j<rank[ic];j++)
							dupfactor[ic]/=factorial[countsame[j]];
						ic+=1;
					}
				}
			}
		}
	}
	NCoefficients=ic;
	if(NCoefficients!=IPar.size()){
		printf("size mismatch\n");
		exit(1);
	}
	printf("NCoefficients=%u\n",NCoefficients);
}

double CSmooth::CalcY(vector<double> &A,double LAMBDA,vector<double> &theta){
	unsigned int ic,ir;
	double answer=0.0,term,rfactor;
	rfactor=GetRFactor(LAMBDA,theta);
	answer=0.0;
	for(ic=0;ic<NCoefficients;ic++){
		term=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/LAMBDA;
		}
		//printf("ic=%3u, A=%10.7f, term=%10.7f dupfactor=%10.7u\n",ic,A[ic],term,dupfactor[ic]);
		answer+=term;
	}
	answer*=rfactor;
	return answer;
}

double CSmooth::CalcY_Remainder(vector<double> &A,double LAMBDA,vector<double> &theta,unsigned int NTrainingPts){
	unsigned int ic,ir;
	double answer=0.0,term,rfactor;
	rfactor=GetRFactor(LAMBDA,theta);
	answer=0.0;
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		term=A[ic]*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/LAMBDA;
		}
		answer+=term;
	}
	answer*=rfactor;
	return answer;
}

double CSmooth::GetRFactor(double LAMBDA,vector<double> &theta){
	unsigned int ir,ipar,NPars=theta.size();
	double r2=0.0,answer;
	if(UseRFactor){
		for(ipar=0;ipar<NPars;ipar++)
			r2+=theta[ipar]*theta[ipar];
		answer=1.0;
		for(ir=1;ir<=MaxRank;ir++){
			answer+=pow(r2/(LAMBDA*LAMBDA),ir)/double(factorial[ir]);
		}
		answer=1.0/sqrt(answer);
	}
	else{
		answer=1.0;
	}
	return answer;

}

double CSmooth::GetM(int ic,double LAMBDA,vector<double> &theta){
	unsigned int ir;
	double answer,rfactor=1.0;
	rfactor=GetRFactor(LAMBDA,theta);
	answer=rfactor*sqrt(double(dupfactor[ic])/double(factorial[rank[ic]]));
	for(ir=0;ir<rank[ic];ir++){
		answer*=theta[IPar[ic][ir]]/LAMBDA;
	}
	return answer;
}

void CSmooth::Copy(CSmooth *smooth){
	MaxRank=smooth->MaxRank;
	NPars=smooth->NPars;
	NCoefficients=smooth->NCoefficients;
	dupfactor=smooth->dupfactor;
	IPar=smooth->IPar;
	rank=smooth->rank;
	factorial=smooth->factorial;
}


