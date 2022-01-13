//#include "emulator.h"
#include "smooth.h"
#include "gslmatrix.h"
using namespace std;
CRandy *CSmooth::randy=NULL;
vector<unsigned int> CSmooth::factorial={};

CSmooth::CSmooth(){
	//
}

CSmooth::CSmooth(unsigned int MaxRank_set,unsigned int NPars_set){
	unsigned int ic,j,isame,ir;
	vector<unsigned int> countsame;
	vector<unsigned int> dummy;
	vector<unsigned int> i;
	NPars=NPars_set;
	MaxRank=MaxRank_set;
	Lambda.resize(NPars);
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
	A.push_back(0);
	dupfactor.push_back(1.0);
	IPar.resize(ic+1);
	IPar[0].resize(rank[ic]);
	
	ic+=1;
	for(i[0]=0;i[0]<NPars;i[0]++){
		A.push_back(0);
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
	printf("Thus far, ic=%u\n",ic);
	
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			A.push_back(0);
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
	printf("Thus far, ic=%u\n",ic);

	
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				A.push_back(0);
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
	printf("Thus far, ic=%u\n",ic);

	
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					A.push_back(0);
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
	printf("Thus far, ic=%u\n",ic);

	
	for(i[0]=0;i[0]<NPars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					for(i[4]=0;i[4]<=i[3];i[4]++){
						A.push_back(0);
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
	if(NCoefficients!=A.size()){
		printf("size mismatch\n");
		exit(1);
	}
}

double CSmooth::CalcY(vector<double> &theta){
	unsigned int ic,ir;
	double answer=0.0,term;
	answer=0.0;
	for(ic=0;ic<NCoefficients;ic++){
		term=A[ic]*dupfactor[ic]/factorial[rank[ic]];
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/Lambda[IPar[ic][ir]];
		}
		answer+=term;
	}
	return answer;
}

void CSmooth::SetLambda(double lambdaset){
	for(unsigned int ipar=0;ipar<NPars;ipar++)
		Lambda[ipar]=lambdaset;
}

void CSmooth::SetA_Zero(){
	for(int ic=0;ic<NCoefficients;ic++){
		A[ic]=0;
	}
}

void CSmooth::SetA_RanGauss(double Amag){
	for(int ic=0;ic<NCoefficients;ic++)
		A[ic]=Amag*randy->ran_gauss();
}

void CSmooth::SetA_RanSech(double Amag){
	double r=randy->ran();
	for(int ic=0;ic<NCoefficients;ic++)
		A[ic]=Amag*2.0*atanh(tan(0.5*PI*(r-0.5)));
}

void CSmooth::SetA_Constant(double Amag){
	for(int ic=0;ic<NCoefficients;ic++)
		A[ic]=Amag;
}

// This adjust first NTrainingPts coefficients to reproduce Yvalues
void CSmooth::CalcAFromTraining(vector<vector<double>> &thetatheta,vector<double> YTrain){
	unsigned int itrain,ic,ir;
	unsigned int NTrainingPts=YTrain.size();
	vector<vector<double>> M;
	if(thetatheta.size()!=NTrainingPts){
		printf("CSmooth:: array size mismatch!!\n");
		exit(1);
	}
	vector<double> YTarget;
	YTarget.resize(NTrainingPts);
	M.resize(NTrainingPts);
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		YTarget[itrain]=YTrain[itrain]-CalcY_Remainder(thetatheta[itrain],NTrainingPts);
		M[itrain].resize(NTrainingPts);
		for(ic=0;ic<NTrainingPts;ic++){
			M[itrain][ic]=0.0;
			for(ic=0;ic<NTrainingPts;ic++){
				M[itrain][ic]=dupfactor[ic]/double(factorial[rank[ic]]);
				for(ir=0;ir<rank[ic];ir++){
					M[itrain][ic]*=thetatheta[itrain][IPar[ic][ir]]/Lambda[IPar[ic][ir]];
				}
			}
		}
	}
	CGSLMatrix_Real *gslmatrix=new CGSLMatrix_Real(NTrainingPts);
	
	gslmatrix->SolveLinearEqs(YTarget,M,A);
	delete gslmatrix;
}

double CSmooth::CalcY_Remainder(vector<double> &theta,unsigned int NTrainingPts){
	unsigned int ic,ir;
	double answer=0.0,term;
	answer=0.0;
	for(ic=NTrainingPts;ic<NCoefficients;ic++){
		term=A[ic]*dupfactor[ic]/factorial[rank[ic]];
		for(ir=0;ir<rank[ic];ir++){
			term*=theta[IPar[ic][ir]]/Lambda[IPar[ic][ir]];
		}
		answer+=term;
	}
	return answer;
}

void CSmooth::Copy(CSmooth *smooth){
	MaxRank=smooth->MaxRank;
	NPars=smooth->NPars;
	NCoefficients=smooth->NCoefficients;
	Lambda=smooth->Lambda;
	A=smooth->A;
	dupfactor=smooth->dupfactor;
	IPar=smooth->IPar;
	rank=smooth->rank;
	factorial=smooth->factorial;
	randy=smooth->randy;
}

double CSmooth::GetLog_AProb(double Amag){
	double answer=0.0;
	for(int ic=0;ic<NCoefficients;ic++){
		answer-=log(1.0+A[ic]*A[ic]/(Amag*Amag));
	}
	return answer;
}

