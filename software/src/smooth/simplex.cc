#include "msu_smooth/simplex.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>

using namespace std;
using namespace NBandSmooth;
using namespace NMSUPratt;

CSimplexSampler::CSimplexSampler(CparameterMap *parmap){
	randy=new Crandy(123);
	string logfilename=parmap->getS("Simplex_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}
	TrainType=parmap->getI("Simplex_TrainType",1);
	RTrain=parmap->getD("Simplex_RTrain",0.9);
	string prior_info_filename="Info/modelpar_info.txt";
	priorinfo=new CPriorInfo(prior_info_filename);
	CModelParameters::priorinfo=priorinfo;
	ModelDirName=parmap->getS("Simplex_ModelRunDirName","modelruns");
	NPars=priorinfo->NModelPars;
}

void CSimplexSampler::SetThetaSimplex(){
	if(TrainType==1)
		SetThetaType1();
	else if(TrainType==2)
		SetThetaType2();
	else if(TrainType==3)
		SetThetaType3();
	else if(TrainType==4)
		SetThetaType4();
	else{
		CLog::Fatal("Inside CSimplexSampler::SetThetaSimplex, TrainType must be 1,2,3, or 4\n");
	}
	CLog::Info("NTrainingPts="+to_string(NTrainingPts)+"\n");
}

void CSimplexSampler::SetThetaType1(){
	unsigned int ipar,itrain,jtrain;
	double R,z,RTrain;
	RTrain=1.0-1.0/double(NPars+1);
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}

	//double R2;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		//R2=0.0;
		for(ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]*=(RTrain/R);
			//R2+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
		}
		//printf("R2=%g\n",R2);
	}
}

void CSimplexSampler::SetThetaType2(){
	unsigned int ipar,itrain,jtrain,N1,n;
	double R,z,RTrain;
	//RTrain=1.0-1.0/double(NPars+1);
	RTrain=0.95;
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}
	N1=NTrainingPts;
	n=N1;
	NTrainingPts+=N1*(N1-1)/2;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=N1;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
	}
	for(itrain=1;itrain<N1;itrain++){
		for(jtrain=0;jtrain<itrain;jtrain++){
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[n][ipar]=0.5*(ThetaTrain[itrain][ipar]+ThetaTrain[jtrain][ipar]);
			}
			n+=1;
		}
	}
	
	double R1,R2;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		if(itrain<=NPars){
			R1=0.0;
			for(ipar=0;ipar<NPars;ipar++)
				R1+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
			R1=sqrt(R1);
		}
	
		if(itrain>NPars){
			R2=0.0;
			for(ipar=0;ipar<NPars;ipar++)
				R2+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
			R2=sqrt(R2);
		}
		//double Rtest=0.0;
		for(ipar=0;ipar<NPars;ipar++){
			if(itrain<=NPars){
				ThetaTrain[itrain][ipar]*=(0.8*RTrain/R1);
			}
			else{
				ThetaTrain[itrain][ipar]*=(RTrain/R2);
			}
			//Rtest+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
		}
	}
}

void CSimplexSampler::SetThetaType3(){
	// This adds points at half-way points of all connections between pts. of simplex, then moves them outward
	unsigned int ipar,itrain,jtrain,N1,n;
	double R,Rprime,z,RTrain;
	RTrain=0.9;
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}

	N1=NTrainingPts;
	n=N1;
	NTrainingPts+=N1*(N1-1)/2;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=N1;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
	}
	for(itrain=1;itrain<N1;itrain++){
		for(jtrain=0;jtrain<itrain;jtrain++){
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[n][ipar]=0.5*(ThetaTrain[itrain][ipar]+ThetaTrain[jtrain][ipar]);
			}
			n+=1;
		}
	}
	Rprime=R*sqrt(double(NPars-1)/double(2*NPars));

	// Scale
	double R2;
	for(itrain=0;itrain<N1;itrain++){
		R2=0.0;
		for(ipar=0;ipar<=NPars;ipar++){
			ThetaTrain[itrain][ipar]*=0.5*(RTrain/R); //(double(NPars-1)/double(NPars))*(RTrain/R);
			R2+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
		}
		printf("Theta^2=%g\n",R2);
	}
	for(itrain=N1;itrain<NTrainingPts;itrain++){
		R2=0.0;
		for(ipar=0;ipar<=NPars;ipar++){
			ThetaTrain[itrain][ipar]*=(RTrain/Rprime);
			R2+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
		}
		printf("Theta^2=%g\n",R2);
	}
	
}

void CSimplexSampler::SetThetaType4(){
	CLog::Info("WARNING: Using Simplex_TrainType=4 sometimes bombs due to\nencountering non-invertible matrices when solving for coefficients\n");
	unsigned int ipar,itrain,jtrain;
	double R,z,RTrain;

	RTrain=1.0-1.0/double(NPars+1);
	ThetaTrain.resize(2*NPars+3);

	for(itrain=0;itrain<NPars+1;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NPars+1;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}
	for(itrain=0;itrain<NPars+1;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]*=(RTrain/R);
		}
	}


	//make reflection points
	for(itrain=NPars+1;itrain<2*NPars+2;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]=-0.8*ThetaTrain[itrain-NPars-1][ipar]+0.1*randy->ran();
		}
	}
	
	// Put last point at origin
	ThetaTrain[2*NPars+2].resize(NPars);
	for(ipar=0;ipar<NPars;ipar++){
		ThetaTrain[2*NPars+2][ipar]=0.0;
		//printf("ThetaTrain[%u]=%g\n",2*NPars+2,ThetaTrain[2*NPars+2][ipar]);
	}
	NTrainingPts=2*NPars+3;

}

void CSimplexSampler::WriteModelPars(){
	FILE *fptr;
	string filename,dirname,command;
	unsigned int itrain,ipar;
	vector<CModelParameters *> modelparameters(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelparameters[itrain]=new CModelParameters();
		for(ipar=0;ipar<NPars;ipar++){
			modelparameters[itrain]->Theta[ipar]=ThetaTrain[itrain][ipar];
		}
		modelparameters[itrain]->TranslateTheta_to_X();
	}
	for(itrain=0;itrain<NTrainingPts;itrain++){
		//command="rm -r -f "+ModelDirName+"/run"+to_string(itrain)+"/mod_parameters.txt";
		//system(command.c_str());
		dirname=ModelDirName+"/run"+to_string(itrain);
		command="mkdir -p "+dirname;
		system(command.c_str());
		filename=dirname+"/mod_parameters.txt";
		fptr=fopen(filename.c_str(),"w");
		for(ipar=0;ipar<NPars;ipar++){
			fprintf(fptr,"%s %g\n",
			priorinfo->parname[ipar].c_str(),modelparameters[itrain]->X[ipar]);
		}
		fclose(fptr);
	}

}
