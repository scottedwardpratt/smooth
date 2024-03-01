#include "msu_smooth/PCA.h"

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

PCA::PCA(CparameterMap *parmapset){
	parmap=parmapset;
	string command="mkdir -p PCA_Info";
	system(command.c_str());
	
	observable_info=new CObservableInfo("Info/observable_info.txt");
	observable_info->ReadExperimentalInfo("Info/experimental_info.txt");
	
	modelruns_dirname=parmap->getS("SmoothEmulator_ModelRunDirName","modelruns");
	string NTrainingStr = parmap->getS("SmoothEmulator_TrainingPts","1");
	
	stringstream ss(NTrainingStr);
	string token;

	while(getline(ss, token, ',')) {
		size_t pos = token.find("-");
		if (pos != string::npos){

			unsigned int start = stoi(token.substr(0, pos));
			unsigned int end = stoi(token.substr(pos+1));

			for (unsigned int i = start; i <= end; i++)
				NTrainingList.push_back(i);
		}
		else {
			NTrainingList.push_back(stoi(token));
		}
	}
	nruns=NTrainingList.size();

	Y.resize(nruns);
	SigmaY.resize(nruns);
	
}

void PCA::CalcTransformationInfo(){
	vector<vector<double>> Ytilde,Y;
	vector<double> Ytildebar;
	char filename[300],obs_name[300];
	double y,sigmay;
	unsigned int iy,jy,irun,nruns=NTrainingList.size();
	Nobs=observable_info->NObservables;
	Eigen::MatrixXd *A;
	FILE *fptr;
	Y.resize(nruns);
	Ytilde.resize(nruns);
	for(irun=0;irun<nruns;irun++){
		Y[irun].resize(Nobs);
		Ytilde[irun].resize(Nobs);
	}
	Ybar.resize(Nobs);
	SigmaY.resize(Nobs);
	for(iy=0;iy<Nobs;iy++){
		Ybar[iy]=0.0;
		SigmaY[iy]=observable_info->SigmaExp[iy];
	}
	
	for(irun=0;irun<nruns;irun++){
		snprintf(filename,300,"%s/run%u/obs.txt",modelruns_dirname.c_str(),NTrainingList[irun]);
		fptr=fopen(filename,"r");
		while(!feof(fptr)){
			fscanf(fptr,"%s %lf %lf",obs_name, &y, &sigmay);
			if(!feof(fptr)){
				iy=observable_info->GetIPosition(obs_name);
				Y[irun][iy]=y;
				Ybar[iy]+=y/double(nruns);
			}
		}
		fclose(fptr);
	}
	
	for(iy=0;iy<Nobs;iy++){
		for(irun=0;irun<nruns;irun++){
			Ytilde[irun][iy]=(Y[irun][iy]-Ybar[iy])/SigmaY[iy];
		}
	}
	
	// Calculate Covariance Matrix
	
	A = new Eigen::MatrixXd(Nobs,Nobs);
	A->setZero();
	for(iy=0;iy<Nobs;iy++){
		for(jy=0;jy<Nobs;jy++){
			for(irun=0;irun<nruns;irun++){
				(*A)(iy,jy)+=(Ytilde[irun][iy]*Ytilde[irun][jy])/double(nruns);
			}
		}
	}
				
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*A);
	eigvals = es.eigenvalues();
	eigvecs = es.eigenvectors();
	delete A;
	
	
	// Write Transformation Info
	ofstream PCAInfoFile("PCA_Info/transformation_info.txt");
	
	if (PCAInfoFile.is_open()){
		PCAInfoFile << "Nobs  Nruns" << endl;
		PCAInfoFile << Nobs <<  " " << nruns << endl;
		PCAInfoFile << "<Y>" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << Ybar[iy] << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "SigmaY" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << SigmaY[iy] << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "EigenValues" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << eigvals(iy) << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "EigenVectors" << endl;
		for(iy=0;iy<Nobs;iy++){
			for(jy=0;jy<Nobs;jy++){
				PCAInfoFile << eigvecs(iy,jy) << " ";
			}
			PCAInfoFile << endl;
		}
	}
	PCAInfoFile.close();

	// Write PCA Observables for Training Pts
	
	vector<vector<double>> Z;
	unsigned int iz;
	string pcaname;
	
	Z.resize(nruns);
	for(irun=0;irun<nruns;irun++){
		Z[irun].resize(nruns);
	}
	
	for(irun=0;irun<nruns;irun++){
		for(iz=0;iz<Nobs;iz++){
			Z[irun][iz]=0.0;
			for(iy=0;iy<Nobs;iy++){
				Z[irun][iz]+=eigvecs(iz,iy)*Ytilde[irun][iy];
			}
		}
		snprintf(filename,300,"%s/run%u/obs_pca.txt",modelruns_dirname.c_str(),NTrainingList[irun]);
		fptr=fopen(filename,"w");
		for(iz=0;iz<Nobs;iz++){
			pcaname="z"+to_string(iz);
			fprintf(fptr,"%s  %g  0.0\n",pcaname.c_str(),Z[irun][iz]);	
		}
		fclose(fptr);
	}
	
	// Write Info File for PCA
	
	double SA0Zsquared,SA0Y,SA0Z;
	fptr=fopen("PCA_Info/observable_info.txt","w");
	for(iz=0;iz<Nobs;iz++){
		SA0Zsquared=0.0;
		for(iy=0;iy<Nobs;iy++){
			SA0Y=observable_info->SigmaA0[iy];
			SA0Zsquared+=SA0Y*SA0Y*eigvecs(iz,iy)*eigvecs(iz,iy);
		}
		SA0Z=sqrt(SA0Zsquared);
		pcaname="z"+to_string(iz);
		fprintf(fptr,"%s %g\n",pcaname.c_str(),SA0Z);
	}
	fclose(fptr);
	
	// write Experimental File for PCA
	vector<double> YExpTilde,ZExp;
	YExpTilde.resize(Nobs);
	ZExp.resize(Nobs);
	fptr=fopen("Info/experimental_info_pca.txt","w");
	for(iy=0;iy<Nobs;iy++){
		YExpTilde[iy]=(observable_info->YExp[iy]-Ybar[iy])/SigmaY[iy];
	}
	for(iz=0;iz<Nobs;iz++){
		ZExp[iz]=0.0;
		for(iy=0;iy<Nobs;iy++){
			ZExp[iz]+=eigvecs(iz,iy)*YExpTilde[iy];
		}
		pcaname="z"+to_string(iz);
		fprintf(fptr,"%s %g %g\n",pcaname.c_str(),ZExp[iz],1.0);
	}
	
	fclose(fptr);
	
}

void PCA::ReadTransformationInfo(){
	char dummy[2000];
	unsigned int iy,jy;
	FILE *fptr=fopen("PCA_Info/transformation_info.txt","r");
	fgets(dummy,2000,fptr);
	fscanf(fptr,"%u %u",&Nobs,&nruns);
	Ybar.resize(Nobs);
	SigmaY.resize(Nobs);
	eigvals.resize(Nobs);
	eigvecs.resize(Nobs,Nobs);
	fscanf(fptr,"%s",dummy);
	for(iy=0;iy<Nobs;iy++){
		fscanf(fptr,"%lf ",&Ybar[iy]);
	}
	fscanf(fptr,"%s",dummy);
	for(iy=0;iy<Nobs;iy++){
		fscanf(fptr,"%lf ",&SigmaY[iy]);
	}
	fscanf(fptr,"%s",dummy);
	for(iy=0;iy<Nobs;iy++){
		fscanf(fptr,"%lf",&eigvals(iy));
	}
	fscanf(fptr,"%s",dummy);
	for(iy=0;iy<Nobs;iy++){
		for(jy=0;jy<Nobs;jy++){
			fscanf(fptr,"%lf ",&eigvecs(iy,jy));
		}
	}
	fclose(fptr);
}

void PCA::TransformZtoY(vector<double> &Z,vector<double> &SigmaZ_emulator,
vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int iy,ky;
	SigmaY_emulator.resize(Nobs);
	for(iy=0;iy<Nobs;iy++){
		SigmaY_emulator.resize(Nobs);
	}
	for(iy=0;iy<Nobs;iy++){
		Y[iy]=0.0;
		for(ky=0;ky<Nobs;ky++){
			Y[iy]+=Z[ky]*eigvecs(ky,iy);
		}
		SigmaY_emulator[iy]=0.0;
		for(ky=0;ky<Nobs;ky++){
			SigmaY_emulator[iy]+=SigmaZ_emulator[ky]*SigmaZ_emulator[ky]*eigvecs(ky,iy)*eigvecs(ky,iy);
		}
		SigmaY_emulator[iy]=SigmaY[iy]*sqrt(SigmaY_emulator[iy]);
	}
}


void PCA::TransformYtoZ(vector<double> &Z,vector<double> &SigmaZ_emulator,
vector<double> &Y,vector<double> &SigmaY_emulator){
	unsigned int iy,ky;
	SigmaZ_emulator.resize(Nobs);
	for(iy=0;iy<Nobs;iy++){
		SigmaZ_emulator.resize(Nobs);
	}
	for(iy=0;iy<Nobs;iy++){
		Z[iy]=0.0;
		for(ky=0;ky<Nobs;ky++){
			Z[iy]+=Y[ky]*eigvecs(iy,ky);
		}
		SigmaZ_emulator[iy]=0.0;
		for(ky=0;ky<Nobs;ky++){
			SigmaZ_emulator[iy]+=(SigmaY_emulator[ky]*SigmaY_emulator[ky]/(SigmaY[iy]*SigmaY[iy]))*eigvecs(iy,ky)*eigvecs(iy,ky);
		}
		SigmaZ_emulator[iy]=sqrt(SigmaZ_emulator[iy]);
	}
}


