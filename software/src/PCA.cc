#include "msu_smooth/PCA.h"

using namespace std;

PCA::PCA(string filename){
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(filename);

	modelruns_dirname=parmap->getS("SmoothEmulator_ModelRunDirName","modelruns");
	string NTrainingStr = parmap->getS("SmoothEmulator_TrainingPts","1");
	
	
	stringstream ss(NTrainingStr);
	string token;

	while(getline(ss, token, ',')) {
		size_t pos = token.find("-");
		if (pos != string::npos) {

			int start = stoi(token.substr(0, pos));
			int end = stoi(token.substr(pos+1));

			for (int i = start; i <= end; i++)
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


void PCA::CalcPCA(){
	vector<vector<double>> YNew;
	vector<double> YBar;
	char filename[300],obs_name[300];
	double y,sigmay;
	int irun;
	FILE *fptr;
		
	for(irun=0;irun<NTrainingList.size();irun++){
		snprintf(filename,300,"%s/run%d/obs.txt",modelruns_dirname.c_str(),NTrainingList[irun]);
		fptr=fopen(filename,"r");
		while(!feof(fptr)){
			fscanf(fptr,"%s %lf %lf",obs_name, &y, &sigmay);
			if(!feof(fptr)){
				Y.at(irun).push_back(y);
				SigmaY.at(irun).push_back(sigmay);
			}

		}
		fclose(fptr);
	}

	YNew.resize(Y.size());

	for(int i = 0 ; i < nruns; i++){

		YNew[i].resize(Y[i].size());

	}

	for (size_t i = 0; i < Y.size(); i++) {
		for (size_t j = 0; j < Y[i].size(); j++) {
			YNew[i][j] = Y[i][j] / 0.5; //SigmaY[i][j] instead of 0.5

		}
	}


	YBar.resize(YNew[0].size());

	for (size_t i = 0; i < YBar.size() ; i++) {
		YBar[i] = 0;
		for (int irun = 0 ; irun < nruns ; irun++){
			YBar[i] += YNew[irun][i] / nruns;

		}

	}

	Eigen::MatrixXd A(YBar.size(), YBar.size());

	for (size_t iy = 0 ; iy < YBar.size() ; iy++) {
		for (size_t ij = 0 ; ij < YBar.size() ; ij++) {
			float avg = 0;
			for (int irun = 0; irun < nruns; irun++) {
				avg += ((YNew[irun][iy] - YBar[iy]) * (YNew[irun][ij] - YBar[ij])) / nruns ;
			}
			A(iy,ij) = avg;

		}
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
	eigvecs = es.eigenvectors();
	Eigen::VectorXd eigvalsVec = es.eigenvalues();

	ofstream eigvalFile("PCA_Info/eigvals.txt");
	if (eigvalFile.is_open())
	{
		eigvalFile << eigvalsVec << '\n';
	}
	eigvalFile.close();

	ofstream eigvecFile("PCA_Info/eigvecs.txt");
	if (eigvecFile.is_open())
	{
		eigvecFile << eigvecs << '\n';
	}
	eigvecFile.close();
}


void PCA::ReadPCA(){
	ifstream eigvalFile("PCA_Info/eigvals.txt");
	string line;
	vector<double> v;

	while (std::getline(eigvalFile, line)){
		v.push_back(stod(line));
	}


	int matrix_size = v.size();

	eigvals.resize(matrix_size, matrix_size);
	eigvals.setZero();



	for (size_t i = matrix_size; i > 0; i--) {
		eigvals(i-1,i-1) = v[matrix_size - i];
	}

	eigvalFile.close();

	ifstream eigvecFile("PCA_Info/eigvecs.txt");
	eigvecs.resize(matrix_size, matrix_size);
	int col = 0;

	while (std::getline(eigvecFile, line)) {
		int row = 0;
		stringstream str_strm;
		str_strm << line;
		string temp_str;
		double temp_dob;
		while(!str_strm.eof()){
			str_strm >> temp_str;
			if (stringstream(temp_str) >> temp_dob) {
				eigvecs(col , matrix_size - 1 - row) = temp_dob;
				row++;
			}
			temp_str = "";
		}
		col++;
	}
}

void PCA::WriteZTraining(){
	int ifile;

	string command="rm -r -f PCA_Info/run*";
	system(command.c_str());

	for (int irun = 0; irun < nruns; irun++) {
		string filename;
		FILE *fptr;

		string dirname= "PCA_Info/run" +to_string(irun);
		command= "mkdir " +dirname;
		system(command.c_str());

		Eigen::VectorXd obsVec(Y[0].size());

		for (size_t i = 0; i < Y[0].size(); i++){
			obsVec(i) = Y[irun][i];
		}

		Eigen::VectorXd PCAVec(Y[0].size());

		PCAVec = eigvecs * obsVec;
		
		ifile=NTrainingList[irun];

		filename=modelruns_dirname+"/run"+to_string(ifile)+"/obs_pca.txt";
		fptr=fopen(filename.c_str(),"w");

		for(size_t i = 0; i < Y[0].size(); i++){
			string obsname = "pca" + to_string(i);
			fprintf(fptr,"%s %s %g\n","dummy" ,obsname.c_str(),PCAVec(i));
		}
		fclose(fptr);
	}
}
