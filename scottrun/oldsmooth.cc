#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "/Users/scottpratt/git/commonutils/software/include/randy.h"
#include "/Users/scottpratt/git/commonutils/software/src/NumMath/randy.cc"

using namespace std;

int main(int argc,char *argv[]){
	vector<vector<vector<vector<vector<double>>>>> A5;
	int Nrank=5,Npars=8,j,isame;
	vector<int> i(Nrank);
	vector<int> factorial(Npars+1);
	double R=1.0,Y,factfact,dY;
	vector<double> theta(Nrank);
	vector<int> count_same(5);
	CRandy *randy=new CRandy(1234);
	
	// set factorial
	factorial[0]=factorial[1]=1.0;
	for(j=2;j<=Npars;j++)
		factorial[j]=double(j)*factorial[j-1];
	
	// A5et Model Coefficients A5
	A5.resize(Npars);
	for(i[0]=0;i[0]<Npars;i[0]++){
		A5[i[0]].resize(i[0]+1);
		for(i[1]=0;i[1]<=i[0];i[1]++){
			A5[i[0]][i[1]].resize(i[1]+1);
			for(i[2]=0;i[2]<=i[1];i[2]++){
				A5[i[0]][i[1]][i[2]].resize(i[2]+1);
				for(i[3]=0;i[3]<=i[2];i[3]++){
					A5[i[0]][i[1]][i[2]][i[3]].resize(i[3]+1);
					for(i[4]=0;i[4]<=i[3];i[4]++){
						//A5[i[0]][i[1]][i[2]][i[3]][i[4]]=randy->ran_gauss();
						A5[i[0]][i[1]][i[2]][i[3]][i[4]]=1.0;
					}
				}
			}
		}
	}
	
	// Set theta
	for(j=0;j<Npars;j++){
		//theta[j]=randy->ran_gauss()/R;
		theta[j]=1.0;
	}
	
	// Calculate Y(theta)
	Y=0.0;
	int nterms=0;
	for(i[0]=0;i[0]<Npars;i[0]++){
		for(i[1]=0;i[1]<=i[0];i[1]++){
			for(i[2]=0;i[2]<=i[1];i[2]++){
				for(i[3]=0;i[3]<=i[2];i[3]++){
					for(i[4]=0;i[4]<=i[3];i[4]++){
						// get duplication factor = factfact
						for(j=1;j<Nrank;j++)
							count_same[j]=0;
						isame=0;
						count_same[isame]=1;
						for(j=1;j<Nrank;j++){
							if(i[j]!=i[j-1])
								isame+=1;
							count_same[isame]+=1;
						}
						factfact=factorial[Nrank];
						for(j=0;j<Nrank;j++)
							factfact=factfact/factorial[count_same[j]];
						dY=factfact;
						printf("%d %d %d %d %d: %g, Y=%g\n",i[0],i[1],i[2],i[3],i[4],factfact,Y);
						for(j=0;j<Nrank;j++){
							dY*=theta[i[j]];
						}
						Y+=dY;	
						nterms+=1;					
					}
				}
			}
		}
	}
	printf("nterms=%d\n",nterms);
	printf("Y=%g\n",Y);
	// Calculate Y(theta)
	Y=0.0;
	for(i[0]=0;i[0]<Npars;i[0]++){
		for(i[1]=0;i[1]<Npars;i[1]++){
			for(i[2]=0;i[2]<Npars;i[2]++){
				for(i[3]=0;i[3]<Npars;i[3]++){
					for(i[4]=0;i[4]<Npars;i[4]++){
						// get duplication factor = factfact
						dY=1.0;
						Y+=dY;						
					}
				}
			}
		}
	}
	printf("Y=%g\n",Y);
}

//Note that if I have n1 of the same theta and n2 of the same theta,,  s.t. n1+n2+n3.. = Nrank, then the number of combinations is Nrank!/(n1!n2!n3!...). Because we are dividing Nrank!, our counting factor is 1/(n1!n2!n3!...)