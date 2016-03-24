//============================================================================
// Name        : simulation.cpp
// Author      : ozan
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <array>
#include <numeric>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ctime>
//#include <chrono>

using namespace std;


struct individual
{
	long int ID;
	vector<long int> neuLocus;
	vector<long int> delLocus;
	vector<long int> benLocus;
	double rel_mu;
	double w;
	//int nOff;
};

struct wlist
{
	long int ID;
	double w;
};

struct mutlist
{
	long int mutInd;
	long int mut;
	bool success;
};

bool acompare(wlist lhs, wlist rhs) { return lhs.w > rhs.w; }



int main() {
	clock_t begin=clock();
	const int t_bot = 50;
	const int t_growth = 25;
	const long int ngenes = 4400;
	const long int l = 4640000;
	double mu_d = 0.9;
	double mu_b = 0.001;
	double mu_n = 9.9;
	double shaped = 0.1;
	double scaled = 0.3;
	double betab = 0.01;
	const int nmutator=20;
	double inc_mut_rate=0.1;
	double mean_e_mutator=0.9;
	double selecthighest=0.2;

	bool mutator_control=true;
	bool increased_mutator=true;
	bool del_dep= true;

	//int run=0;



	long int initialN=1;
	//int tbot=1;
	double w_next=1;
	double mu_next=1;


	double mu_neu=mu_n/l;
	double mu_del=mu_d/l;
	double mu_ben=mu_b/l;

	double mu_scaled=1;

	//random_device rd;

	minstd_rand0 gen (clock());
	normal_distribution<double> dist1(mu_neu,0.1*mu_neu);
	normal_distribution<double> dist2(mu_del,0.1*mu_del);
	normal_distribution<double> dist3(mu_ben,0.1*mu_ben);
	gamma_distribution<double> dist4(shaped,scaled);
	exponential_distribution<double> dist5(betab);
	uniform_int_distribution<long int> dist6(0,ngenes-1);
	normal_distribution<double> dist7(mean_e_mutator,0.1*mean_e_mutator);
	poisson_distribution<int> dist8(inc_mut_rate);
	uniform_int_distribution<int> dist9(0,nmutator-1);
	uniform_int_distribution<long int> dist14(0,l-1);
	uniform_int_distribution<int> dist17(0,1000);

	array<array<bool,1000>,1000> success;
	for (int i00=0;i00<500;++i00){
		for (int i01=0;i01<1000;++i01){
			success[i00][i01]=0;
		}
	}
	for (int i00=500;i00<1000;++i00){
		for (int i01=0;i01<1000;++i01){
			success[i00][i01]=1;
		}
	}

	for (int i00=0;i00<500;++i00){
		for (int i01=0;i01<=i00;++i01){
			int place=dist17(gen);
			while (success[i00][place]==1) place=dist17(gen);
			success[i00][place]=1;
		}
	}
	for (int i00=500;i00<1000;++i00){
		for (int i01=0;i01<=1000-i00;++i01){
			int place=dist17(gen);
			while (success[i00][place]==0) place=dist17(gen);
			success[i00][place]=0;
		}
	}

	array<array<double,t_growth>,t_bot> muNeu;
	array<array<double,t_growth>,t_bot> muDel;
	array<array<double,t_growth>,t_bot> muBen;

	for (int i00=0; i00<t_bot ; ++i00){
		for (int i01=0; i01<t_growth; ++i01 ){
			muNeu[i00][i01]=dist1(gen);
			muDel[i00][i01]=dist2(gen);
			muBen[i00][i01]=dist3(gen);
		}
	}

	array<double,ngenes> s_del;
	array<double,ngenes> s_ben;
	for (int i00=0;i00<ngenes;++i00){
		s_del[i00]=dist4(gen);
		s_ben[i00]=dist5(gen);
	}

	array<int,nmutator> mutator;
	array<double,nmutator> mut_e;
	if (mutator_control){
		for (int i00=0;i00<nmutator;++i00){
			mutator[i00]=dist6(gen);
			mut_e[i00]=dist7(gen);
			s_del[mutator[i00]]=0;
			s_ben[mutator[i00]]=0;
		}
	}

	array<array<int,t_growth>,t_bot> result_nOff;
	array<array<double,t_growth>,t_bot> result_wMean;
	array<double,t_bot> result_muScaled;

	vector<int> neuLoci;
	vector<int> delLoci;
	vector<int> benLoci;


	individual initialind={0,{},{},{},1,1};


	vector<int> nOff;
	vector<double> wadj;
	vector<individual> pop;
	long int sumOff;
	long int nhighest;
	long int poolsize; //= new long int;
	long int one=1;
	individual tmpf;
	long int luckyInd;
	vector<wlist> w_list;
	vector<individual> popFinal;

	for (int i1=0;i1<t_bot;++i1){
		clock_t lap0 = clock();
		long int n=initialN;
		pop.clear();
		for (long int i00=0;i00<n;++i00){
			pop.push_back(initialind);
		}
		pop[0].rel_mu=mu_next;
		pop[0].w=w_next;

		if (mutator_control && increased_mutator){
			if (del_dep){
				inc_mut_rate=max(0.0001,1-w_next);
			}
			int tmp=dist8(gen);
			for (int i00=0;i00<tmp;++i00){
				int mut;
				mut=dist9(gen);
				if(!(find(delLoci.begin(),delLoci.end(),mut)!=delLoci.end())){
					pop[0].delLocus.push_back(mut);
					for (int i01=0;i01<nmutator;++i01){
						if(mut==mutator[i01]){
							pop[0].rel_mu=pop[0].rel_mu*mut_e[i01];
						}
					}
				}
			}
		}

		nOff.clear();
		wadj.clear();

		wadj.push_back(pop[0].w*2.5-0.1);
		nOff.push_back(round(wadj[0]));

		if (nOff[0]<0){
		nOff[0]=0;
		}
		if (nOff[0]>2){
		nOff[0]=2;
		}
		sumOff=accumulate(nOff.begin(),nOff.end(),0);

		if (sumOff<=0){
			clock_t end = clock();
			double es = double(end - begin);
			cout << "Extinction at bottleneck" << endl;
			cout << "t_bot= " << i1 << endl;
			cout << "time= " << es << endl;
			exit(0);
		}

		result_muScaled[i1]=mu_scaled;

		for (int i2=0;i2<t_growth;++i2){
			clock_t lap = clock();
			for (long int i00=0;i00<n;++i00){
				if (nOff[i00]>0){
					for (int i01=0;i01<nOff[i00]-1;++i01){
						pop.push_back(pop[i00]);
					}
				}
				if (nOff[i00]==0){
					pop.erase(pop.begin()+i00);
				}
			}

			n=sumOff;

			for (long int i00=0;i00<n;++i00){
				pop[i00].ID=i00;
			}

			{
			poisson_distribution<long int> dist10(n*l*muNeu[i1][i2]);
			poisson_distribution<long int> dist11(n*l*muDel[i1][i2]);
			poisson_distribution<long int> dist12(n*l*muBen[i1][i2]);

			long int nNeu=dist10(gen);
			long int nDel=dist11(gen);
			long int nBen=dist12(gen);

			uniform_int_distribution<long int> dist13(0,n-1);
			{
				{
				mutlist neu[100000];
				while (nNeu>100000){
					int rmu;
					for (long int i00=0;i00<100000;++i00){
						neu[i00].mutInd=dist13(gen);
						neu[i00].mut=dist6(gen);
						//binomial_distribution<long int> dist15(1,pop[del[i00].mutInd].rel_mu);
						//del[i00].success=dist15(gen);
						rmu=1000*pop[neu[i00].mutInd].rel_mu;
						rmu=max(0,min(1000,rmu));
						neu[i00].success=success[rmu][dist17(gen)];
					}
					for (long int i00=0;i00<100000;++i00){
						if (neu[i00].success){
							pop[neu[i00].mutInd].neuLocus.push_back(neu[i00].mut);
						}
					}
					nNeu-=100000;
				}
				}
				mutlist neu[nNeu];
				int rmu;
				for (long int i00=0;i00<nNeu;++i00){
					neu[i00].mutInd=dist13(gen);
					neu[i00].mut=dist14(gen);
					//binomial_distribution<long int> dist15(1,pop[neu[i00].mutInd].rel_mu);
					rmu=1000*pop[neu[i00].mutInd].rel_mu;
					rmu=max(0,min(1000,rmu));
					neu[i00].success=success[rmu][dist17(gen)];
				}
				for (long int i00=0;i00<nNeu;++i00){
					if (neu[i00].success){
						pop[neu[i00].mutInd].neuLocus.push_back(neu[i00].mut);
					}
				}
			}

			{
				{
				mutlist del[100000];
				while (nDel>100000){
					int rmu;
					for (long int i00=0;i00<100000;++i00){
						del[i00].mutInd=dist13(gen);
						del[i00].mut=dist6(gen);
						//binomial_distribution<long int> dist15(1,pop[del[i00].mutInd].rel_mu);
						//del[i00].success=dist15(gen);
						rmu=1000*pop[del[i00].mutInd].rel_mu;
						rmu=max(0,min(1000,rmu));
						del[i00].success=success[rmu][dist17(gen)];
					}
					for (long int i00=0;i00<100000;++i00){
						if (del[i00].success && !(find(delLoci.begin(),delLoci.end(),del[i00].mut)!=delLoci.end())){
							pop[del[i00].mutInd].delLocus.push_back(del[i00].mut);
							pop[del[i00].mutInd].w-=s_del[del[i00].mut];
							for (int i01=0;i01<nmutator;++i01){
								if(del[i00].mut==mutator[i01]){
									(pop[del[i00].mutInd].rel_mu)*=mut_e[i01];
								}
							}
						}
					}
					nDel-=100000;
				}
				}
				mutlist del[nDel];
				int rmu;
				for (long int i00=0;i00<nDel;++i00){
					del[i00].mutInd=dist13(gen);
					del[i00].mut=dist6(gen);
					//binomial_distribution<long int> dist15(1,pop[del[i00].mutInd].rel_mu);
					//del[i00].success=dist15(gen);
					rmu=1000*pop[del[i00].mutInd].rel_mu;
					rmu=max(0,min(1000,rmu));
					del[i00].success=success[rmu][dist17(gen)];
				}
				for (long int i00=0;i00<nDel;++i00){
					if (del[i00].success && !(find(delLoci.begin(),delLoci.end(),del[i00].mut)!=delLoci.end())){
						pop[del[i00].mutInd].delLocus.push_back(del[i00].mut);
						pop[del[i00].mutInd].w-=s_del[del[i00].mut];
						for (int i01=0;i01<nmutator;++i01){
							if(del[i00].mut==mutator[i01]){
								(pop[del[i00].mutInd].rel_mu)*=mut_e[i01];
							}
						}
					}
				}
			}

			{
			mutlist ben[nBen];
			int rmu;
			for (long int i00=0;i00<nBen;++i00){
				ben[i00].mutInd=dist13(gen);
				ben[i00].mut=dist6(gen);
				//binomial_distribution<long int> dist15(1,pop[ben[i00].mutInd].rel_mu);
				//ben[i00].success=dist15(gen);
				rmu=1000*pop[ben[i00].mutInd].rel_mu;
				rmu=max(0,min(1000,rmu));
				ben[i00].success=success[rmu][dist17(gen)];
			}
			for (long int i00=0;i00<nBen;++i00){
				if (ben[i00].success && !(find(benLoci.begin(),benLoci.end(),ben[i00].mut)!=benLoci.end())){
					pop[ben[i00].mutInd].benLocus.push_back(ben[i00].mut);
					pop[ben[i00].mutInd].w+=s_ben[ben[i00].mut];
					/*
					for (int i01=0;i01<nmutator;++i01){
						if(ben[i00].mut==mutator[i01]){
							pop[ben[i00].mutInd].rel_mu/=mut_e[i01];
						}
					}*/
				}
			}
			}
			}

			wadj.clear();
			nOff.clear();
			{
				double wtmp=0;
				for (int i00=0;i00<n;++i00){
					double tmp0=pop[i00].w*2.5-0.1;
					wadj.push_back(tmp0);
					wtmp+=pop[i00].w;
					if (wadj[i00]<=0){
						wadj[i00]=0.000001;
					}
					int tmp1=round(wadj[i00]);
					nOff.push_back(tmp1);

					if (nOff[i00]<0){
					nOff[i00]=0;
					}
					if (nOff[i00]>2){
					nOff[i00]=2;
					}
				}

				sumOff=accumulate(nOff.begin(),nOff.end(),0);

				result_nOff[i1][i2]=sumOff;
				result_wMean[i1][i2]=wtmp/n;
			}
			if (sumOff<=0){
				cout << "Extinction at growth" << endl;
				cout << "t_bot= " << i1 << endl;
				cout << "t_growth= " << i2 << endl;

				clock_t end = clock();
				double es = double(end - begin);
				cout << "time= " << es << endl;
				exit(0);
			}
			if(n>1000000){

				clock_t end = clock();
				double es = double(end - lap);
				cout << "Finished Gen " << i2 << " in Bot " << i1 << " n= " << n << " time= " << es << endl;
				cout << "Population over a million" << endl;
				break;
			}
			{

			clock_t end = clock();
			double es = double(end - lap);
			cout << "Finished Gen " << i2 << " in Bot " << i1 << " n= " << n << " time= " << es << endl;
			}
		}


		w_list.clear();
		{
			wlist ind;
			for (long int i00=0;i00<n;++i00){
				ind.ID=pop[i00].ID;
				ind.w=pop[i00].w;
				w_list.push_back(ind);
			}


		sort(w_list.begin(),w_list.end(),acompare);


		nhighest = round(n*selecthighest);
		poolsize = max(one,nhighest);

		uniform_int_distribution<long int> dist18(0,poolsize-1);


			long int tmp = dist18(gen);
			luckyInd = w_list[tmp].ID;

		}
		tmpf=pop[luckyInd];

		for (unsigned long int i00=0;i00<tmpf.neuLocus.size();++i00){
			neuLoci.push_back(tmpf.neuLocus[i00]);
		}
		for (unsigned long int i00=0;i00<tmpf.delLocus.size();++i00){
			delLoci.push_back(tmpf.delLocus[i00]);
		}
		for (unsigned long int i00=0;i00<tmpf.benLocus.size();++i00){
			benLoci.push_back(tmpf.benLocus[i00]);
		}
		w_next=tmpf.w;
		mu_next=tmpf.rel_mu;

		//tbot=i1+1;
		popFinal.push_back(tmpf);

		clock_t end = clock();
		double es = double(end - lap0);
		cout << "finished Bot " << i1 << " time= " << es << endl;
	}

	clock_t end = clock();
	double es = double(end - begin);
	cout << "End of simulation" << " time= " << es << endl;

	return 0;
}



