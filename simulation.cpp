//============================================================================
// Name        : simulation.cpp
// Author      : ozan kiratli - ozankiratli@gmail.com
// Version     : 0.1
// Copyright   : Contact author before scientific use.
// Description : Individual based haploid population simulation with bottlenecking feature
//============================================================================

#include <iostream>
#include <array>
#include <numeric>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <string>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "variables.h"

using namespace std;

struct individual
{
	long int ID;
	vector<long int> neuLocus;
	vector<long int> delLocus;
	vector<long int> benLocus;
	double rel_mu;
	double w;
	int nOff;
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

struct resultsum
{
	long int nNeut;
	double rNeut;
	long int nDele;
	double rDele;
	long int nBene;
	double rBene;
	long int nMuta;
	double muScale;
	double w;
};

string fileformat=".out";

bool acompare(wlist lhs, wlist rhs) { return lhs.w > rhs.w; }


const int t_bot = T_bot;
const int t_growth = T_growth;
const long int ngenes = Ngenes;
const long int l = L;
double mu_d = Mu_d;
double mu_b = Mu_b;
double mu_n = Mu_n;
double shaped = Shaped;
double scaled = Scaled;
double betab = Betab;
const int nmutator = Nmutator;
double inc_mut_rate = Inc_mut_rate;
double mean_e_mutator = Mean_e_mutator;
double selecthighest = Selecthighest;

bool mutator_control = Mutator_control;
bool increased_mutator = Increased_mutator;
bool del_dep = Del_dep;

int run = Run;


void generate_results(vector<individual> popFinal,
		array<array<int,t_growth>,t_bot> result_nOff,
		array<array<double,t_growth>,t_bot> result_wMean,
		array<double,t_bot> result_muScaled,
		array<int,nmutator> mutator,
		array<double,ngenes> s_del,
		array<double,ngenes> s_ben,
		array<double,nmutator> mut_e,
		string fileheader,
		int tbot){

	string filename;
	vector< vector<long int>> result_genome;
	vector< vector<long int>> result_deleterious;
	vector< vector<long int>> result_beneficial;
	vector<double> result_fitness;
	for (unsigned int i00=0;i00<popFinal.size();++i00){
		result_genome.push_back({});
		result_genome[i00]=popFinal[i00].neuLocus;
		result_deleterious.push_back({});
		result_deleterious[i00]=popFinal[i00].delLocus;
		result_beneficial.push_back({});
		result_beneficial[i00]=popFinal[i00].benLocus;
		result_fitness.push_back(popFinal[i00].w);
	}

	vector< vector <int>> result_mutator;
	if(mutator_control){
		for (unsigned int i00=0;i00<popFinal.size();++i00){
			result_mutator.push_back({});
			for (int i01=0;i01<nmutator;++i01){
				for (unsigned long int i02=0;i02<popFinal[i00].benLocus.size();++i02){
					if (mutator[i01]==popFinal[i00].benLocus[i02]){
						result_mutator[i00].push_back(mutator[i01]);
					}
				}
				for (unsigned long int i02=0;i02<popFinal[i00].delLocus.size();++i02){
					if (mutator[i01]==popFinal[i00].delLocus[i02]){
						result_mutator[i00].push_back(mutator[i01]);
					}
				}
			}
		}
	}

	{
	filename=fileheader+"_popfinal"+fileformat;

	ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc

	for (unsigned int i00=0;i00<popFinal.size();++i00){
		fout << popFinal[i00].ID << " | ";

		for (unsigned long int i01=0;i01<popFinal[i00].neuLocus.size();++i01){
			fout << popFinal[i00].neuLocus[i01] << " , ";
		}
		fout << " | ";
		for (unsigned long int i01=0;i01<popFinal[i00].delLocus.size();++i01){
			fout << popFinal[i00].delLocus[i01] << " , ";
		}
		fout << " | ";
		for (unsigned long int i01=0;i01<popFinal[i00].benLocus.size();++i01){
			fout << popFinal[i00].benLocus[i01] << " , ";
		}

		fout << " | " << popFinal[i00].w << " | " <<
				popFinal[i00].nOff << " | " <<
				popFinal[i00].rel_mu << " | " << endl;
	}
	fout.close();
	}

	{
	filename=fileheader+"_genome"+fileformat;
	ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc
	for (unsigned int i00=0;i00<popFinal.size();++i00){
		for (unsigned long int i01=0;i01<result_genome[i00].size();++i01){
			fout << result_genome[i00][i01] << " , " ;
		}
		fout << endl;
	}
	fout.close();
	}

	{
	filename=fileheader+"_deleterious"+fileformat;
	ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc
	for (unsigned int i00=0;i00<popFinal.size();++i00){
		for (unsigned long int i01=0;i01<result_deleterious[i00].size();++i01){
			fout << result_deleterious[i00][i01] << " , " ;
		}
		fout << endl;
	}
	fout.close();
	}

	{
	filename=fileheader+"_beneficial"+fileformat;
	ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc
	for (unsigned int i00=0;i00<popFinal.size();++i00){
		for (unsigned long int i01=0;i01<result_beneficial[i00].size();++i01){
			fout << result_beneficial[i00][i01] << " , " ;
		}
		fout << endl;
	}
	fout.close();
	}

	{
	if (mutator_control){
		filename=fileheader+"_mutator"+fileformat;
		ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc
		for (unsigned int i00=0;i00<popFinal.size();++i00){
			for (unsigned long int i01=0;i01<result_mutator[i00].size();++i01){
				fout << result_mutator[i00][i01] << " , " ;
			}
			fout << endl;
		}
		fout.close();
	}
	}

	{
	filename=fileheader+"_meanfitness"+fileformat;
	ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc
	for (unsigned int i00=0;i00<t_bot;++i00){
		for (unsigned  int i01=0;i01<t_growth;++i01){
			fout << result_wMean[i00][i01] << " , " ;
		}
		fout << endl;
	}
	fout.close();
	}

	{
	filename=fileheader+"_offspringnum"+fileformat;
	ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc
	for (unsigned int i00=0;i00<t_bot;++i00){
		for (unsigned  int i01=0;i01<t_growth;++i01){
			fout << result_nOff[i00][i01] << " , " ;
		}
		fout << endl;
	}
	fout.close();
	}

	array<long int, t_bot> result_nNeu;
	array<double, t_bot> result_rNeu;
	array<long int, t_bot> result_nDel;
	array<double, t_bot> result_rDel;
	array<long int, t_bot> result_nBen;
	array<double, t_bot> result_rBen;
	array<long int, t_bot> result_nMut;
	array<double, t_bot> result_w;
	for (int i00=0;i00<tbot;++i00){
		result_nNeu[i00]=result_genome[i00].size();
		result_rNeu[i00]=result_nNeu[i00]/(double)l;
		result_nDel[i00]=result_deleterious[i00].size();
		result_rDel[i00]=result_nDel[i00]/(double)ngenes;
		result_nBen[i00]=result_beneficial[i00].size();
		result_rBen[i00]=result_nBen[i00]/(double)ngenes;
		if (mutator_control){
			result_nMut[i00]=result_mutator[i00].size();
		}
		else{
			result_nMut[i00]=0;
		}
		result_w[i00]=result_fitness[i00];
	}

	{
	filename=fileheader+"_summary"+fileformat;
	ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc
	for (int i00=0;i00<tbot;++i00){
		fout << " | " << result_nNeu[i00] << " | "
				<< result_rNeu[i00] << " | "
				<< result_nDel[i00] << " | "
				<< result_rDel[i00] << " | "
				<< result_nBen[i00] << " | "
				<< result_rBen[i00] << " | "
				<< result_nMut[i00] << " | "
				<< result_muScaled[i00] << " | "
				<< result_w[i00] << " | "
				<< endl;
	}
	fout.close();
	}

	{
	filename=fileheader+"_coefficients"+fileformat;
	ofstream fout(filename.c_str());  // default mode is ios::out | ios::trunc
	for (int i00=0;i00<nmutator;++i00){
		fout << mutator[i00] << " , ";
	}
	fout << endl;
	for (int i00=0;i00<nmutator;++i00){
		fout << mut_e[i00] << " , ";
	}
	fout << endl;
	for (int i00=0;i00<ngenes;++i00){
		fout << s_del[i00] << " , " ;
	}
	fout << endl;
	for (int i00=0;i00<ngenes;++i00){
		fout << s_ben[i00] << " , " ;
	}
	fout << endl;
	fout.close();
	}
	return;
}

int main() {
	clock_t begin=clock();

	string fileheader;
	if (mutator_control && increased_mutator && del_dep){
		fileheader="deldepn";
	}
	if (mutator_control && increased_mutator && !del_dep){
		fileheader="incmttr";
	}
	if (mutator_control && !increased_mutator){
		fileheader="mutator";
	}
	if (!mutator_control){
		fileheader="nomuttr";
	}

	fileheader += "_";
	fileheader += to_string(run);

	long int initialN=1;
	int tbot=1;
	double w_next=1;
	double mu_next=1;

	double lambdab=1/betab;

	double mu_neu=mu_n/l;
	double mu_del=mu_d/l;
	double mu_ben=mu_b/l;

	random_device rd;
	//mt19937 gen (rd());
	minstd_rand0 gen (rd());
	normal_distribution<double> dist1(mu_neu,0.1*mu_neu);
	normal_distribution<double> dist2(mu_del,0.1*mu_del);
	normal_distribution<double> dist3(mu_ben,0.1*mu_ben);
	gamma_distribution<double> dist4(shaped,scaled);
	exponential_distribution<double> dist5(lambdab);
	uniform_int_distribution<long int> dist6(0,ngenes-1);
	normal_distribution<double> dist7(mean_e_mutator,0.01*mean_e_mutator);
	poisson_distribution<int> dist8(inc_mut_rate);
	uniform_int_distribution<int> dist9(0,nmutator-1);
	uniform_int_distribution<long int> dist14(0,l-1);
	uniform_int_distribution<int> dist17(0,1000);

	array<array<bool,1001>,1001> success;
	for (int i00=0;i00<500;++i00){
		for (int i01=0;i01<1001;++i01){
			success[i00][i01]=0;
		}
	}
	for (int i00=500;i00<1001;++i00){
		for (int i01=0;i01<1001;++i01){
			success[i00][i01]=1;
		}
	}

	for (int i00=0;i00<500;++i00){
		for (int i01=0;i01<i00;++i01){
			int place=dist17(gen);
			while (success[i00][place]==1) place=dist17(gen);
			success[i00][place]=1;
		}
	}
	for (int i00=500;i00<1001;++i00){
		for (int i01=0;i01<1001-i00;++i01){
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

	vector<long int> neuLoci;
	vector<long int> delLoci;
	vector<long int> benLoci;

	individual initialind={0,{},{},{},1,1,2};

	long int sumOff;
	long int nhighest;
	long int poolsize;
	long int one=1;
	individual tmpf;
	long int luckyInd;
	vector<wlist> w_list;
	vector<individual> popFinal;

	for (int i1=0;i1<t_bot;++i1){
		clock_t lap0 = clock();
		long int n=initialN;
		{
		vector<individual> pop;
		for (long int i00=0;i00<n;++i00){
			pop.push_back(initialind);
		}
		pop[0].rel_mu=mu_next;
		pop[0].w=w_next;
		pop[0].nOff=round(pop[0].w*2.5-0.1);
		pop[0].nOff=max(0,min(3,pop[0].nOff));

		if (mutator_control && increased_mutator){
			if (del_dep){
				inc_mut_rate=max(0.0001,1-w_next);
			}
			int tmp=dist8(gen);
			for (int i00=0;i00<tmp;++i00){
				int mut;
				mut=dist9(gen);
				if(!(find(benLoci.begin(),benLoci.end(),mut)!=benLoci.end())){
					pop[0].benLocus.push_back(mut);
					for (int i01=0;i01<nmutator;++i01){
						if(mut==mutator[i01]){
							pop[0].rel_mu=pop[0].rel_mu*mut_e[i01];
						}
					}
				}
			}
		}

		sumOff=0;
		for (int i00=0;i00<n;++i00){
			sumOff=pop[i00].nOff;
		}
		if (sumOff<=0){
			generate_results(popFinal,result_nOff, result_wMean,
				result_muScaled, mutator,s_del, s_ben,mut_e,
				fileheader,tbot);
			clock_t end = clock();
			double es = double(end - begin);
			cout << "Extinction at bottleneck" << endl;
			cout << "t_bot= " << i1 << endl;
			cout << "time= " << es << endl;
			exit(0);
		}

		result_muScaled[i1]=mu_next;

		for (int i2=0;i2<t_growth;++i2){
			clock_t lap = clock();
			for (long int i00=0;i00<n;++i00){
				if (pop[i00].nOff>0){
					for (int i01=0;i01<pop[i00].nOff-1;++i01){
						pop.push_back(pop[i00]);
					}
				}
				if (pop[i00].nOff==0){
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
			uniform_int_distribution<long int> dist13(0,n-1);
			long int nNeu=dist10(gen);
			long int nDel=dist11(gen);
			long int nBen=dist12(gen);

			{
				{
				while (nNeu>100000){
					{
					int rmu;
					mutlist neu;
					for (long int i00=0;i00<100000;++i00){
						neu.mutInd=dist13(gen);
						neu.mut=dist6(gen);
						rmu=1000*pop[neu.mutInd].rel_mu;
						rmu=max(0,rmu);
						while (rmu>1000){
							pop[neu.mutInd].neuLocus.push_back(neu.mut);
							neu.mut=dist6(gen);
							rmu-=1000;
						}
						neu.success=success[rmu][dist17(gen)];
						if (neu.success){
							pop[neu.mutInd].neuLocus.push_back(neu.mut);
						}
					}
					}
					nNeu-=100000;
				}
				}
				{
				mutlist neu;
				int rmu;
				for (long int i00=0;i00<nNeu;++i00){
					neu.mutInd=dist13(gen);
					neu.mut=dist14(gen);
					rmu=1000*pop[neu.mutInd].rel_mu;
					rmu=max(0,rmu);
					while (rmu>1000){
						pop[neu.mutInd].neuLocus.push_back(neu.mut);
						neu.mut=dist6(gen);
						rmu-=1000;
					}
					neu.success=success[rmu][dist17(gen)];
					if (neu.success){
						pop[neu.mutInd].neuLocus.push_back(neu.mut);
					}
				}
				}
			}

			{
				{

				while (nDel>100000){
					{
					mutlist del;
					int rmu;
					for (long int i00=0;i00<100000;++i00){
						del.mutInd=dist13(gen);
						del.mut=dist6(gen);
						rmu=1000*pop[del.mutInd].rel_mu;
						rmu=max(0,rmu);
						while (rmu>1000){
							pop[del.mutInd].delLocus.push_back(del.mut);
							del.mut=dist6(gen);
							rmu-=1000;
						}
						del.success=success[rmu][dist17(gen)];
						if (del.success && !(find(delLoci.begin(),delLoci.end(),del.mut)!=delLoci.end())){
							pop[del.mutInd].delLocus.push_back(del.mut);
							pop[del.mutInd].w-=s_del[del.mut];
							pop[del.mutInd].nOff=round(pop[del.mutInd].w*2.5-0.1);
							pop[del.mutInd].nOff=max(0,min(3,pop[del.mutInd].nOff));
							for (int i01=0;i01<nmutator;++i01){
								if(del.mut==mutator[i01]){
									(pop[del.mutInd].rel_mu)/=mut_e[i01];
								}
							}
							}
					}
					}
					nDel-=100000;
				}
				}
				{
				mutlist del;
				int rmu;
				for (long int i00=0;i00<nDel;++i00){
					del.mutInd=dist13(gen);
					del.mut=dist6(gen);
					rmu=1000*pop[del.mutInd].rel_mu;
					rmu=max(0,rmu);
					while (rmu>1000){
						pop[del.mutInd].delLocus.push_back(del.mut);
						del.mut=dist6(gen);
						rmu-=1000;
					}
					del.success=success[rmu][dist17(gen)];
					if (del.success && !(find(delLoci.begin(),delLoci.end(),del.mut)!=delLoci.end())){
						pop[del.mutInd].delLocus.push_back(del.mut);
						pop[del.mutInd].w-=s_del[del.mut];
						pop[del.mutInd].nOff=round(pop[del.mutInd].w*2.5-0.1);
						pop[del.mutInd].nOff=max(0,min(3,pop[del.mutInd].nOff));
						for (int i01=0;i01<nmutator;++i01){
							if(del.mut==mutator[i01]){
								(pop[del.mutInd].rel_mu)/=mut_e[i01];
							}
						}
					}
				}
				}
			}

			{
			mutlist ben;
			int rmu;
			for (long int i00=0;i00<nBen;++i00){
				ben.mutInd=dist13(gen);
				ben.mut=dist6(gen);
				rmu=1000*pop[ben.mutInd].rel_mu;
				while (rmu>1000){
					pop[ben.mutInd].neuLocus.push_back(ben.mut);
					ben.mut=dist6(gen);
					rmu-=1000;
				}
				ben.success=success[rmu][dist17(gen)];
				if (ben.success && !(find(benLoci.begin(),benLoci.end(),ben.mut)!=benLoci.end())){
					pop[ben.mutInd].benLocus.push_back(ben.mut);
					pop[ben.mutInd].w+=s_ben[ben.mut];
					pop[ben.mutInd].nOff=round(pop[ben.mutInd].w*2.5-0.1);
					pop[ben.mutInd].nOff=max(0,min(3,pop[ben.mutInd].nOff));
					for (int i01=0;i01<nmutator;++i01){
						if(ben.mut==mutator[i01]){
							pop[ben.mutInd].rel_mu*=mut_e[i01];
						}
					}
				}
			}

			}
			}

			{
				sumOff=0;
				double wtmp=0;
				for (int i00=0;i00<n;++i00){
					sumOff+=pop[i00].nOff;
					wtmp+=pop[i00].w;
				}

				result_nOff[i1][i2]=sumOff;
				result_wMean[i1][i2]=wtmp/n;
			}

			if (sumOff<=0){
				tbot=i1;
				generate_results(popFinal,result_nOff, result_wMean,
						result_muScaled, mutator,s_del, s_ben,mut_e,
						fileheader,tbot);

				cout << "Extinction at growth" << endl;
				cout << "t_bot= " << i1+1 << endl;
				cout << "t_growth= " << i2+1 << endl;

				clock_t end = clock();
				double es = double(end - begin);
				cout << "time= " << es << endl;
				exit(0);
			}
			if(n>1000000){

				clock_t end = clock();
				double es = double(end - lap);
				cout << "Finished Gen " << i2+1 << " in Bot " << i1+1 << " n= " << n << " time= " << es << endl;
				cout << "Population over a million" << endl;
				break;
			}
			{

			clock_t end = clock();
			double es = double(end - lap);
			cout << "Finished Gen " << i2+1 << " in Bot " << i1+1 << " n= " << n << " time= " << es << endl;
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
		}
		for (unsigned long int i00=0;i00<tmpf.neuLocus.size();++i00){
			neuLoci.push_back(tmpf.neuLocus[i00]);
		}
		tmpf.neuLocus=neuLoci;
		for (unsigned long int i00=0;i00<tmpf.delLocus.size();++i00){
			delLoci.push_back(tmpf.delLocus[i00]);
		}
		tmpf.delLocus=delLoci;
		for (unsigned long int i00=0;i00<tmpf.benLocus.size();++i00){
			benLoci.push_back(tmpf.benLocus[i00]);
		}
		tmpf.benLocus=benLoci;
		w_next=tmpf.w;
		mu_next=tmpf.rel_mu;

		tbot=i1+1;
		popFinal.push_back(tmpf);

		clock_t end = clock();
		double es = double(end - lap0);
		cout << "finished Bot " << i1 << " time= " << es << endl;
	}
	generate_results(popFinal,result_nOff, result_wMean,
			result_muScaled, mutator,s_del, s_ben,mut_e,
			fileheader,tbot);
	clock_t end = clock();
	double es = double(end - begin);
	cout << "End of simulation" << " time= " << es << endl;

	return 0;
}
