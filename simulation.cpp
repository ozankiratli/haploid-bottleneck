//============================================================================
// Name        : simulation.cpp
// Author      : ozan kiratli - ozankiratli@gmail.com
// Version     : 0.1
// Copyright   : Your copyright notice
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
#include <chrono>
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
	long int nAMut;
	double muScale;
	double w;
};

inline double fastPow(double a, double b) {
  int e = (int) b;
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;

  double r = 1.0;
  while (e) {
    if (e & 1) {
      r *= a;
    }
    a *= a;
    e >>= 1;
  }
  return r * u.d;
}



string fileformat=".out";

bool acompare(wlist lhs, wlist rhs) { return lhs.w > rhs.w; }

const int t_bot = T_bot;
const int t_growth = T_growth;
const long int n_growth = N_growth;
const long int ngenes = Ngenes;
const long int l = L;
double mu_d = Mu_d;
double mu_b = Mu_b;
double mu_n = Mu_n;
double shaped = Shaped;
double scaled = Scaled;
double betab = Betab;
double w_sel = W_sel;
const int nmutator = Nmutator;
double mean_e_mutator = Mean_e_mutator;
double var_e_mutator = Var_e_mutator;
double selecthighest = Selecthighest;
bool mutator_control = Mutator_control;

int run = Run;
array <double,4> w_to_off = {WtoOffa,WtoOffb,WtoOffc,WtoOffd};

inline int wtooffn(double fitness, double rand){
    double n1;
    double li1=1.0;
    n1=rand+w_to_off[0];
    for (int i=1;i<4;i++){
        li1*=fitness;
        n1+=w_to_off[i]*li1;
    }
    n1=max(0.,min(3.,n1));
    int nfinal=round(n1);
    return nfinal;
}

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
	vector< vector <int>> result_antimutator;
	if(mutator_control){
		for (unsigned int i00=0;i00<popFinal.size();++i00){
			result_antimutator.push_back({});
			for (int i01=0;i01<nmutator;++i01){
				for (unsigned long int i02=0;i02<popFinal[i00].benLocus.size();++i02){
					if (mutator[i01]==popFinal[i00].benLocus[i02]){
						result_antimutator[i00].push_back(mutator[i01]);
					}
				}
			}
		}
		for (unsigned int i00=0;i00<popFinal.size();++i00){
			result_mutator.push_back({});
			for (int i01=0;i01<nmutator;++i01){
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
	array<long int, t_bot> result_nAMt;
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
			result_nAMt[i00]=result_antimutator[i00].size();
		}
		else{
			result_nMut[i00]=0;
			result_nAMt[i00]=0;
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
					<< result_nAMt[i00] << " | "
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
	auto begin = chrono::high_resolution_clock::now();
	string fileheader;
	if (mutator_control){
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
    	normal_distribution<double> rand_Neumu(mu_neu,0.1*mu_neu);
    	normal_distribution<double> rand_Delmu(mu_del,0.1*mu_del);
    	normal_distribution<double> rand_Benmu(mu_ben,0.1*mu_ben);
    	gamma_distribution<double> rand_deleffects(shaped,scaled);
    	exponential_distribution<double> rand_beneffects(lambdab);
    	uniform_int_distribution<long int> rand_ngenes(0,ngenes-1);
    	normal_distribution<double> rand_mutatoreffect(mean_e_mutator,var_e_mutator);
    	uniform_int_distribution<int> rand_nMutator(0,nmutator-1);
    	uniform_int_distribution<long int> rand_genome(0,l-1);
    	uniform_int_distribution<int> rand_success_bool(0,1000);
	normal_distribution<double> rand_offvar(0,w_sel);

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
            		int place=rand_success_bool(gen);
            		while (success[i00][place]==1) place=rand_success_bool(gen);
			success[i00][place]=1;
		}
	}
	for (int i00=500;i00<1001;++i00){
		for (int i01=0;i01<1001-i00;++i01){
            		int place=rand_success_bool(gen);
            		while (success[i00][place]==0) place=rand_success_bool(gen);
			success[i00][place]=0;
		}
	}
	array<array<double,t_growth>,t_bot> muNeu;
	array<array<double,t_growth>,t_bot> muDel;
	array<array<double,t_growth>,t_bot> muBen;
	for (int i00=0; i00<t_bot ; ++i00){
		for (int i01=0; i01<t_growth; ++i01 ){
            		muNeu[i00][i01]=rand_Neumu(gen);
            		muDel[i00][i01]=rand_Delmu(gen);
            		muBen[i00][i01]=rand_Benmu(gen);
		}
	}

	array<double,ngenes> s_del;
	array<double,ngenes> s_ben;
	for (int i00=0;i00<ngenes;++i00){
        	s_del[i00]=rand_deleffects(gen);
        	s_ben[i00]=rand_beneffects(gen);
	}

	array<int,nmutator> mutator;
	array<double,nmutator> mut_e;
	if (mutator_control){
		for (int i00=0;i00<nmutator;++i00){
        		mutator[i00]=rand_ngenes(gen);
            		mut_e[i00]=rand_mutatoreffect(gen);
			s_del[mutator[i00]]=0;
			s_ben[mutator[i00]]=0;
		}
	}
	else{
		for (int i00=0;i00<nmutator;++i00){
            		mutator[i00]=rand_ngenes(gen);
			mut_e[i00]=1;
		}
	}


	array<array<int,t_growth>,t_bot> result_nOff;
	array<array<double,t_growth>,t_bot> result_wMean;
	array<double,t_bot> result_muScaled;

	vector<long int> neuLoci;
	vector<long int> delLoci;
	vector<long int> benLoci;
    	auto delLocibegin=delLoci.begin();
    	auto delLociend=delLoci.end();
    	auto benLocibegin=benLoci.begin();
    	auto benLociend=benLoci.end();

	individual initialind={0,{},{},{},1,1,2};

	long int sumOff;
	long int nhighest;
	long int poolsize;
	long int one=1;
	vector<individual> popFinal;

	for (int i1=0;i1<t_bot;++i1){
        auto lap0 = chrono::high_resolution_clock::now();
        auto mutatorbegin=mutator.begin();
        auto mutatorend=mutator.end();
		long int n=initialN;
		{
			vector<individual> pop;
			for (long int i00=0;i00<n;++i00){
				pop.push_back(initialind);
			}
			pop[0].rel_mu=mu_next;
			pop[0].w=w_next;
            		double rndnum=rand_offvar(gen);
        		pop[0].nOff=wtooffn(pop[0].w,rndnum);
           		sumOff=0;
			for (int i00=0;i00<n;++i00){
				sumOff=pop[i00].nOff;
			}
            		if (sumOff<=0){
                		generate_results(popFinal,result_nOff, result_wMean,
                			result_muScaled, mutator,s_del, s_ben,mut_e,
                        		fileheader,tbot);
                		auto end = chrono::high_resolution_clock::now();
                		auto es = end - begin;
                		cout << "Extinction at bottleneck" << endl;
                		cout << "t_bot= " << i1 << endl;
                		cout << "time= " << chrono::duration_cast<chrono::seconds>(es).count() << " seconds" << endl;
                		exit(0);
            		}

			for (int i2=0;i2<t_growth;++i2){
                		auto lap = chrono::high_resolution_clock::now();
				for (long int i00=0;i00<n;++i00){
					if (pop[i00].nOff>0){
						for (int i01=0;i01<pop[i00].nOff-1;++i01){
							pop.push_back(pop[i00]);
						}
					}
					if (pop[i00].nOff<=0){
						pop.erase(pop.begin()+i00);
					}
				}
				n=pop.size();
				for (long int i00=0;i00<n;++i00){
					pop[i00].ID=i00;
				}
				{
                    			poisson_distribution<long int> rand_nNeu(n*l*muNeu[i1][i2]);
                    			poisson_distribution<long int> rand_nDel(n*l*muDel[i1][i2]);
                    			poisson_distribution<long int> rand_nBen(n*l*muBen[i1][i2]);
                    			uniform_int_distribution<long int> rand_npop(0,n-1);
                    			long int nNeu=rand_nNeu(gen);
                    			long int nDel=rand_nDel(gen);
                    			long int nBen=rand_nBen(gen);
                    			// Neutral Mutations
					{
						{
							while (nNeu>100000){
								{
									int rmu;
									mutlist neu;
									for (long int i00=0;i00<100000;++i00){
                                        					neu.mutInd=rand_npop(gen);
                                        					neu.mut=rand_ngenes(gen);
										rmu=1000*pop[neu.mutInd].rel_mu;
										rmu=max(0,rmu);
										while (rmu>1000){
											pop[neu.mutInd].neuLocus.push_back(neu.mut);
                                            						neu.mut=rand_ngenes(gen);
											rmu-=1000;
										}
                                        					neu.success=success[rmu][rand_success_bool(gen)];
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
                                				neu.mutInd=rand_npop(gen);
                                				neu.mut=rand_genome(gen);
								rmu=1000*pop[neu.mutInd].rel_mu;
								rmu=max(0,rmu);
								while (rmu>1000){
									pop[neu.mutInd].neuLocus.push_back(neu.mut);
        			                   			neu.mut=rand_ngenes(gen);
									rmu-=1000;
								}
                                				neu.success=success[rmu][rand_success_bool(gen)];
								if (neu.success){
									pop[neu.mutInd].neuLocus.push_back(neu.mut);
								}
							}
						}
					}
					//Deleterious Mutations
					{
						{

                            				while (nDel>100000){
                                				{
                                    					mutlist mutt;
                                    					int rmu;
                                    					for (long int i00=0;i00<100000;++i00){
                                        					mutt.mutInd=rand_npop(gen);
                                        					mutt.mut=rand_ngenes(gen);
                                        					rmu=1000*pop[mutt.mutInd].rel_mu;
                                        					rmu=max(0,rmu);
                                        					while (rmu>1000){
                                            						if (find(delLocibegin,delLociend,mutt.mut)==delLociend){
                                                						pop[mutt.mutInd].delLocus.push_back(mutt.mut);
                                                						pop[mutt.mutInd].w*=(1-s_del[mutt.mut]);
                                                						double rndnum=rand_offvar(gen);
                                                						pop[mutt.mutInd].nOff=wtooffn(pop[mutt.mutInd].w,rndnum);
                                                						auto found=find(mutatorbegin,mutatorend,mutt.mut);
                                                						if(found!=mutatorend){
                                                    							int pos=distance(mutatorbegin,found);
                                                    							pop[mutt.mutInd].rel_mu/=mut_e[pos];
                                                    							pop[mutt.mutInd].rel_mu=max(0.0001,pop[mutt.mutInd].rel_mu);
                                                						}
                                            						}
                                            						mutt.mut=rand_ngenes(gen);
                                            						rmu-=1000;
                                        					}
                                        					mutt.success=success[rmu][rand_success_bool(gen)];
                                        					if (mutt.success && (find(delLocibegin,delLociend,mutt.mut)==delLociend)){
                                            						pop[mutt.mutInd].delLocus.push_back(mutt.mut);
                                            						pop[mutt.mutInd].w*=(1-s_del[mutt.mut]);
                                        						double rndnum=rand_offvar(gen);
                                            						pop[mutt.mutInd].nOff=wtooffn(pop[mutt.mutInd].w,rndnum);
                                            						auto found=find(mutatorbegin,mutatorend,mutt.mut);
                                            						if(found!=mutatorend){
                                                						int pos=distance(mutatorbegin,found);
                                                						pop[mutt.mutInd].rel_mu/=mut_e[pos];
                                        			   				pop[mutt.mutInd].rel_mu=max(0.0001,pop[mutt.mutInd].rel_mu);
                                            						}
                                        					}
                                    					}
                                				}
                                				nDel-=100000;
                            				}
						}
						{
                            				mutlist mutt;
                            				int rmu;
                            				for (long int i00=0;i00<nDel;++i00){
                                				mutt.mutInd=rand_npop(gen);
                                				mutt.mut=rand_ngenes(gen);
                                				rmu=1000*pop[mutt.mutInd].rel_mu;
                                				rmu=max(0,rmu);
                                				while (rmu>1000){
                                    					if (find(delLocibegin,delLociend,mutt.mut)==delLociend){
                                        				pop[mutt.mutInd].delLocus.push_back(mutt.mut);
                                        				pop[mutt.mutInd].w*=(1-s_del[mutt.mut]);
                                        				double rndnum=rand_offvar(gen);
                                        				pop[mutt.mutInd].nOff=wtooffn(pop[mutt.mutInd].w,rndnum);
                                        				auto found=find(mutatorbegin,mutatorend,mutt.mut);
                                        				if(found!=mutatorend){
                                            					int pos=distance(mutatorbegin,found);
                                            					pop[mutt.mutInd].rel_mu/=mut_e[pos];
                                            					pop[mutt.mutInd].rel_mu=max(0.0001,pop[mutt.mutInd].rel_mu);
                                        				}
                                    				}
                                    					mutt.mut=rand_ngenes(gen);
                                    					rmu-=1000;
                                				}
                                				mutt.success=success[rmu][rand_success_bool(gen)];
                                				if (mutt.success && (find(delLocibegin,delLociend,mutt.mut)==delLociend)){
                                    					pop[mutt.mutInd].delLocus.push_back(mutt.mut);
                                    					pop[mutt.mutInd].w*=(1-s_del[mutt.mut]);
                                    					double rndnum=rand_offvar(gen);
                                    					pop[mutt.mutInd].nOff=wtooffn(pop[mutt.mutInd].w,rndnum);
                                    					auto found=find(mutatorbegin,mutatorend,mutt.mut);
                                    					if(found!=mutatorend){
                                        					int pos=distance(mutatorbegin,found);
                                        					pop[mutt.mutInd].rel_mu/=mut_e[pos];
                                        					pop[mutt.mutInd].rel_mu=max(0.0001,pop[mutt.mutInd].rel_mu);
                                    					}
                                				}
                            				}
                        			}	
                    			}
					//Beneficial Mutations
					{
						{
                        				while (nBen>100000){
                            					{
                                					mutlist mutt;
                                					int rmu;
                                					for (long int i00=0;i00<100000;++i00){
                                    						mutt.mutInd=rand_npop(gen);
                                    						mutt.mut=rand_ngenes(gen);
                                    						rmu=1000*pop[mutt.mutInd].rel_mu;
                                    						while (rmu>1000){
                                        						if (find(benLocibegin,benLociend,mutt.mut)==benLociend){
                                        							pop[mutt.mutInd].benLocus.push_back(mutt.mut);
                                            							pop[mutt.mutInd].w*=(1+s_ben[mutt.mut]);
                                            							double rndnum=rand_offvar(gen);
                                            							pop[mutt.mutInd].nOff=wtooffn(pop[mutt.mutInd].w,rndnum);
                                            							auto found=find(mutatorbegin,mutatorend,mutt.mut);
                                            							if(found!=mutatorend){
                                                							int pos=distance(mutatorbegin,found);
                                                							pop[mutt.mutInd].rel_mu*=mut_e[pos];
                                                							pop[mutt.mutInd].rel_mu=max(0.0001,pop[mutt.mutInd].rel_mu);
                                            							}
                                        						}
                                        						mutt.mut=rand_ngenes(gen);
                                        						rmu-=1000;
                                    						}
                                    						mutt.success=success[rmu][rand_success_bool(gen)];
                                    						if (mutt.success && (find(benLocibegin,benLociend,mutt.mut)==benLociend)){
                                        						pop[mutt.mutInd].benLocus.push_back(mutt.mut);
                                        						pop[mutt.mutInd].w*=(1+s_ben[mutt.mut]);
                                        						double rndnum=rand_offvar(gen);
                                        						pop[mutt.mutInd].nOff=wtooffn(pop[mutt.mutInd].w,rndnum);
                                        						auto found=find(mutatorbegin,mutatorend,mutt.mut);
                                        						if(found!=mutatorend){
                                            							int pos=distance(mutatorbegin,found);
                                            							pop[mutt.mutInd].rel_mu*=mut_e[pos];
                                            							pop[mutt.mutInd].rel_mu=max(0.0001,pop[mutt.mutInd].rel_mu);
                                        						}	
                                    						}
                                					}
                        					}
                        					nBen-=100000;
                					}
        					}
                    				{
                        				mutlist mutt;
                        				int rmu;
                        				for (long int i00=0;i00<nBen;++i00){
                            					mutt.mutInd=rand_npop(gen);
                        	 				mutt.mut=rand_ngenes(gen);
                            					rmu=1000*pop[mutt.mutInd].rel_mu;
                            					while (rmu>1000){
                                					if (find(benLocibegin,benLociend,mutt.mut)==benLociend){
                                    						pop[mutt.mutInd].benLocus.push_back(mutt.mut);
                                    						pop[mutt.mutInd].w*=(1+s_ben[mutt.mut]);
                                    						double rndnum=rand_offvar(gen);
                                    						pop[mutt.mutInd].nOff=wtooffn(pop[mutt.mutInd].w,rndnum);
                                    						auto found=find(mutatorbegin,mutatorend,mutt.mut);
                                    						if(found!=mutatorend){
                                        						int pos=distance(mutatorbegin,found);
                                        						pop[mutt.mutInd].rel_mu*=mut_e[pos];
                                        						pop[mutt.mutInd].rel_mu=max(0.0001,pop[mutt.mutInd].rel_mu);
                                    						}
                                					}
                                					mutt.mut=rand_ngenes(gen);
                                					rmu-=1000;
                            					}
                            					mutt.success=success[rmu][rand_success_bool(gen)];
                            					if (mutt.success && (find(benLocibegin,benLociend,mutt.mut)==benLociend)){
                                					pop[mutt.mutInd].benLocus.push_back(mutt.mut);
                                					pop[mutt.mutInd].w*=(1+s_ben[mutt.mut]);
                                					double rndnum=rand_offvar(gen);
                                					pop[mutt.mutInd].nOff=wtooffn(pop[mutt.mutInd].w,rndnum);
                                					auto found=find(mutatorbegin,mutatorend,mutt.mut);
                                					if(found!=mutatorend){
                                    						int pos=distance(mutatorbegin,found);
                                    						pop[mutt.mutInd].rel_mu*=mut_e[pos];
                                    						pop[mutt.mutInd].rel_mu=max(0.0001,pop[mutt.mutInd].rel_mu);
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
                    			auto end = chrono::high_resolution_clock::now();
                    			auto es = end - begin;
                    			cout << "time= " << chrono::duration_cast<chrono::seconds>(es).count() << " seconds" << endl;
                			exit(0);
                		}
                		if(n>n_growth){
                    			auto end = chrono::high_resolution_clock::now();
                    			auto es = end - lap;
                    			cout << "Finished Gen " << i2+1 << " in Bot " << i1+1 << " n= " << n << " time= " << chrono::duration_cast<chrono::milliseconds>(es).count() << " milliseconds" << endl;
                    			cout << "Population over " << n_growth << endl;
                    			break;
                		}
                		{
                    			auto end = chrono::high_resolution_clock::now();
                    			auto es = end - lap;
                    			cout << "Finished Gen " << i2+1 << " in Bot " << i1+1 << " n= " << n << " time= " << chrono::duration_cast<chrono::milliseconds>(es).count() << " milliseconds" << endl;
                		}
			}

			{
				individual tmpf;
				long int luckyInd;
				{
					vector<wlist> w_list;
					wlist ind;
					for (long int i00=0;i00<n;++i00){
						ind.ID=pop[i00].ID;
						ind.w=pop[i00].w;
						w_list.push_back(ind);
					}
					sort(w_list.begin(),w_list.end(),acompare);
					nhighest = round(n*selecthighest);
					poolsize = max(one,nhighest);
                    			uniform_int_distribution<long int> rand_indpool(0,poolsize-1);
                			long int tmp = rand_indpool(gen);
					luckyInd = w_list[tmp].ID;
				}
				tmpf=pop[luckyInd];

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
                		delLocibegin=delLoci.begin();
                		delLociend=delLoci.end();
                		benLocibegin=benLoci.begin();
                		benLociend=benLoci.end();
				{
					double tmprm=1.0;
                    			for (int i00=0;i00<nmutator;++i00){
                        			auto founddel = find(delLocibegin,delLociend,mutator[i00]);
                        			if(founddel!=delLociend){
                            				tmprm/=mut_e[i00];
                            				tmprm=max(0.0001,tmprm);
                        			}
                        			auto foundben = find(benLocibegin,benLociend,mutator[i00]);
                        			if(foundben!=benLociend){
                            				tmprm*=mut_e[i00];
                            				tmprm=max(0.0001,tmprm);
                        			}
					}
					tmpf.rel_mu=tmprm;
				}
				w_next=tmpf.w;
				mu_next=tmpf.rel_mu;
				tbot=i1+1;
				popFinal.push_back(tmpf);
				result_muScaled[i1]=mu_next;
			}
		}

        	auto end = chrono::high_resolution_clock::now();
        	auto es = end - lap0;
        	cout << "finished Bot " << i1+1 << " time= " << chrono::duration_cast<chrono::seconds>(es).count() << " seconds" << endl;
    	}
    	generate_results(popFinal,result_nOff, result_wMean,
            result_muScaled, mutator,s_del, s_ben,mut_e,
            fileheader,tbot);
    	auto end = chrono::high_resolution_clock::now();
    	auto es = end - begin;
	cout << "End of simulation" << " time= " << chrono::duration_cast<chrono::seconds>(es).count() << " seconds" <<endl;

	return 0;
}
