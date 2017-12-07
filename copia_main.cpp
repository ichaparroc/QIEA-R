/*ISRAEL CHAPARRO
QIEA-R*/

//#include <iostream.h>

#include <iostream>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <algorithm>

#define DBG 1

using namespace std;

///Mersenne Twister
typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val;           // populate somehow
MyRNG rng;                   // e.g. keep one global instance (per thread)
void initialize()
{
  rng.seed(seed_val);
}
//Mersenne Twister

pair<double,int> get_min(vector<double>&fit_t) //get min and position
{
	double min=fit_t[0];
	int idx=0;

	for(int i=1;i<(int)fit_t.size();i++)
		if(fit_t[i]<min)
		{
			min=fit_t[i];
			idx=i;
		}

	pair<double,int> ans;
	ans.first = min;
	ans.second = idx;

	return ans;
}

void Q_initialization(vector<vector<vector<double>>>&Q_t,double U_LIM,double L_LIM)
{
	//Same Domain
	for(int i=0;i<(int)Q_t.size();i++)
		for(int j=0;j<(int)Q_t[0].size();j++)
		{
			Q_t[i][j][0]=U_LIM-L_LIM; //sigma -- pulse width
			Q_t[i][j][1]=L_LIM+(U_LIM-L_LIM)/2.0; //mu -- square pulse
		}
}

void C_generation(vector<vector<double>>&C_t, vector<vector<vector<double>>>&Q_t)
{
	uniform_real_distribution<double> dist(0,1);
	for(int i=0;i<(int)Q_t.size();i++)
		for(int j=0;j<(int)Q_t[0].size();j++)
			C_t[i][j]=dist(rng)*Q_t[i][j][0]+(Q_t[i][j][1]-Q_t[i][j][0]/2.0);
}


/*
		//inicializa Qt con m individuos de n genes.
		//while t <= T do
			//generar Et observando los individuos Qt
			//if t = 1 then
			//	Ct = Et
			//else
			//recombinar(Et,Ct) -> Et
			evaluar Et
			seleccionar Ct<-k mejores individuos de Et U Ct
			end if
			actualizar Qt+1 con los m mejores individuos de Ct
			t = t + 1
		end while
*/

double evaluation(vector<double>&ind)
{
	double sum=0.0,a,b;
	for(int i=0;i<(int)ind.size();i++)
		sum+=pow(ind[i],2);
	a=pow(sin(sqrt(sum)),2)-0.5;
	b=pow(1+0.001*sum,2);
	return 0.5+(a/b);
}

void evaluation_2_fit(vector<double>&fit_t, vector<vector<double>>&C_t)
{
	for(int i=0;i<(int)C_t.size();i++)
	{
		fit_t[i]=evaluation(C_t[i]);
	}
}

bool sort_evaluation(vector<double>&ind_1, vector<double>&ind_2)
{
	return (evaluation(ind_1)<evaluation(ind_2));
}

void print_Q(vector<vector<vector<double>>>&Q_t)
{
	for(int i=0;i<(int)Q_t.size();i++)
	{
		cout<<endl<<"Q("<<i<<")";
		for(int j=0;j<(int)Q_t[0].size();j++)
			cout<<endl<<"Q("<<i<<")("<<j<<"): o="<<Q_t[i][j][0]<<" u="<<Q_t[i][j][1];
	}
}


void print_C(vector<vector<double>>&C_t)
{
	for(int i=0;i<(int)C_t.size();i++)
	{
		cout<<endl<<"C("<<i<<"):";
		for(int j=0;j<(int)C_t[0].size();j++)
			cout<<" "<<C_t[i][j];
	}
}


void print_fit(vector<double>&fit_t)
{
	cout<<endl<<"Fit:";
	for(int i=0;i<(int)fit_t.size();i++)
		cout<<" "<<fit_t[i];
}

int main(int argc,char **argv)
{

	int T_POB=20,N_DIM=2,U_LIM=10,L_LIM=-10,N_GEN=100;
	double P_CRO=0.5,P_UD=1;

	double best_fit_t,best_fit_t_1,best_global;
	int best_pos=0; //elitism does that the best are in the first position ever

	vector<vector<vector<double>>>Q_t(T_POB,vector<vector<double>>(N_DIM,vector<double>(2)));

	vector<vector<double>>C_t(T_POB,vector<double>(N_DIM));
	vector<vector<double>>E_t(T_POB,vector<double>(N_DIM));

	vector<vector<double>>E_U_C(T_POB*2,vector<double>(N_DIM));

	vector<double>fit_t(T_POB);
	vector<double>fit_t_1(T_POB);

	vector<double>best_ind(N_DIM);

	uniform_int_distribution<int>dist_T_POB(0,T_POB-1); 	//To get index of classical individuals
	uniform_real_distribution<double>dist_01(0,1);	//To get prob of crossover

	//t=1
	int t=0;

	//Initialize Q(t) with m individuals from n genes
	Q_initialization(Q_t,U_LIM,L_LIM);

	if(DBG) {cout<<endl<<"Q inicial";} print_Q(Q_t);


	//while t<=T
	while(t<N_GEN)
	{
		if(DBG) cout<<endl<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::t="<<t;

		//generate E(t) observing the Q(t) individuals
		C_generation(E_t,Q_t);

		if(DBG) {cout<<endl<<"E::::::::::"; print_C(E_t);}

		//if t=1 then
		if(t==0)
		{
			//C(t)=E(t)
			C_t=E_t;

			evaluation_2_fit(fit_t,C_t);

			//eliminar a futuro;
			print_fit(fit_t);

			sort(C_t.begin(),C_t.end(),sort_evaluation);


			evaluation_2_fit(fit_t,C_t);


			//eliminar a futuro
			print_C(C_t);

			best_fit_t=fit_t[0];

			if(DBG) {print_fit(fit_t); cout<<"best: "<<best_fit_t;}
		}
		//else
		else
		{
			//recombine {E(t),C(t)} -> E(t) -- Using Differential Evolution: crossover

			int father,mother;	//parents
			vector<double> x_c1(N_DIM);//To save crossover1
            vector<double> x_c2(N_DIM);//To save crossover2

			for(int i=0;i<T_POB;i++)
			{
				double random=dist_01(rng);
				if(random<=P_CRO)
				{
					do
						father=dist_T_POB(rng);
					while(father==i); //bestpos ?

					do
						mother=dist_T_POB(rng);
					while(mother==i||mother==father); //bestpos ?

					for(int j=0;j<N_DIM;j++)
					{

	                    double alpha=dist_01(rng);
	                    x_c1[j]=alpha*E_t[father][j]+(1-alpha)*C_t[mother][j];
	                    x_c2[j]=alpha*C_t[mother][j]+(1-alpha)*E_t[father][j];
					}

					E_t[father]=x_c1; //for all dimensions :)
	                E_t[mother]=x_c2; //for all dimensions :)
				}
			}

			if(DBG) {cout<<endl<<"recombinar(E,C) -> E"; print_C(E_t);}

			//evaluate E(t)

			fit_t_1=fit_t;

			evaluation_2_fit(fit_t,E_t);

			best_fit_t_1=best_fit_t;

			best_fit_t=get_min(fit_t).first;

			if(DBG) {print_fit(fit_t); cout<<endl<<"best:"<<best_fit_t;}

			//select C(t) <- k best individuals of E(t) U C(t)

			int idx=0;
			for(int i=0;i<(int)E_t.size();i++)
			{
				E_U_C[idx]=E_t[i];
				idx++;
			}
			for(int i=0;i<(int)C_t.size();i++)
			{
				E_U_C[idx]=C_t[i];
				idx++;
			}

			sort(E_U_C.begin(),E_U_C.end(),sort_evaluation);

			for(int i=0;i<(int)C_t.size();i++)
				C_t[i]=E_U_C[i];

			best_global=evaluation(C_t[0]);

			if(DBG) cout<<endl<<"BEST GLOBAL: "<<best_global;

			if(DBG) {cout<<endl<<"C <- E U C : "; print_C(C_t);}
		}

		//update Q(t+1) with the individuals of C(t)
        double random=dist_01(rng);

        if(random<P_UD && t!=0)
        {
            //Pulse width adjustment (sigma) -- Rule 1/5th
            double delta=dist_01(rng);
            int counter = 0;

            for(int i=0;i<T_POB;i++) //Josimar Chire's Strategy (he have a extitive gap generational)
				if(fit_t[i]<best_fit_t_1)
					counter++;
            double phi=(double)counter/T_POB;
			cout<<endl<<"phi"<<phi;
			double mul_fact = 1.0;

            if(phi<0.2)
                mul_fact=delta;
            else if(phi>0.2)
                mul_fact=1.0/delta;

			cout<<endl<<"mul_fact="<<mul_fact;

            for(int i=0;i<T_POB;i++)
                for( int j=0;j<N_DIM;j++)
                    Q_t[i][j][0]*=mul_fact;

        	//Moving Square pulse (mu)
            double lambda=dist_01(rng);

            for (int i=0;i<T_POB;i++)
                for (int j=0;j<N_DIM;j++)
                    Q_t[i][j][1]+=lambda*(C_t[best_pos][j]-Q_t[i][j][1]);
        }

		//t=t+1
		t++;
	}
}
