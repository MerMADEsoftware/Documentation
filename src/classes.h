#pragma once
//#define _HAS_TR1_NAMESPACE 1
//#if _HAS_TR1_NAMESPACE
//#endif

//#define _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING 1
//#if _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
//#endif

#ifndef CLASSES_H
#define CLASSES_H


#include <iostream>
#include <fstream>  //stream class to both read and write from/to files
#include <sstream>
#include <string>
#include <math.h> //for atan2 and pow
#include <random> //for distributions
//#include <dinkumware/random>
#include <iterator>
#include <vector>
#include <time.h>
#include <ctime>
#include <queue>
#include <algorithm>

using namespace std;

//declaration
class Patch;
class Cell;
class Population;
class Subpopulation;
class Individual;
class Landscape;
class Model;

//random number generator
//if debugging
//tr1::mt19937 gen(666);
extern mt19937 eng;
//tr1::normal_distribution<double> norm_dist;
//tr1::poisson_distribution<double> pois_dist;
extern uniform_real_distribution<double> unif_dist;

struct mat333 { double neigh[3][3][3]; }; //declaring it here so that i can use it in various classes
//void fill_mat333(mat333* mat, float base, float phi);

extern double M_PI;
extern double M2_PI;
//definition


class Population{
private:

public:
	Population();
	~Population();

	struct pop_info {
		float R; //intrinsic growth rate
		float bc; //competition coefficient, type of density dependence 
		short max_age; //lifespan in years
		short Nrepro; //number of reproductive events per year
		short RepInt; //interval between reproductive events; how many repro events does it skip before reproducing again
		float PRep; //probability of reproducing 
		float juv_size, juv_size_sd; //juv_size is in cm, if iiv_juv_size==1, each indiv will sample its size from a normal dist with this sd
		short juv_rel;
		short iiv_juv_size; //inter-individual variation?
		short survsched; //survival schedule, 0=between reproductive events, 1=annually (if there is only one season per year, then 0 and 1 are the same)
		short Hsize; //harem size, 1=monogamous, >1=harem . this is required if Repro=2
		short fecdensdep, devdensdep, survdensdep; //is density dependence true in fecundity, development or survival
		float devdenscoef=-9, survdenscoef=-9; //only required if the related density dependence is true
		float fishing_mort = -9;

		short DepthPref; //0= set min and max, 1=temperature dependent depth preference
		short DepthSex; //0=depth preference not sex dependent, 1= sex dependent
		short DepthStage; //0=depth preference is not stage dependent, 1=yes
		short Depthmin, Depthmax; //this is for the whole population so encompasses male and female preferences, should it be sex-dependent

		int total_births = 0; //this is to keep track of ID numbers
		struct evo_info {
			//emigration
			int ep_evo = 0; float ep_mut_size = -9, ep_mut = -9, alpha_mut_size = -9, beta_mut_size = -9; //does emigration rate evolve? if yes (ep_evo==1), what is mutation size and rate?
			//transfer
			int active_evo = 0; float active_mut = -9, active_mut_size = -9; //does this parameter evolve? if yes, at what probability and by how much does it mutate?
			int dv_evo = 0; float dv_mut = -9, dv_mut_size = -9; //minimum size at diel vertical migration
																			 //growth
			int growth_evo = 0; float growth_mut = -9, m_mut_size = -9, b_mut_size = -9; //does the growth function evolve? if so, what is the mutation rate, and size for slope and intercept
			float Linf_mut_size = -9, GK_mut_size = -9, Ti_mut_size = -9;
																						 //settlement
			int comp_evo = 0; float comp_mut = -9, comp_mut_size = -9; //does competency time/size evolve? if so, at what probability and by how much?
			int S0_evo = 0; float S0_mut = -9, S0_mut_size = -9, alphaS_mut_size = -9, betaS_mut_size = -9; //does settlement prob evolve? if so, at what prob and by how much? if settlement is density dependent, how much do alpha and beta change if a mutation occurs?
		} pop_evo;

	} pop_dyn;
	
	struct emig_info {
		float f_emig_prob = -9, f_D0 = -9, f_alpha = -9, f_beta = -9; //set emig prob if emigration is density independent; emigration distribution parameters if density dependent (BUT NOT STAGE DEPENDENT)
		float m_emig_prob = -9, m_D0 = -9, m_alpha = -9, m_beta = -9;  //in case emigration is sex dependent but not stage dependent
		int densdep=0,  stagedep=0, sexdep=0, indvar=0, emigstage=0; //default emigration stage is 0
		int ep_evo=0; float ep_mut_size=-9, ep_mut=-9, alpha_mut_size=-9, beta_mut_size=-9; //does emigration rate evolve? if yes (ep_evo==1), what is mutation size and rate?
		float iiv_sd=-9;
	}; // pop_emig;

	struct growth_info {
		int method = 0; //can be linear (0), gompertz (1) or von bertalanffy (2)
		int sex_dep = 0, stage_dep = 0;; //sex dependent growth? 0=no, 1=yes
		float m_m = -9, m_m_sd = -9, f_m = -9, f_m_sd = -9, f_b = -9, m_b = -9; //slope of growth function, intercept of growth function
		int growth_evo = 0; float growth_mut = -9, m_mut_size = -9, b_mut_size = -9; //does the growth function evolve? if so, what is the mutation rate, and size for slope and intercept
		float m_Linf = -9, m_Linf_sd = -9, m_G_K = -9, f_Linf = -9, f_Linf_sd = -9, f_G_K = -9, Linf_mut_size = -9, GK_mut_size = -9, Ti_sd = -9;
		int Ti = -9; //this is the inflection point of the gompertz function, equates to earliest settlement date 

	}; // pop_grow;

	struct trans_info {
		float f_rho=-9, f_rho_sd=-9, f_SL=-9, f_comp_SL=-9, f_SL_sd=-9; //if organism is passive, SL=1, rho=.9. 
		float m_rho=-9, m_rho_sd=-9, m_SL=-9, m_comp_SL=-9, m_SL_sd=-9;

		int stagedep=0, sexdep=0, distmort=0; //true false
		float f_mortprob=-9, f_slope=-9, f_inflpoint=-9; //if distmort==1, slope/inflpoint are the parameters for the reaction norm. else constant mortality probability
		float m_mortprob=-9, m_slope=-9, m_inflpoint=-9;
		int dp=0; int memory=1; //directional persistence and memory
		float m_min_active_time=-9, f_min_active_time=-9, m_min_active_size=-9, f_min_active_size=-9, f_active_sd=-9, m_active_sd=-9;//minimum time/size at which indivs switch to active movement,
		//float m_m, m_m_sd, f_m, f_m_sd, f_b, m_b; //slope of growth function, intercept of growth function
		//int growth_evo; float growth_mut, m_mut_size, b_mut_size; //does the growth function evolve? if so, what is the mutation rate, and size for slope and intercept
		int active_evo=0; float active_mut=-9, active_mut_size=-9; //does this parameter evolve? if yes, at what probability and by how much does it mutate?
		struct growth_info pop_grow;
		int diel_vert = 0; //is behaviour size dependent or time dependent? is the growth function temperature dependent? do they undergo diel vertical migration?
			//if behaviour is size dependent, a growth function is applied during dispersal, otherwise time limits in days are applied
		int dv_active = 0; //does dvm apply to the active stage of dispersal?
		float min_dv_time=-9, min_dv_size = -9, min_dv_sd=-9; //if diel_vert==1, at what point do indivs start dvm? allows flexibility for passive->passive with dvm->active dispersal. is there iiv?
		int down_bias = 0; //is there a downward bias in dispersal movements, ie tracking the seafloor
		int buoy_min=-9, buoy_max=-9, dv_range=-9; //what are the upper and deep limits of buoyancy, and what is the diel vertical migration range?
		//this is calculated as buoy_min+dv_range and buoy_max-dv_range
	}; // pop_trans;

	struct sett_info {
		int stagedep=0, sexdep=0, densdep=0; //true false, inter individual variability
		int f_findmate=0, m_findmate=0; //
		float f_pld = -9, f_pld_sd = -9, f_comp = -9, f_comp_size = -9, f_comp_sd = -9, m_pld = -9, m_pld_sd = -9, m_comp = -9, m_comp_size = -9, m_comp_sd = -9; //pelagic larval duration in days; time to competency in days, must be less than pld, only after this point is settlement possible
		//comp is when behaviour is time-dependent, comp_size is when behaviour is size-dependent
		int comp_evo=0; float comp_mut=-9, comp_mut_size=-9; //does competency time/size evolve? if so, at what probability and by how much?
		//the sd is for when there is inter indiv variability
		float f_S0=-9, f_alphaS=-9, f_betaS=-9, f_sett_sd=-9; //if densdep=1, these are the parameters for the probability distribution
		float m_S0=-9, m_alphaS=-9, m_betaS=-9, m_sett_sd=-9; //the sd applies if there is interindiv variation and will be applied to all parameters. if no iiv, then =-9
		int S0_evo=0; float S0_mut=-9, S0_mut_size=-9, alphaS_mut_size=-9, betaS_mut_size=-9; //does settlement prob evolve? if so, at what prob and by how much? if settlement is density dependent, how much do alpha and beta change if a mutation occurs?
		int buffer = -9; //0=no, yes=1
		int buffer_xycells = 0; //detection distance in number of cells
		int buffer_zlayers = 0; //detection distance in depth layers (since this might be different than xy resolution)
		int buffer_cond = -9; //0= suitable habitat only, 1= presence of conspecifics, 2= density threshold of conspecifics
		float buffer_thresh = -9; //if buffer_cond=2, what is the dnesity threshold of conspecifics in ind/ha
		int buffer_capture = 0; //0=no, 1=yes, will it be captured by the buffer? settle automatically if within the buffer zone
	};// pop_sett;

	struct depth_info { //if depth preference is sex dependent, this struct will hold sex-specific parameters
		int Depthmin=-9, Depthmax=-9;
	};
	
	
	
	float **stage_mat; //this is to create a dynamic array to hold all the stage structured info 
					// needs to be dynamic because i dont know how many stages i have yet!
	
	struct stage_infob {

		int stage;
		int m_min_age = 0, f_min_age = 0, m_max_age, f_max_age;
		float m_fec = 0, f_fec = 0; //fecundity
		float m_surv = 0, f_surv = 0; //survival probability
		float m_trans = 0, f_trans = 0; //transition probability to next stage
		int m_min_depth, m_max_depth, f_min_depth, f_max_depth; //depth preferences

		struct emig_info s_emig;
		struct trans_info s_trans;
		struct sett_info s_sett;


		struct depth_info { //if depth preference is sex dependent, this struct will hold sex-specific parameters
			int Depthmin, Depthmax;
		};
	};
	
	//struct stage_info {
	//	int stage;
	//	int m_min_age = 0, f_min_age = 0, m_max_age, f_max_age;
	//	float m_fec = 0, f_fec = 0; //fecundity
	//	float m_surv = 0, f_surv = 0; //survival probability
	//	float m_trans = 0, f_trans = 0; //transition probability to next stage
	//	int m_min_depth, m_max_depth, f_min_depth, f_max_depth; //depth preferences
	//	//emigration
	//	float m_ep = 0, f_ep = 0; 
	//	float m_D0 = 0, m_alpha = 0, m_beta = 0; 
	//	float f_D0 = 0, f_alpha = 0, f_beta = 0;
	//	int ep_evo; float ep_mut_size, ep_mut, alpha_mut_size, beta_mut_size; //does emigration rate evolve? if yes (ep_evo==1), what is mutation size and rate?
	//	//transfer
	//	float f_rho, f_rho_sd, f_SL, f_comp_SL, f_SL_sd, f_mortprob, f_slope, f_inflpoint;
	//	float m_rho, m_rho_sd, m_SL, m_comp_SL, m_SL_sd, m_mortprob, m_slope, m_inflpoint;
	//	int dp, memory;
	//	float  m_min_active_time, f_min_active_time, m_min_active_size, f_min_active_size, f_active_sd, m_active_sd; //minimum size at which indivs switch to active movement
	//	//float m_m, m_m_sd, f_m, f_m_sd, f_b, m_b; //slope of growth function, intercept of growth function
	//	//int growth_evo; float growth_mut, m_mut_size, b_mut_size; //does the growth function evolve? if so, what is the mutation rate, and size for slope and intercept
	//	struct growth_info {
	//		int method = 0; //can be linear (0), gompertz (1) or von bertalanffy (2)
	//		int sex_dep = 0; //sex dependent growth? 0=no, 1=yes
	//		float m_m, m_m_sd, f_m, f_m_sd, f_b, m_b; //slope of growth function, intercept of growth function
	//		int growth_evo; float growth_mut, m_mut_size, b_mut_size; //does the growth function evolve? if so, what is the mutation rate, and size for slope and intercept
	//		float m_Linf = -9, m_Linf_sd = -9, m_G_K = -9, f_Linf = -9, f_Linf_sd = -9, f_G_K = -9, Linf_mut_size = -9, GK_mut_size = -9;
	//	}pop_grow;
	//	int active_evo; float active_mut, active_mut_size; //does this parameter evolve? if yes, at what probability and by how much does it mutate?
	//	int diel_vert = 0; //is behaviour size dependent or time dependent? is the growth function temperature dependent? do they undergo diel vertical migration?
	//		//if behaviour is size dependent, a growth function is applied during dispersal, otherwise time limits in days are applied
	//	int dv_active = 0; //does dvm apply to the active stage of dispersal?
	//	float min_dv_time, min_dv_size = -9; //if diel_vert==1, at what point do indivs start dvm? allows flexibility for passive->passive with dvm->active dispersal
	//	int down_bias = 0; //is there a downward bias in dispersal movements, ie tracking the seafloor
	//	int buoy_min, buoy_max, dv_range; //what are the upper and deep limits of buoyancy, and what is the diel vertical migration range?
	//	//this is calculated as buoy_min+dv_range and buoy_max-dv_range

	//	//settlement
	//	float f_S0 = -9, f_alphaS = -9, f_betaS = -9, f_sett_sd = -9;
	//	float m_S0 = -9, m_alphaS = -9, m_betaS = -9, m_sett_sd = -9;
	//	int S0_evo; float S0_mut, S0_mut_size, alphaS_mut_size, betaS_mut_size; //does settlement prob evolve? if so, at what prob and by how much? if settlement is density dependent, how much do alpha and beta change if a mutation occurs?

	//	int f_findmate, m_findmate;
	//	int f_pld = -9, f_pld_sd = -9, f_comp = -9, f_comp_size = -9, f_comp_sd = -9, m_pld = -9, m_pld_sd = -9, m_comp = -9, m_comp_size = -9, m_comp_sd = -9;
	//	int comp_evo; float comp_mut, comp_mut_size; //does competency time/size evolve? if so, at what probability and by how much?

	//	int buffer = -9; //0=no, yes=1
	//	int buffer_xycells = 0; //detection distance in number of cells
	//	int buffer_zlayers = 0; //detection distance in depth layers (since this might be different than xy resolution)
	//	int buffer_cond = -9; //0= suitable habitat only, 1= presence of conspecifics, 2= density threshold of conspecifics
	//	float buffer_thresh = -9; //if buffer_cond=2, what is the dnesity threshold of conspecifics in ind/ha
	//};

	Landscape* lands; //create a pointer to a landscape
	Model* mod; //create a pointer to the model
	vector<Subpopulation*> subpops; //vector of all active subpopulations
	//vector<stage_info*> stages; //vector of all stage structures
	vector<stage_infob*> stages_b;
	vector<depth_info*> sex_depths; //if depth preference is sex dependent BUT NOT stage-dependent

	struct matrix {
		int total; //this will count down until total=0, that is when dispersal phase ends? 
		vector<Individual*> matrix_indivs;
	}pop_matrix;

	void fill_stage_info(string transmatfile);
	void set_stages();
	void set_emig_info();
	void set_transfer_info();
	void set_settle_info();
	void set_evo_info();
	//stage_info* get_stage_info(int s); //s=required stage
	depth_info* get_depth_info(int s); //s= required sex, 0=male, 1=female
	//emig_info* get_pop_emig();
	//trans_info* get_pop_trans();
	//sett_info* get_pop_sett();
	void set_pop();
	float get_attribute(int attribute);
	void create_subpop(Cell* pcell); //overloaded create_subpop function so that i can pass in either a cell or a patch pointer
	void create_subpop(Patch* ppatch);
	void initialise();

	void one_year(ofstream& popfile, ofstream& indfile, ofstream& indmovfile, ofstream& rangefile, int rep, int year);
	//write outputs functions
	//float mean_evolution(string param);
	void write_popfile(ofstream&popfile, int rep, int year, int repseason);
	void write_rangefile(ofstream& rangefile, int rep, int year, int repseason);
	void write_indivfile(ofstream& indfile, int rep, int year, int repseason);
	void write_indmovfile(ofstream& indmovfile, int rep, int year, int repseason, bool lastyear);

	/*void write_popfile(int rep, int year, int repseason);
	void write_rangefile( int rep, int year, int repseason);
	void write_indivfile(int rep, int year, int repseason);
	void write_indmovfile(int rep, int year, int repseason, bool lastyear);*/

	//reset for next rep
	void reset();
	//delete subpops
	void delete_subpops();
};

class Subpopulation:public Population {
private:
	
	
public:
	Subpopulation();
	~Subpopulation();
	//Subpopulation(Cell* c, Model* m, Landscape* la, pop_info* pd, matrix* pm, emig_info* e, trans_info* t, sett_info*s, vector<stage_info*> st);
	//Subpopulation(Patch* p, Model* m, Landscape* la, pop_info* pd, matrix* pm, emig_info* e, trans_info* t, sett_info*s, vector<stage_info*> st);
	Subpopulation(Cell* c, Model* m, Landscape* la, pop_info* pd, matrix* pm,  vector<stage_infob*> st);
	Subpopulation(Patch* p, Model* m, Landscape* la, pop_info* pd, matrix* pm,  vector<stage_infob*> st);

	//inherited objects:
	//Landscape* lands;
	//Model* mod;

	pop_info * pop_dyn; //pointer to population information structure
	matrix* pop_matrix; //pointer to the matrix structure which stores indivs that are dispersing
	//emig_info* pop_einfo; //pointer to an emig_info structure containing the population emigration information
	//trans_info* pop_tinfo; //pointer to the struct containing transfer information
	//sett_info* pop_sinfo; //pointer to the struct containing settlement information
	//vector<stage_info*> pop_stageinfo;
	
	Patch* patch_num;
	Cell* cell_num;
	int size = 0; //number of individuals
	int nmales = 0, nfemales = 0; //these are the numbers of SEXUALLY MATURE males and females

	vector<Individual*> rep_females; //reproductively mature females
	vector<Individual*> rep_males;  //reproductively mature males

	vector<Individual*> indivs; //all individuals
	vector<Individual*> juv_indivs; //offspring produced in reproduction in a NON STAGE STRUCTURED MODEL. these will then be sent onto either transit_indivs or surv_indivs vectors
	vector<Individual*> surv_indivs; //surviving individuals after all processes. this vector will overwrite the indivs vector for the next generation
	vector<stage_infob*> pop_stageinfo;

	void init_age(Individual* pind, int s);
	void init_size(Individual* pind);
	int assign_sex();
	void init_position(Individual* pind);
	int init_number();
	void init_indivs();

	void leave(Individual* pind);
	void enter(Individual* pind);
	Population::stage_infob* get_stage_info(int s);
	
	void reproduction();
	void juv_release(Individual* pind);
	void emigration();
	void survival();
	void development(bool aging);
	void clean_up();

	vector<float> summarise();
	vector<float> mean_evolution(string param, bool first);
	void delete_indivs();

};

class Individual:public Subpopulation {
private:

public:
	Individual(); //default constructor
	~Individual();
	//Individual(Model* mo, pop_info* p, Landscape* la,  int x, int IDnum); //non stage structured individual constructor
	Individual(Model* mo, pop_info* p, Landscape* la,  stage_infob* st, int x, int IDnum, bool invent, vector<stage_infob*>* emigst); //stage structured constructor
	//Individual(Model* mo, pop_info* p, Landscape* la, emig_info* e, trans_info* t, sett_info* s, int x, int IDnum); //non stage structured individual constructor
	//Individual(Model* mo, pop_info* p, Landscape* la, emig_info* e, trans_info* t, sett_info* s, stage_infob* st, int x, int IDnum); //stage structured constructor

	//inherited objects:
	//pop_info * pop_dyn
	//Landscape* lands;
	//emig_info* pop_einfo; //pointer to an emig_info structure containing the population emigration information
	//trans_info* pop_tinfo; //pointer to the struct containing transfer information
	//matrix* pop_matrix; //pointer to the matrix struct to add in-transit individuals to matrix
	//Model* mod;

	int ID;
	int age = 0; //in years?
		//bool rep; //has it reproduced?
	float size = 0, size_at_birth=0; //in mm
	float m, b; //the slope and intercept of the linear growth function, applied during transfer. this can be temp dependent if temp info has been given
	bool sett; //has it settled?
	bool emig; //has it emigrated?
	bool dvm; //is it old/big enough to undergo diel vertical migration?
	bool dead; //has it run out of time to settle and died?
	int gen = 0; //what generation is it in?
	int sex; //0=male, 1=female
	bool skip_breeding;
	int rep_int; //this will increment every time an individual skips a reproductive event. when rep_int > propulation RepInt, skip_breeding=false
	int status=0; //0=just born, 1=disperser, 2=local recruit, 3=settled in suitable habitat, 4=forced to settle on unsuitable habitat and died,
					// 5=died by dispersal-related mortality, 6= died by annual mortality, 7=died by exceeding max age, 8=died by being absorbed, 
					// 9=forced settlement(passive) on suitable habitat, 10= ran out of time; 11= fishing mortality, 12= settlement by buffer_capture


	stage_infob* pstage; //pointer to a stage structure containing transition matrix information
	vector<stage_infob*>* subpop_stages;
	
	struct disp_info { //for the individual to store its own dispersal information, what it actually expresses!
		struct parent_info {
			//if it's an asexual model, they will inherit only one allele
			//if it's a sexual model, they will inherit 2 alleles: one from the female and one from the male
			//if it's a sexual model where the parameter is sex-dependent, then they will inherit 4 alleles: female allele for female offspring, ...male offspring, male allele for female offspring,...male offspring

			//top row of these arrays is always female allele for female, male allele for female
			//bottom row is female allele for male, male allele for male
			//NB: if the model is sex independent, then only female allele for female and male allele for female is used. so it's a 1x2 matri 
			//NB if it's an asexual model then only element [0][0] will be initialised to anything other than -9
			//so the full matrix looks like this:
			//			Mother		Father
			// female	
			// male
			float size_at_birth[2][2] = { -9,-9,-9,-9 }, emig_prob[2][2] = { -9,-9,-9,-9 }, D0[2][2] = { -9,-9,-9,-9 }, alpha[2][2] = { -9,-9,-9,-9 }, beta[2][2] = { -9,-9,-9,-9 },
				rho[2][2] = { -9,-9,-9,-9 }, rho_m[2][2] = { -9,-9,-9,-9 }, SL[2][2] = { -9,-9,-9,-9 }, SL_m[2][2] = { -9,-9,-9,-9 }, active_SL[2][2] = { -9,-9,-9,-9 },
				min_active_size[2][2] = { -9,-9,-9,-9 }, min_active_time[2][2] = { -9,-9,-9,-9 }, min_dv_size[2][2] = { -9,-9,-9,-9 }, min_dv_time[2][2] = { -9,-9,-9,-9 },
				m[2][2] = { -9,-9,-9,-9 }, Linf[2][2] = { -9,-9,-9,-9 }, G_K[2][2] = { -9,-9,-9,-9 }, Ti[2][2] = { -9,-9,-9,-9 },
				S0[2][2] = { -9,-9,-9,-9 }, alphaS[2][2] = { -9,-9,-9,-9 }, betaS[2][2] = { -9,-9,-9,-9 }, comp_size[2][2] = { -9,-9,-9,-9 }, comp_time[2][2] = { -9,-9,-9,-9 };
				//pld[2][2] = { -9,-9,-9,-9 };

		};// parent;
		parent_info* parent = nullptr; //if it's an asexual model, i dont need to initialise one
		void assign_parent() {
			parent = new parent_info;
		}
		int mode = 0; //0=passive, 1=hybrid, 2=active;
		int phase = 0; //0=passive, 1=active;
		//emigration
		float emig_prob = -9, D0 = -9, alpha = -9, beta = -9; 
		/*int ep_evo; float ep_mut_size, ep_mut, alpha_mut_size, beta_mut_size;*/
		//transfer
		float rho, rho_m=-9, rho_sd=-9, SL=0, SL_m=-9, SL_sd=-9, mortprob=-9, slope=-9, inflpt=-9; //rho_m and SL_m are the slopes for size-dependent parameters, these are randomly sampled so there are weak and strong swimmers and those that choose to follow current and dont
		//rho_sd and SL_sd are for inter indiv var. if ==-9, then there is no iiv 
		float buoy_min=-9, buoy_max=-9; //what are the boundaries for larval buoyancy? if this isnt applicable, value stays -9 and indivs can move to all boundaries of the seascape
		int diel_vert=0; //do individuals undergo diel vertical migration? 0=no, 1=yes, default=no
		int dv_range = 0; //if diel_vert=1, in what range do the indivs gather? ie if dv_range=5, and buoy_min=3, then indivs would gather between 3 and 8 m during the night
		int dv_active = 0; //does dvm apply to the active stage of dispersal? 0=no, 1=yes
		float min_dv_time, min_dv_size = -9, min_dv_sd=-9; //if diel_vert==1, at what point do indivs start dvm? allows flexibility for passive->passive with dvm->active dispersal. is ther iiv?
		float active_SL = -9; //if indiv is a hybrid (passive then active) disperser, this is the SL once it passes its competency period
		int down_bias = 0; //is there a downward bias during dispersal movements? this would simulate tracking the seafloor
		int min_active_time=0; //this is the minimum number of days after which indiv can disperse actively
		float min_active_size=0; //if behaviour is size_dependent, this is the minimum size after which indivs can disperse actively
		float active_sd = -9;
		/*int active_evo; float active_mut, active_mut_size;*/ //does this parameter evolve? if yes, at what probability and by how much does it mutate?
		//float m=-9, b=-9; //if growth function is applied, m= slope of linear function, b=intercept of linear function. this is only applied if behaviour is size-dependent
		//float m_sd = -9; //if there is inter indiv var in growth rate, this is the sd
		//int growth_evo; float growth_mut, m_mut_size, b_mut_size; //does the growth function evolve? if so, what is the mutation rate, and size for slope and intercept
		struct growth_info {
			int method = 0; //can be linear (0), gompertz (1) or von bertalanffy (2)
			int sex_dep = 0, stage_dep=0; //sex dependent growth? 0=no, 1=yes
			float m = -9, b = -9; //if growth function is applied, m= slope of linear function, b=intercept of linear function. this is only applied if behaviour is size-dependent
			float m_sd = -9; //if there is inter indiv var in growth rate, this is the sd
			/*int growth_evo; float growth_mut, m_mut_size, b_mut_size; *///does the growth function evolve? if so, what is the mutation rate, and size for slope and intercept
			float Linf = -9, Linf_sd = -9, G_K = -9, Linf_mut_size = -9, GK_mut_size = -9, Ti_sd = -9;
			int Ti = -9; //this is the inflection point of the gompertz function, equates to earliest settlement date 
		}grow_info;
		//float temp_SL = 0; //this is for when it detects a new cell, so it can adjust its speed when moving vertically
		//settlement
		float S0 = -9, alphaS = -9, betaS = -9, sett_sd = -9;
		/*int S0_evo; float S0_mut, S0_mut_size, alphaS_mut_size, betaS_mut_size;*/ //does settlement prob evolve? if so, at what prob and by how much? if settlement is density dependent, how much do alpha and beta change if a mutation occurs?
		int pld = -9, pld_sd = -9, comp_time = -9, comp_size = -9, comp_sd = -9; //comp_time here is minimum time that indiv is able to settle; if behaviour is size-dependent, then use comp_size (min size)
		/*int comp_evo; float comp_mut, comp_mut_size;*/ //does competency time/size evolve? if so, at what probability and by how much?
		//pld is the maximum time they can disperse for. if still dispersing after this time, they die
		int findmate;
		bool comp = false; //is it competent for settlement? this can be size or time dependent
		//bool buffer_capture = false; //are they considered settled if they are within the buffer distance of suitable habitat? if true, movement ends and settlement is immediate
			//this would apply to both active and passive but only after competency
		//if settlement has a detection buffer
		int goal_bias = 0; //this is if the individual is within a buffer cell, therefore has a goal bias towards suitable habitat
		float goal_rho = .99; //if goal_bias==1, use this rho instead of user-defined, so that indiv moves towards goal
		Cell* goal_focus = nullptr;
		vector<int> tested_sett;
		int memory, dp;
		float disp_time = 0; //this counts in numbers. individuals will attempt dispersal until this number meets their limitation, if applicable
		float disp_distance = 0; //this keeps track of the distance travelled (2D, in x,y plane) in case of distance-dependent mortality
		float per_step_disp_dist = 0;
		bool stuck = false;
		int status = 0;
		//status=1 : individual has been absorbed, dead
		//status=2 : an active individual has met a reflective boundary, resample direction
		//status=3 : an active individual has met a suitable habitat cell, check settlement
		//status=4 : a passive individual has met a suitable habitat cell, forced settlement
		//status=5 : a passive individual has met an unsuitable habitat cell, dead
		//status=6 : a passive individual has hit the surface, continue with z coordinate=0
		//status=7 : valid move
		//status 8 : considering buffer capture
	}dinfo;

	//emig_info* pop_einfo; //pointer to an emig_info structure containing the population emigration information
	//trans_info* pop_tinfo; //pointer to the struct containing transfer information
	//sett_info* pop_sinfo; //pointer to the struct containing settlement information
	int natal_patch;
	int natal_cell;
	Patch* current_patch=nullptr;
	Cell* current_cell=nullptr;
	vector<float> m_x_movements;
	vector<float> m_y_movements;
	vector<float> m_z_movements;
	struct locn { float x, y, z; } current_pos;
	std::queue<locn> memory;
	//vector<locn> movements;
	//vector<float> current_pos; //m_x, m_y, m_z

	void disp_gene_expression();
	void parameter_iiv();
	//void inheritance(Individual* parent);
	void inheritance(Individual* mother, Individual* father);
	bool check_repro();
	float get_localK();
	int get_localN();
	int reproduction();

	void set_emig_info(bool invent);
	bool check_emig();
	void calc_SL_rho();
	void set_growth(bool invent);
	void set_trans_info(bool invent);
	//void update_pos(int new_x, int new_y, int new_z);
	//bool is_buffer();
	//void assess_destination(float x, float y, float z);
	//void transfer_pos();
	void set_settle_info(bool invent);
	bool check_settle(Cell* cell_sett, float pot_x, float pot_y, float pot_z);
	//bool check_settle();
	bool transfer_to_settle();
	void dispersal_mort();
	bool survival();
	bool fishing_mort();
	void stage_development();
	void larval_growth();

	//mat333 calc_dp();
	//vector<float> assess_move();

	//movement model functions
	float calc_cell_mag(int row, int col, int layer, float x, float y, float z);
	void update_pos(float new_x, float new_y, float new_z, int row, int col, int layer, bool mem, bool mem_empty, bool c_cell, bool c_pos);
	vector<int> calc_r_c_l(float x, float y, float z);
	float calc_distance(float x, float y, float z);
	vector<float> evade(float phi);
	bool check_cell(int row, int col, int layer, float x, float y, float z);
	//void passive_traverse(float x, float y, float z);
	//void traverse1(float x, float y, float z);
	void traverse2(float x, float y, float z);
	void traverse3(float x, float y, float z);
	//void traverse(float x, float y, float z);
	void leave_natal();
	vector<float> calc_move();
	vector<float> calc_move2();
	bool buffer_capture();

};

class Landscape {
private:
	int nrows;
	int ncols;
	short nlayers;
	int ncells;
	short nhabitats;
	short npatches;
	float xllcorner;
	float yllcorner;
	float resolution;
	float no_data;
	float x_min, x_max;
	float y_min, y_max;
	float z_min, z_max;
	short dint;


public:
	Landscape();
	~Landscape();
	bool cell_based;
	bool patch_based;
	bool sp_dist;
	short dynamic = 0; //is landscape dynamic? 0=no, 1=yes
	int patch_extinct = 0, p_ext_method = -9, p_ext_int = -9, p_ext_burnin = -9; //is there regular local patch extinction? which method: 0=random, 1=specific; at what interval? how many years between events?
	float p_ext_p_prop = -9, p_ext_i_prop = -9; //if p_ext_method==0, what proportion of patches?
	int p_ext_patch = -9; //if p_ext_method==1, which patch?
	


	Cell ****raster_3d = 0;

	Model* mod;
	vector<Patch*> patch_vector; //vector of pointers to patches of length npatches
	vector<Cell*> suitable_cells; //based on habitat type
	vector<Cell*> suitable_depth; //based on habitat type and depth preferences of species
	vector<Cell*> spdist_cells; //based on species distribution data
	vector<float> car_caps; //vector of carrying capacities corresponding to each habitat type, in inds/ha
	
	struct dynland { //this is a struct for dynamic landscapes, each struct will contain the files for that year's new landscape data
		vector<string> u_files; //this vector will be nlayers long
		vector<string> v_files;
		vector<string> w_files;
		//vector<string> hab_files;
		//vector<string> p_files;
	};
	struct dyn_info { //this holds all the info needed for dynamic landscapes
		string dyn_files;
		short num=0; //how many sets of new landscape is being read in?
		short burn_in=0; //after how many years does first new input start?
		short dyn_interval=0; //how many years between consecutive inputs?
		short dyn_hab = 0, dyn_hydro = 0, dyn_patches = 0; //what is changing? habitats? hydrodynamics? patches?
		
		vector<dynland> dyn_land_changes; //this is a vector of structs, .size() will be number of changes made
	}dynamic_info;
	
	//get functions
	float get_land_att(int x);
	
	//other functions
	void create_raster();
	void create_patch_vector();
	void delete_raster();
	void delete_patch_vector();
	void check_landscape();
	//set the x and y metre coordinates
	void fill_md();
	//set values for the other attributes in the cell
	//void set_value(ifstream& file, int layer, int attribute);
	//determine if cell is nodata
	void set_NoData(int layer, string input_loc, string filename);
	//call other functions to fill the landscape
	void fill_cells(string input_loc, vector<string> files_to_fill);
	dynland declare_dynland(string u, string v, string w, string h, string p);
	void read_dyn();
	void do_dyn(int change_num);
	//reads control file and calls above functions to fill cells with information
	void get_info();  

	//void fill_neigh();
};

class Patch /*:public Landscape*/ {
public:
	Patch();
	~Patch();
	int patch_ID=-9;
	bool in_spdist;
	//int current_density=0; //how many individuals are present?
	//vector<int> densities; //record of the densities over time
	struct limits {
		int xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0; //these are in rows,columns, layers. but remember that x=columns and y=rows
	} patchlimits; //initialise a limits struct called patchlimits

	struct K_info{
		float K = 0; //carrying capacity in inds/ha
		int size_cells=0; // size measured in cells. this will affect recruitment
		float size_ha=0; //how big is it in ha? float because it might be 4.53 hectares
		int max_inds=0; //what is the maximum number of individuals this area can hold?
	};

	vector<Cell*> included_cells;
	vector<Cell*> indepth_cells;
	vector<Cell*> in_spdist_cells;
	Subpopulation* p_subpop=nullptr;

	K_info patch_info, indepth_info; //the patch might be bigger than the depth range
	void set_limits(int row, int column, int depth, int dint);
	bool check_limits(int minX, int maxX, int minY, int maxY, int minZ, int maxZ);
	void calc_K(string whichstruct);

};

class Cell/* :public Landscape */{
private:

public:
	Cell();
	~Cell();
	Cell(Model* modp);
	Model* mod;
	//Landscape* lands = nullptr;
	int cell_ID;
	float max_inds; //K*ha
	//int current_density;
	short habtype,spdist, r, c, l;//col, row, layer
	short bufferpos = -9999; //if 0, then it is not within the user-define detection distance of suitable habitat
	short bufferneg = -9999;
	//int cost;
	float res,u, v, w; //east, north, upwards;
	//short 
	//float 
	float cell_speed, cell_angle, cell_temp; // in degrees C
	float min_depth, max_depth, seafloor_depth;
	float K;  //in case of a cell-based model, each cell will have a carrying capacity, in inds/ha
	float ha; //size in ha
	vector<float> midpoint;
	
	vector<Cell*> buffer_focus; //if buffer==1, then this will hold the pointers to the focus cells that are suitable habitat in the vicinity. if there is more than one, an indiv can randomly choose 
	//vector<Cell*> buffer_natal; 
	/*struct open_water_neigh {
		vector<Cell*> down_neigh;
		vector<Cell*> forward_neigh;
		vector<Cell*> up_neigh;
	}open_water_neigh;*/
	//Cell* ow_neigh[3][3][3];
	//vector<Cell*> open_water_neigh;
	struct mat_index { int row, col, layer; };
	vector<mat_index> possibilities;
	struct four_corners {
		vector<float> llc, ulc, lrc, urc; //lower left corner, upper left corner, etc
	} cell_corners;

	//mat333 costs, current;
	
	Subpopulation* c_subpop=nullptr;
	Patch* patch=nullptr;

	bool check_limits(int minX, int maxX, int minY, int maxY, int minZ, int maxZ);
	void calc_midpoint();
	void calc_maxinds();
	void assess_pos_buffer(int detect_dist, int det_layers, Landscape* lands);
	void find_ow_neigh(Landscape* lands);
	//void assess_neg_buffer(int detect_dist, int det_layers, int natal, Landscape* lands);

	//float find_mean_cost(int minrow, int maxrow, int mincol, int maxcol, int minlay, int maxlay, Landscape* lands);
	//vector<mat333> calc_weightings(int xy_pr, int z_pr, float rho, Landscape* lands); //direction indiv was travelling in, directional persistence, and perceptual range
};

class Model {

public:
	Model();
	~Model();
	short nsims=-9; //how many sims in this batch
	short nreps=-9; //how many replicates of each sim
	int nyears=-9; //how many years should the model run for
	short batchnum=-9;
	short patchmodel; //0=no, 1=yes
	short nhabs=-9;
	short speciesdist; //0=not from species distribution, 1=from species distribution
	int distres; //what resolution is the species distribution?
	short stagestruct; //0=not stagestructured, 1=stage structured
	short stages=1; //how many stages are there
	short repro; //what type of reproduction? 0=asexual, only females modelled. 1=sexual, both males and females
	short temp_dep = 0; //0=no, 1=yes, is this model influenced by temperature? if yes, then it will be implemented into the growth function
	short size_or_time = 1; //0=size, 1=time. is behaviour during dispersal size or time dependent? if size, then a growth function is applied and behaviour changes are defined by minimum size. 
		//otherwise time limits (in days) are used
	short dynamic = 0; //are landscapes dynamic? 0=no, 1=yes
	string grow_method;
	float propmales;
	
	//info files
	string wrkdir;
	string paramfile;
	string landfile;
	string dynfile;
	string initfile;
	string popfile;
	string stagestructfile;
	string growthfile;
	string emigfile;
	string transfile;
	string setfile;
	string evofile;

	//output options
	short Outpop, Outind, Outindmov, Outrange, outpop_int, outind_int, outrange_int; //output population, individual and range files? at what intervals?
	

	//initialisation parameters
	short seedtype; //0=free initialisation, 1=from species distribution, 2=specific patch
	short freetype; // 0= random with a proportion of suitable cells, 1= all suitable cells
	float freeprop; //if freetype=0, what proportion of suitable cells
	int minX, maxX, minY, maxY, minZ, maxZ;
	short sptype; float spprop; //0=random from a proportion, 1=all suitable, proportion if sptype=0
	short initdens; //0=at K, 1= at half K, 2= specific density
	float indsha; //what density if initdens=2
	short whichpatch; //if seedtype=2, which patch to initialise
	short initage;
	vector<float> stage_props = { 0 };



	void read_control();
	void read_init();
};


#endif CLASSES_H