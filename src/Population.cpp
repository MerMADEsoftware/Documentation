#include "classes.h"

Population::Population() {

}

Population::~Population()
{
}

/////setting functions//////////
void Population::fill_stage_info(string transmatfile) { //this fills reproductive and transition information, max age, max depth preferences
	int nrow, ncol;

	//how many stages structures do we need
	for (int s = 0; s < mod->stages; s++) {
		stage_infob* stage = new stage_infob;
		stage->stage = s;
		stages_b.push_back(stage);
	}

	//the number of columns depends on whether we are looking at females and males or just females
	if (mod->repro > 0) { nrow = (2 * mod->stages) - 1; ncol = (2 * mod->stages) + 3; }
	else { nrow = mod->stages; ncol = mod->stages + 3; }

	//create a matrix the right size to hold all the information
	stage_mat = new float*[nrow];
	for (int i = 0; i < nrow; i++) {
		stage_mat[i] = new float[ncol];
	}

	//open the transition matrix file
	ifstream transmat(mod->wrkdir + "/Inputs/" + transmatfile);
	string header;
	//copy the transition matrix into the dynamic array
	for (int h = 0; h < ncol + 1; h++) { transmat >> header; } //+1 is because i have a column header "Transition" in the way
	for (int s = 0; s < nrow; s++) {
		transmat >> header;
		for (int c = 0; c < ncol; c++) {
			float val;
			transmat >> val;
			stage_mat[s][c]=val;
		}
	}


	//assign the values to the appropriate stage 
	ncol -= 1; nrow -= 1; //subtract 1 from ncol and nrow so i can index stage_mat properly (make ncol=last column of stage_mat)
	for (int s = 0; s < mod->stages; s++) {
		int row, col;
		if (mod->repro > 0) { //if it's a sexual model, with males and females
			if (s == 0) { //because the only variables that have value for stage 0 is transition probability, age and depth
				stages_b[s]->m_trans = 1; stages_b[s]->f_trans = 1;  //because they can't stay in stage 0, transition prob =1
				stages_b[s]->m_surv = stage_mat[1][0]; //the value given in the matrix is surv*transition but since transition=1, surv= value given
				stages_b[s]->f_surv = stage_mat[2][1];
				stages_b[s]->stage = 0;
				stages_b[s]->m_max_age = 1; stages_b[s]->f_max_age = 1; //stage 0 dont remain stage 0 so max age is 1 timestep
				if (pop_dyn.DepthStage == 1) { stages_b[s]->m_min_depth = stage_mat[0][ncol - 1]; stages_b[s]->m_max_depth = stage_mat[0][ncol]; } //if depth pref is stage dependent
			}
			else {
				row = (s * 2) - 1; col = s * 2;
				stages_b[s]->stage = s;
				stages_b[s]->m_fec = stage_mat[0][col];
				stages_b[s]->f_fec = stage_mat[0][col + 1];

				stages_b[s]->m_min_age = stage_mat[row][ncol - 2];
				stages_b[s]->f_min_age = stage_mat[row + 1][ncol - 2];
				if (pop_dyn.DepthStage == 1) { //if depth preference is stage dependent
					stages_b[s]->m_min_depth = stage_mat[row][ncol - 1]; stages_b[s]->m_max_depth = stage_mat[row][ncol];
					stages_b[s]->f_min_depth = stage_mat[row + 1][ncol - 1]; stages_b[s]->f_max_depth = stage_mat[row + 1][ncol];
				}
				if (s == mod->stages - 1) { //if it's the oldest stage
					stages_b[s]->m_max_age = pop_dyn.max_age; //of species
					stages_b[s]->f_max_age = pop_dyn.max_age; //of species
					//stages[s]->m_surv = stage_mat[row][col]; // because they only have one number coding for survival, i dont need to add this to anything like the other stages
					//stages[s]->f_surv = stage_mat[row + 1][col + 1];
					stages_b[s]->m_surv = stage_mat[row][col];
					stages_b[s]->m_trans = 0;
					stages_b[s]->f_surv = stage_mat[row + 1][col + 1];
					stages_b[s]->f_trans = 0;
				}
				else {
					stages_b[s]->m_max_age = stage_mat[row + 2][ncol - 2];//take the min age of the next stage
					stages_b[s]->f_max_age = stage_mat[row + 3][ncol - 2];
					//survival: need to add transition prob + retention prob together (even if retention prob might be 0)
					stages_b[s]->m_surv = stage_mat[row][col] + stage_mat[row + 2][col];
					stages_b[s]->f_surv = stage_mat[row + 1][col + 1] + stage_mat[row + 3][col + 1];
					//transition prob needs to take into account survival prob and transition prob, so it is not the value that is given! need to calculate it
					stages_b[s]->m_trans = (stage_mat[row + 2][col]) / stages_b[s]->m_surv;
					stages_b[s]->f_trans = (stage_mat[row + 3][col + 1]) / stages_b[s]->f_surv;

				}
			}
		}
		else { //if it's a female only model
			stages_b[s]->stage = s;
			stages_b[s]->f_fec = stage_mat[0][s];
			stages_b[s]->f_min_age = stage_mat[s][ncol - 2];
			if (pop_dyn.DepthStage == 1) { //if depth preference is stage dependent
				stages_b[s]->f_min_depth = stage_mat[s][ncol - 1];
				stages_b[s]->f_max_depth = stage_mat[s][ncol];
			}
			if (s == 0) { //if it's the first stage
				stages_b[s]->f_max_age = 0; //can't stay in stage 0 for more than 1 reproductive season
				stages_b[s]->f_surv = stage_mat[s + 1][0];
				stages_b[s]->f_trans = 1; //they have to transition
			}
			else if (s == mod->stages - 1) {  //if its the last stage
				stages_b[s]->f_max_age = pop_dyn.max_age; //of species
				stages_b[s]->f_surv = stage_mat[s][s]; //only has one value coding for survival so don't need to add anything to it
			}
			else {
				stages_b[s]->f_max_age = stage_mat[s + 1][ncol - 2];//of next stage
				stages_b[s]->f_surv = stage_mat[s][s] + stage_mat[s + 1][s]; //need to add transition prob + retention prob together (even if retention prob might be 0)
				stages_b[s]->f_trans = (stage_mat[s + 1][s]) / stages_b[s]->f_surv; //last stage doesnt have transition prob so this only applies here

			}
		}

	}
}

//void Population::fill_stage_info(string transmatfile) { //this fills reproductive and transition information, max age, max depth preferences
//	int nrow, ncol;
//
//	//how many stages structures do we need
//	for (int s = 0; s < mod->stages; s++) {
//		stage_info* stage = new stage_info;
//		stages.push_back(stage);
//	}
//
//	//the number of columns depends on whether we are looking at females and males or just females
//	if (mod->repro > 0) { nrow = (2 * mod->stages) - 1; ncol = (2 * mod->stages) + 3; }
//	else { nrow = mod->stages; ncol = mod->stages + 3; }
//	
//	//create a matrix the right size to hold all the information
//	stage_mat = new float*[nrow];
//	for (int i = 0; i < nrow; i++) {
//		stage_mat[i] = new float [ncol];
//	}
//
//	//open the transition matrix file
//	ifstream transmat(mod->wrkdir + "/Inputs/" + transmatfile);
//	string header;
//	//copy the transition matrix into the dynamic array
//	for (int h = 0; h < ncol + 1; h++) {transmat >> header;} //+1 is because i have a column header "Transition" in the way
//	for (int s = 0; s < nrow; s++) { 
//		transmat >> header;
//		for (int c = 0; c < ncol; c++) {
//			transmat >> stage_mat[s][c];
//		}
//	}
//
//	//assign the values to the appropriate stage 
//	ncol -= 1; nrow -= 1; //subtract 1 from ncol and nrow so i can index stage_mat properly (make ncol=last column of stage_mat)
//	for (int s = 0; s < mod->stages; s++) {
//		int row, col;
//		if (mod->repro > 0) { //if it's a sexual model, with males and females
//			if (s == 0) { //because the only variables that have value for stage 0 is transition probability, age and depth
//				stages[s]->m_trans = 1; stages[s]->f_trans = 1;  //because they can't stay in stage 0, transition prob =1
//				stages[s]->m_surv = stage_mat[1][0]; //the value given in the matrix is surv*transition but since transition=1, surv= value given
//				stages[s]->f_surv = stage_mat[2][1];
//				stages[s]->stage = 0;
//				stages[s]->m_max_age = 1; stages[s]->f_max_age = 1; //stage 0 dont remain stage 0 so max age is 1 timestep
//				if (pop_dyn.DepthStage == 1) { stages[s]->m_min_depth = stage_mat[0][ncol - 1]; stages[s]->m_max_depth = stage_mat[0][ncol]; } //if depth pref is stage dependent
//			}
//			else {
//				row = (s * 2) - 1; col = s * 2;
//				stages[s]->stage = s;
//				stages[s]->m_fec = stage_mat[0][col];
//				stages[s]->f_fec = stage_mat[0][col + 1];
//				
//				stages[s]->m_min_age = stage_mat[row][ncol-2];
//				stages[s]->f_min_age = stage_mat[row + 1][ncol-2];
//				if (pop_dyn.DepthStage == 1) { //if depth preference is stage dependent
//					stages[s]->m_min_depth = stage_mat[row][ncol - 1]; stages[s]->m_max_depth = stage_mat[row][ncol];
//					stages[s]->f_min_depth = stage_mat[row+1][ncol - 1]; stages[s]->f_max_depth = stage_mat[row+1][ncol];
//				}
//				if (s == mod->stages-1) { //if it's the oldest stage
//					stages[s]->m_max_age = pop_dyn.max_age; //of species
//					stages[s]->f_max_age = pop_dyn.max_age; //of species
//					//stages[s]->m_surv = stage_mat[row][col]; // because they only have one number coding for survival, i dont need to add this to anything like the other stages
//					//stages[s]->f_surv = stage_mat[row + 1][col + 1];
//					stages[s]->m_surv = stage_mat[row][col]; 
//					stages[s]->m_trans = 0;
//					stages[s]->f_surv = stage_mat[row + 1][col + 1];
//					stages[s]->f_trans = 0;
//				}
//				else{ 
//					stages[s]->m_max_age = stage_mat[row + 2][ncol - 2];//take the min age of the next stage
//					stages[s]->f_max_age = stage_mat[row + 3][ncol - 2];
//					//survival: need to add transition prob + retention prob together (even if retention prob might be 0)
//					stages[s]->m_surv = stage_mat[row][col] + stage_mat[row + 2][col]; 
//					stages[s]->f_surv = stage_mat[row + 1][col + 1] + stage_mat[row + 3][col + 1];
//					//transition prob needs to take into account survival prob and transition prob, so it is not the value that is given! need to calculate it
//					stages[s]->m_trans = (stage_mat[row + 2][col]) / stages[s]->m_surv;
//					stages[s]->f_trans = (stage_mat[row + 3][col + 1])/stages[s]->f_surv;
//					
//				}
//
//				
//			}
//		}
//		else { //if it's a female only model
//			stages[s]->stage = s;
//			stages[s]->f_fec = stage_mat[0][s];
//			stages[s]->f_min_age = stage_mat[s][ncol-2];
//			if (pop_dyn.DepthStage == 1) { //if depth preference is stage dependent
//				stages[s]->f_min_depth = stage_mat[s][ncol - 1];
//				stages[s]->f_max_depth = stage_mat[s][ncol];
//			}
//			if (s == 0) { //if it's the first stage
//				stages[s]->f_max_age = 0; //can't stay in stage 0 for more than 1 reproductive season
//				stages[s]->f_surv = stage_mat[s + 1][0];
//				stages[s]->f_trans = 1; //they have to transition
//			}
//			else if (s == mod->stages - 1) {  //if its the last stage
//				stages[s]->f_max_age = pop_dyn.max_age; //of species
//				stages[s]->f_surv = stage_mat[s][s]; //only has one value coding for survival so don't need to add anything to it
//			}
//			else { 
//				stages[s]->f_max_age = stage_mat[s + 1][ncol-2];//of next stage
//				stages[s]->f_surv = stage_mat[s][s] + stage_mat[s + 1][s]; //need to add transition prob + retention prob together (even if retention prob might be 0)
//				stages[s]->f_trans = (stage_mat[s + 1][s])/stages[s]->f_surv; //last stage doesnt have transition prob so this only applies here
//				
//			} 
//		}
//
//	}
//}

void Population::set_stages() {
	if (mod->stages > 1 && mod->stagestruct == 1) {
		int sim;
		ifstream structfile(mod->wrkdir + "/Inputs/" + mod->stagestructfile);
		string header, transmatfile;
		int nheaders = 3; //change this when i start doing simulations!
		for (int h = 0; h < nheaders; h++) {
			structfile >> header;
		}
		structfile >> sim >> pop_dyn.DepthStage >> transmatfile;
		fill_stage_info(transmatfile);
	}
}

void Population::set_emig_info() { //this is only written for one sim at a time right now!
	int sim, holder;
	string header;
	ifstream emigfile(mod->wrkdir + "/Inputs/" + mod->emigfile);
	int nheaders = 13;
	for (int h = 0; h < nheaders; h++) {
		emigfile >> header;
	}
	emigfile >> sim >> stages_b[0]->s_emig.densdep >> stages_b[0]->s_emig.stagedep >> stages_b[0]->s_emig.sexdep >> stages_b[0]->s_emig.indvar >> stages_b[0]->s_emig.iiv_sd >> stages_b[0]->s_emig.emigstage >> holder >> holder;
	//holder holder is because ill step through the stages and sexes in the for loops
	//cout << "in set_emig_info, ";
	for (int s = 0; s < mod->stages; s++) {
		//cout << "doing stage " << s ;
		if (s != 0) {
			stages_b[s]->s_emig.densdep = stages_b[0]->s_emig.densdep; stages_b[s]->s_emig.stagedep = stages_b[0]->s_emig.stagedep; 
			stages_b[s]->s_emig.sexdep = stages_b[0]->s_emig.sexdep; stages_b[s]->s_emig.indvar = stages_b[0]->s_emig.indvar; 
			stages_b[s]->s_emig.iiv_sd = stages_b[0]->s_emig.iiv_sd; stages_b[s]->s_emig.emigstage = stages_b[0]->s_emig.emigstage;
			//holder holder is because ill step through the stages and sexes in the for loops
		}
		if (stages_b[0]->s_emig.stagedep == 1) { //this means there is more than one stage that emigrates, so there is a line for every stage (or two, if sex dependent too)
			if (s != 0) { for (int h = 0; h < 9; h++) { emigfile >> header; } } //to get to the next line for this stage
			if (stages_b[0]->s_emig.sexdep == 1) { //if it's sex dependent
				emigfile >> stages_b[s]->s_emig.m_emig_prob >> stages_b[s]->s_emig.m_D0 >> stages_b[s]->s_emig.m_alpha >> stages_b[s]->s_emig.m_beta;
				for (int h = 0; h < 9; h++) { emigfile >> header; } //then skip the headers until the female information
			}
			emigfile >> stages_b[s]->s_emig.f_emig_prob >> stages_b[s]->s_emig.f_D0 >> stages_b[s]->s_emig.f_alpha >> stages_b[s]->s_emig.f_beta;
			if (stages_b[0]->s_emig.sexdep == 0 && mod->repro > 0) { //if not sex dependent but model is sexual, then male values=female values
				stages_b[s]->s_emig.m_D0 = stages_b[s]->s_emig.f_D0; stages_b[s]->s_emig.m_alpha = stages_b[s]->s_emig.f_alpha;
				stages_b[s]->s_emig.m_beta = stages_b[s]->s_emig.f_beta;
			}
		}
		else { //there is only one stages that emigrates, so there is only one line (or two if sex dependent)
			if (s != stages_b[0]->s_emig.emigstage) { //if this isn't the stage that emigrates, give it -9 for values
				stages_b[s]->s_emig.m_emig_prob = stages_b[s]->s_emig.f_emig_prob = -9;
				stages_b[s]->s_emig.m_D0 = stages_b[s]->s_emig.f_D0 = -9; 
				stages_b[s]->s_emig.m_alpha = stages_b[s]->s_emig.f_alpha = -9;
				stages_b[s]->s_emig.m_beta = stages_b[s]->s_emig.f_beta = -9;
			}
			else {
				if (stages_b[0]->s_emig.sexdep == 1) { //if it's sex dependent
					emigfile >> stages_b[s]->s_emig.m_emig_prob >> stages_b[s]->s_emig.m_D0 >> stages_b[s]->s_emig.m_alpha >> stages_b[s]->s_emig.m_beta;
					for (int h = 0; h < 9; h++) { emigfile >> header; } //then skip the headers until the female information
				}
				emigfile >> stages_b[s]->s_emig.f_emig_prob >> stages_b[s]->s_emig.f_D0 >> stages_b[s]->s_emig.f_alpha >> stages_b[s]->s_emig.f_beta;
				//cout << "emigration prob of stage " << s << "is " << stages_b[s]->s_emig.f_emig_prob << endl;
				if (stages_b[0]->s_emig.sexdep == 0 && mod->repro > 0) { //if not sex dependent but model is sexual, then male values=female values
					stages_b[s]->s_emig.m_emig_prob = stages_b[s]->s_emig.f_emig_prob; stages_b[s]->s_emig.m_D0 = stages_b[s]->s_emig.f_D0; stages_b[s]->s_emig.m_alpha = stages_b[s]->s_emig.f_alpha;
					stages_b[s]->s_emig.m_beta = stages_b[s]->s_emig.f_beta; 
				}
			}
		}
		//cout << " emig prob is " << stages_b[s]->s_emig.f_emig_prob << endl;
	}
	//cout << "emig prob of stage 0" << stages_b[0]->s_emig.m_emig_prob << endl;
	//cout << "emig prob of stage 1" << stages_b[1]->s_emig.m_emig_prob << endl;
	emigfile.close();
	//cout << "population stageinfo vector size is " << stages_b.size() << ", and stage 0 emig prob is " << stages_b[0]->s_emig.f_emig_prob << endl;
}

//void Population::set_emig_info() { //this is only written for one sim at a time right now!
//	int sim, holder;
//	string header;
//	ifstream emigfile(mod->wrkdir + "/Inputs/" + mod->emigfile);
//	int nheaders = 18;
//	for (int h = 0; h < nheaders; h++) {
//		emigfile >> header;
//	}
//	emigfile >> sim >> pop_emig.densdep >> pop_emig.stagedep >> pop_emig.sexdep >> pop_emig.indvar >> pop_emig.iiv_sd >> pop_emig.emigstage >> holder >> holder; //holder holder is because ill step through the stages and sexes in the for loops
//	if (pop_emig.stagedep == 1) { //if emigration is stage dependent
//		for (int s = 0; s < mod->stages; s++) { //for each stage
//			if (pop_emig.sexdep == 1) { //if it's sex dependent, fill the male information 
//				emigfile >> stages[s]->m_ep >> stages[s]->ep_evo >> stages[s]->ep_mut_size >> stages[s]->ep_mut >> stages[s]->m_D0 >> stages[s]->m_alpha >> stages[s]->m_beta >> stages[s]->alpha_mut_size >> stages[s]->beta_mut_size;
//				for (int h = 0; h < 9; h++) { emigfile >> header; } //then skip the headers until the female information
//			}
//			emigfile >> stages[s]->f_ep >> stages[s]->ep_evo >> stages[s]->ep_mut_size >> stages[s]->ep_mut >> stages[s]->f_D0 >> stages[s]->f_alpha >> stages[s]->f_beta >> stages[s]->alpha_mut_size >> stages[s]->beta_mut_size; //fill in female information
//			for (int h = 0; h < 9; h++) { emigfile >> header; } //skip necessary spaces until next stage information
//		}
//	}
//	else { //stage independent
//		if (pop_emig.sexdep == 1) { //if it's sex dependent, get male information
//			emigfile >> pop_emig.m_emig_prob >> pop_emig.ep_evo >> pop_emig.ep_mut_size >> pop_emig.ep_mut >> pop_emig.m_D0 >> pop_emig.m_alpha >> pop_emig.m_beta >> pop_emig.alpha_mut_size >> pop_emig.beta_mut_size;
//			for (int h = 0; h < 9; h++) { emigfile >> header; } //skip necessary spaces until next stage information
//		}
//		emigfile >> pop_emig.f_emig_prob >> pop_emig.ep_evo >> pop_emig.ep_mut_size >> pop_emig.ep_mut >> pop_emig.f_D0 >> pop_emig.f_alpha >> pop_emig.f_beta >> pop_emig.alpha_mut_size >> pop_emig.beta_mut_size; //get female information
//	}
//
//	emigfile.close();
//
//}



void Population::set_transfer_info() {
	int sim, holder;
	string header;
	ifstream transfile(mod->wrkdir + "/Inputs/" + mod->transfile);
	int nheaders = 28;
	for (int h = 0; h < nheaders; h++) {
		transfile >> header;
	}

	
	transfile >> sim >> stages_b[0]->s_trans.stagedep >> stages_b[0]->s_trans.sexdep >> stages_b[0]->s_trans.distmort >> holder >> holder;
	//i want all stages to hold info like stage_dep, sex_dep, etc
	for (int st = 0; st < mod->stages; st++) { //if not stage structured, then this will only happen once
		if (st != 0) {
			stages_b[st]->s_trans.stagedep = stages_b[0]->s_trans.stagedep; stages_b[st]->s_trans.sexdep = stages_b[0]->s_trans.sexdep;
			stages_b[st]->s_trans.distmort = stages_b[0]->s_trans.distmort;
		}
		if (stages_b[0]->s_trans.stagedep == 1) { //this means there is more than one stage that emigrates, so there is a line for every stage (or two, if sex dependent too)
			if (stages_b[0]->s_trans.sexdep == 1) { //if it's sex dependent
				transfile >> stages_b[st]->s_trans.m_rho >> stages_b[st]->s_trans.m_rho_sd >> stages_b[st]->s_trans.m_SL >> stages_b[st]->s_trans.m_comp_SL >> stages_b[st]->s_trans.m_SL_sd;
				for (int ho = 0; ho < 14; ho++) { //because these things will be assigned in the female version, i don't need to write it twice
					transfile >> holder;
				}
				transfile >> stages_b[st]->s_trans.m_mortprob >> stages_b[st]->s_trans.m_slope >> stages_b[st]->s_trans.m_inflpoint;
				stages_b[st]->s_trans.m_mortprob /= 24; //mortality is per day, get it to per hour
				for (int h = 0; h < 6; h++) { transfile >> holder; } //get it to the next line for females

			}
			transfile >> stages_b[st]->s_trans.f_rho >> stages_b[st]->s_trans.f_rho_sd >> stages_b[st]->s_trans.f_SL >> stages_b[st]->s_trans.f_comp_SL >> stages_b[st]->s_trans.f_SL_sd >>
				stages_b[st]->s_trans.f_min_active_time >> stages_b[st]->s_trans.f_min_active_size >> stages_b[st]->s_trans.f_active_sd >> 
				stages_b[st]->s_trans.buoy_min >> stages_b[st]->s_trans.buoy_max >>
				stages_b[st]->s_trans.diel_vert >> stages_b[st]->s_trans.dv_range >> stages_b[st]->s_trans.dv_active >> stages_b[st]->s_trans.min_dv_time >> stages_b[st]->s_trans.min_dv_size >>
				stages_b[st]->s_trans.min_dv_sd >> stages_b[st]->s_trans.down_bias >>
				stages_b[st]->s_trans.dp >> stages_b[st]->s_trans.memory >> stages_b[st]->s_trans.f_mortprob >> stages_b[st]->s_trans.f_slope >> stages_b[st]->s_trans.f_inflpoint;
			stages_b[st]->s_trans.f_mortprob /= 24; //mortality is per day, make it per hour
			if (stages_b[st]->s_trans.f_min_active_time != -9) { stages_b[st]->s_trans.f_min_active_time *= 24; } //this is in days, make it in hours
			if (stages_b[st]->s_trans.min_dv_time != -9) { stages_b[st]->s_trans.min_dv_time *= 24; } //this is in days make it in hours
			if (stages_b[0]->s_trans.sexdep == 0 && mod->repro > 0) { //if it's sex independent but model is still sexual, then male values==female values
				stages_b[st]->s_trans.m_rho = stages_b[st]->s_trans.f_rho; stages_b[st]->s_trans.m_rho_sd = stages_b[st]->s_trans.f_rho_sd; stages_b[st]->s_trans.m_SL = stages_b[st]->s_trans.f_SL;
				stages_b[st]->s_trans.m_comp_SL = stages_b[st]->s_trans.f_comp_SL; stages_b[st]->s_trans.m_SL_sd = stages_b[st]->s_trans.f_SL_sd;
				stages_b[st]->s_trans.m_min_active_time = stages_b[st]->s_trans.f_min_active_time; stages_b[st]->s_trans.m_min_active_size = stages_b[st]->s_trans.f_min_active_size;
				stages_b[st]->s_trans.m_active_sd = stages_b[st]->s_trans.f_active_sd;
				stages_b[st]->s_trans.m_mortprob = stages_b[st]->s_trans.f_mortprob; stages_b[st]->s_trans.m_slope = stages_b[st]->s_trans.f_slope; stages_b[st]->s_trans.m_inflpoint = stages_b[st]->s_trans.f_inflpoint;
				
			}
			for (int h = 0; h < 6; h++) {
				transfile >> header;
			} //go to next stage
		}
		else {
			if (st != stages_b[0]->s_emig.emigstage) { //if this is not the stage that emigrates,
				stages_b[st]->s_trans.m_rho = stages_b[st]->s_trans.f_rho = -9; stages_b[st]->s_trans.m_rho_sd = stages_b[st]->s_trans.f_rho_sd = -9; stages_b[st]->s_trans.m_SL = stages_b[st]->s_trans.f_SL = -9;
				stages_b[st]->s_trans.m_comp_SL = stages_b[st]->s_trans.f_comp_SL = -9; stages_b[st]->s_trans.m_SL_sd = stages_b[st]->s_trans.f_SL_sd = -9;
				stages_b[st]->s_trans.m_mortprob = stages_b[st]->s_trans.f_mortprob = -9; stages_b[st]->s_trans.m_slope = stages_b[st]->s_trans.f_slope = -9; stages_b[st]->s_trans.m_inflpoint = stages_b[st]->s_trans.f_inflpoint = -9;
				stages_b[st]->s_trans.buoy_min = -9; stages_b[st]->s_trans.buoy_max = -9; stages_b[st]->s_trans.diel_vert = -9; stages_b[st]->s_trans.dv_range = -9; stages_b[st]->s_trans.dv_active = -9;
				stages_b[st]->s_trans.min_dv_time = -9; stages_b[st]->s_trans.min_dv_size = -9; stages_b[st]->s_trans.min_dv_sd = -9; stages_b[st]->s_trans.down_bias = -9; stages_b[st]->s_trans.dp = -9; stages_b[st]->s_trans.memory = -9;
			}
			else {
				if (stages_b[0]->s_trans.sexdep == 1) { //if it's sex dependent
					transfile >> stages_b[st]->s_trans.m_rho >> stages_b[st]->s_trans.m_rho_sd >> stages_b[st]->s_trans.m_SL >> stages_b[st]->s_trans.m_comp_SL >> stages_b[st]->s_trans.m_SL_sd;
					for (int ho = 0; ho < 14; ho++) { //because these things will be assigned in the female version, i don't need to write it twice
						transfile >> holder;
					}
					transfile >> stages_b[st]->s_trans.m_mortprob >> stages_b[st]->s_trans.m_slope >> stages_b[st]->s_trans.m_inflpoint;
					stages_b[st]->s_trans.m_mortprob /= 24; //mortality is 
					for (int h = 0; h < 6; h++) { transfile >> holder; } //get it to the next line for females

				}
				transfile >> stages_b[st]->s_trans.f_rho >> stages_b[st]->s_trans.f_rho_sd >> stages_b[st]->s_trans.f_SL >> stages_b[st]->s_trans.f_comp_SL >> stages_b[st]->s_trans.f_SL_sd >>
					stages_b[st]->s_trans.f_min_active_time >> stages_b[st]->s_trans.f_min_active_size >> stages_b[st]->s_trans.f_active_sd >>
					stages_b[st]->s_trans.buoy_min >> stages_b[st]->s_trans.buoy_max >>
					stages_b[st]->s_trans.diel_vert >> stages_b[st]->s_trans.dv_range >> stages_b[st]->s_trans.dv_active >> stages_b[st]->s_trans.min_dv_time >> stages_b[st]->s_trans.min_dv_size >> 
					stages_b[st]->s_trans.min_dv_sd >> stages_b[st]->s_trans.down_bias >>
					stages_b[st]->s_trans.dp >> stages_b[st]->s_trans.memory >> stages_b[st]->s_trans.f_mortprob >> stages_b[st]->s_trans.f_slope >> stages_b[st]->s_trans.f_inflpoint;
				stages_b[st]->s_trans.f_mortprob /= 24; //mortality is
				if (stages_b[st]->s_trans.f_min_active_time != -9) { stages_b[st]->s_trans.f_min_active_time *= 24; }
				if (stages_b[st]->s_trans.min_dv_time != -9) { stages_b[st]->s_trans.min_dv_time *= 24; }
				
				if (stages_b[0]->s_trans.sexdep == 0 && mod->repro > 0) { //if it's sex independent but model is still sexual, then male values==female values
					stages_b[st]->s_trans.m_rho = stages_b[st]->s_trans.f_rho; stages_b[st]->s_trans.m_rho_sd = stages_b[st]->s_trans.f_rho_sd; stages_b[st]->s_trans.m_SL = stages_b[st]->s_trans.f_SL;
					stages_b[st]->s_trans.m_comp_SL = stages_b[st]->s_trans.f_comp_SL; stages_b[st]->s_trans.m_SL_sd = stages_b[st]->s_trans.f_SL_sd;
					stages_b[st]->s_trans.m_min_active_time = stages_b[st]->s_trans.f_min_active_time; stages_b[st]->s_trans.m_min_active_size = stages_b[st]->s_trans.f_min_active_size;
					stages_b[st]->s_trans.m_active_sd = stages_b[st]->s_trans.f_active_sd;
					stages_b[st]->s_trans.m_mortprob = stages_b[st]->s_trans.f_mortprob; stages_b[st]->s_trans.m_slope = stages_b[st]->s_trans.f_slope; stages_b[st]->s_trans.m_inflpoint = stages_b[st]->s_trans.f_inflpoint;
					
				}
			}
		}
		
	}
	cout << "checking dvm time: " << stages_b[0]->s_trans.min_dv_time <<", and min active " << stages_b[0]->s_trans.f_min_active_time << ", and mort prob: " << stages_b[0]->s_trans.f_mortprob << endl;
	transfile.close();
	if (mod->size_or_time==0) { //size dependent, therefore need growth information
		ifstream growthfile(mod->wrkdir + "/Inputs/" + mod->growthfile);
		int nheaders = 12; //this is the maximum number of columns, but some columns are omitted, depending on method used
		for (int h = 0; h < nheaders; h++) {
			growthfile >> header;
		}

		growthfile >> sim >> stages_b[0]->s_trans.pop_grow.stage_dep >> stages_b[0]->s_trans.pop_grow.sex_dep >> stages_b[0]->s_trans.pop_grow.method;
		for (int st = 0; st < mod->stages; st++) {
			if (st != 0) {
				stages_b[st]->s_trans.pop_grow.sex_dep = stages_b[0]->s_trans.pop_grow.sex_dep; stages_b[st]->s_trans.pop_grow.method = stages_b[0]->s_trans.pop_grow.method; 
				stages_b[st]->s_trans.pop_grow.stage_dep = stages_b[0]->s_trans.pop_grow.stage_dep;
			}
			if (stages_b[0]->s_trans.pop_grow.stage_dep == 1) { //means there will be a line for each stage (or two, if sex dependent)
				if (stages_b[0]->s_trans.pop_grow.sex_dep == 1) { //sex dependent
					growthfile >> stages_b[st]->s_trans.pop_grow.m_m >> stages_b[st]->s_trans.pop_grow.m_m_sd >> stages_b[st]->s_trans.pop_grow.m_b
						>> stages_b[st]->s_trans.pop_grow.m_Linf >> stages_b[st]->s_trans.pop_grow.m_Linf_sd >> stages_b[st]->s_trans.pop_grow.Ti >> stages_b[st]->s_trans.pop_grow.Ti_sd >> stages_b[st]->s_trans.pop_grow.m_G_K;
						
					for (int h = 0; h < 4; h++) { growthfile >> holder; } //get it to the next line for females
					//then females
				}
				growthfile >> stages_b[st]->s_trans.pop_grow.f_m >> stages_b[st]->s_trans.pop_grow.f_m_sd >> stages_b[st]->s_trans.pop_grow.f_b
					>> stages_b[st]->s_trans.pop_grow.f_Linf >> stages_b[st]->s_trans.pop_grow.f_Linf_sd >> stages_b[st]->s_trans.pop_grow.Ti >> stages_b[st]->s_trans.pop_grow.Ti_sd >> stages_b[st]->s_trans.pop_grow.f_G_K;
					
				if (stages_b[0]->s_trans.pop_grow.sex_dep == 0 && mod->repro > 0) { //if not sex dependent but model is still sexual, male values=female values
					stages_b[st]->s_trans.pop_grow.m_m = stages_b[st]->s_trans.pop_grow.f_m; stages_b[st]->s_trans.pop_grow.m_m_sd = stages_b[st]->s_trans.pop_grow.f_m_sd;
					stages_b[st]->s_trans.pop_grow.m_b = stages_b[st]->s_trans.pop_grow.f_b;
					stages_b[st]->s_trans.pop_grow.m_Linf = stages_b[st]->s_trans.pop_grow.f_Linf; stages_b[st]->s_trans.pop_grow.m_Linf_sd = stages_b[st]->s_trans.pop_grow.f_Linf_sd;
					stages_b[st]->s_trans.pop_grow.m_G_K = stages_b[st]->s_trans.pop_grow.f_G_K;
				}
				for (int h = 0; h < 4; h++) { growthfile >> holder; } //go to next stage
			}
			else {//there is only one stage dispersing
				if (st != stages_b[0]->s_emig.emigstage) { //if this is not the stage that emigrates,
					stages_b[st]->s_trans.pop_grow.m_m = stages_b[st]->s_trans.pop_grow.f_m=-9; stages_b[st]->s_trans.pop_grow.m_m_sd = stages_b[st]->s_trans.pop_grow.f_m_sd=-9;
					stages_b[st]->s_trans.pop_grow.m_b = stages_b[st]->s_trans.pop_grow.f_b=-9;
					stages_b[st]->s_trans.pop_grow.m_Linf = stages_b[st]->s_trans.pop_grow.f_Linf=-9; stages_b[st]->s_trans.pop_grow.m_Linf_sd = stages_b[st]->s_trans.pop_grow.f_Linf_sd=-9;
					stages_b[st]->s_trans.pop_grow.m_G_K = stages_b[st]->s_trans.pop_grow.f_G_K=-9;
				}
				else {
					if (stages_b[0]->s_trans.pop_grow.sex_dep == 1) { //sex dependent
						growthfile >> stages_b[st]->s_trans.pop_grow.m_m >> stages_b[st]->s_trans.pop_grow.m_m_sd >> stages_b[st]->s_trans.pop_grow.m_b
							>> stages_b[st]->s_trans.pop_grow.m_Linf >> stages_b[st]->s_trans.pop_grow.m_Linf_sd >> stages_b[st]->s_trans.pop_grow.Ti >> stages_b[st]->s_trans.pop_grow.Ti_sd >> stages_b[st]->s_trans.pop_grow.m_G_K;

						for (int h = 0; h < 5; h++) { growthfile >> holder; } //get it to the next line for females
						//then females
					}
					growthfile >> stages_b[st]->s_trans.pop_grow.f_m >> stages_b[st]->s_trans.pop_grow.f_m_sd >> stages_b[st]->s_trans.pop_grow.f_b
						>> stages_b[st]->s_trans.pop_grow.f_Linf >> stages_b[st]->s_trans.pop_grow.f_Linf_sd >> stages_b[st]->s_trans.pop_grow.Ti >> stages_b[st]->s_trans.pop_grow.Ti_sd >> stages_b[st]->s_trans.pop_grow.f_G_K;

					if (stages_b[0]->s_trans.pop_grow.sex_dep == 0 && mod->repro > 0) { //if not sex dependent but model is still sexual, male values=female values
						stages_b[st]->s_trans.pop_grow.m_m = stages_b[st]->s_trans.pop_grow.f_m; stages_b[st]->s_trans.pop_grow.m_m_sd = stages_b[st]->s_trans.pop_grow.f_m_sd;
						stages_b[st]->s_trans.pop_grow.m_b = stages_b[st]->s_trans.pop_grow.f_b;
						stages_b[st]->s_trans.pop_grow.m_Linf = stages_b[st]->s_trans.pop_grow.f_Linf; stages_b[st]->s_trans.pop_grow.m_Linf_sd = stages_b[st]->s_trans.pop_grow.f_Linf_sd;
						stages_b[st]->s_trans.pop_grow.m_G_K = stages_b[st]->s_trans.pop_grow.f_G_K;
					}
				}
				
			}
			
		}
		//cout << "Ti is " << stages_b[0]->s_trans.pop_grow.Ti << endl;
		//cout << "min active size is " << stages_b[0]->s_trans.m_min_active_size << endl;
		//cout << "active sd is " << stages_b[0]->s_trans.m_active_sd << endl;
		
		growthfile.close();
	}
	
}

//void Population::set_transfer_info() {
//	int sim, holder;
//	string header;
//	ifstream transfile(mod->wrkdir + "/Inputs/" + mod->transfile);
//	int nheaders = 30;
//	for (int h = 0; h < nheaders; h++) {
//		transfile >> header;
//	}
//
//	transfile >> sim >> pop_trans.stagedep >> pop_trans.sexdep >> pop_trans.distmort >> holder >> holder; //add in the holders because i dont need to read the stage and sex values. if either are important i go through all stages and all sexes
//	if (pop_trans.stagedep == 1) { //stage dependent transfer parameters
//		for (int s = 0; s < mod->stages; s++) {
//			if (pop_trans.sexdep == 1) { //sex dependent
//				transfile >> stages[s]->m_rho >> stages[s]->m_rho_sd >> stages[s]->m_SL >> stages[s]->m_comp_SL >> stages[s]->m_SL_sd >>
//					stages[s]->m_min_active_time >> stages[s]->m_min_active_size >> stages[s]->m_active_sd >> stages[s]->active_evo >> stages[s]->active_mut >> stages[s]->active_mut_size >>
//					stages[s]->buoy_min >> stages[s]->buoy_max >>
//					stages[s]->diel_vert >> stages[s]->dv_range >> stages[s]->dv_active >> stages[s]->min_dv_time >> stages[s]->min_dv_size >> stages[s]->down_bias >>
//					stages[s]->dp >> stages[s]->memory >> stages[s]->m_mortprob >> stages[s]->m_slope >> stages[s]->m_inflpoint;
//				stages[s]->m_mortprob /= 24; //mortality is 
//				for (int h = 0; h < 6; h++) { transfile >> holder; } //get it to the next line for females
//			}
//			transfile >> stages[s]->f_rho >> stages[s]->f_rho_sd >> stages[s]->f_SL >> stages[s]->f_comp_SL >> stages[s]->f_SL_sd >>
//				stages[s]->f_min_active_time >> stages[s]->f_min_active_size >> stages[s]->f_active_sd >> stages[s]->active_evo >> stages[s]->active_mut >> stages[s]->active_mut_size >>
//				stages[s]->buoy_min >> stages[s]->buoy_max >>
//				stages[s]->diel_vert >> stages[s]->dv_range >> stages[s]->dv_active >> stages[s]->min_dv_time >> stages[s]->min_dv_size >> stages[s]->down_bias >>
//				stages[s]->dp >> stages[s]->memory >> stages[s]->f_mortprob >> stages[s]->f_slope >> stages[s]->f_inflpoint;
//			for (int h = 0; h < 6; h++) { transfile >> holder; } //get it to the next stage
//		}
//	}
//
//	else { //stage independent
//		if (pop_trans.sexdep == 1) { //sex dependent
//			transfile >> pop_trans.m_rho >> pop_trans.m_rho_sd >> pop_trans.m_SL >> pop_trans.m_comp_SL >> pop_trans.m_SL_sd >>
//				pop_trans.m_min_active_time >> pop_trans.m_min_active_size >> pop_trans.m_active_sd >> pop_trans.active_evo >> pop_trans.active_mut >> pop_trans.active_mut_size >>
//				pop_trans.buoy_min >> pop_trans.buoy_max >>
//				pop_trans.diel_vert >> pop_trans.dv_range >> pop_trans.dv_active >> pop_trans.min_dv_time >> pop_trans.min_dv_size >> pop_trans.down_bias >>
//				pop_trans.dp >> pop_trans.memory >> pop_trans.m_mortprob >> pop_trans.m_slope >> pop_trans.m_inflpoint;
//			for (int h = 0; h < 6; h++) { transfile >> holder; } //get it to the next line for females
//		}
//		transfile >> pop_trans.f_rho >> pop_trans.f_rho_sd >> pop_trans.f_SL >> pop_trans.f_comp_SL >> pop_trans.f_SL_sd >>
//			pop_trans.f_min_active_time >> pop_trans.f_min_active_size >> pop_trans.f_active_sd >> pop_trans.active_evo >> pop_trans.active_mut >> pop_trans.active_mut_size >>
//			pop_trans.buoy_min >> pop_trans.buoy_max >>
//			pop_trans.diel_vert >> pop_trans.dv_range >> pop_trans.dv_active >> pop_trans.min_dv_time >> pop_trans.min_dv_size >> pop_trans.down_bias >>
//			pop_trans.dp >> pop_trans.memory >> pop_trans.f_mortprob >> pop_trans.f_slope >> pop_trans.f_inflpoint;
//	}
//	transfile.close();
//	if (mod->size_or_time == 0) { //size dependent, therefore need growth information
//		ifstream growthfile(mod->wrkdir + "/Inputs/" + mod->growthfile);
//		int nheaders = 15; //this is the maximum number of columns, but some columns are omitted, depending on method used
//		for (int h = 0; h < nheaders; h++) {
//			growthfile >> header;
//		}
//
//		growthfile >> sim >> pop_trans.pop_grow.sex_dep >> pop_trans.pop_grow.method >> pop_trans.pop_grow.growth_evo >> pop_trans.pop_grow.growth_mut;
//		//if (pop_trans.stagedep == 1) { //even though growth parameters only apply to one stage, need to do this because it's saved in a different location
//		//	for (int s = 0; s < mod->stages; s++) {
//		//		if (pop_trans.pop_grow.sex_dep == 1) { //if it's sex dependent, do males first
//		//			growthfile >> stages[s]->pop_grow.m_m >> stages[s]->pop_grow.m_m_sd >> stages[s]->pop_grow.m_b >> stages[s]->pop_grow.m_mut_size >> stages[s]->pop_grow.b_mut_size
//		//				>> stages[s]->pop_grow.m_Linf >> stages[s]->pop_grow.m_Linf_sd >> stages[s]->pop_grow.m_G_K >> stages[s]->pop_grow.Linf_mut_size >> stages[s]->pop_grow.GK_mut_size;
//		//			for (int h = 0; h < 5; h++) { growthfile >> holder; } //get it to the next line for females
//		//		}
//		//		//then females
//		//		growthfile >> stages[s]->pop_grow.f_m >> stages[s]->pop_grow.f_m_sd >> stages[s]->pop_grow.f_b >> stages[s]->pop_grow.m_mut_size >> stages[s]->pop_grow.b_mut_size
//		//			>> stages[s]->pop_grow.f_Linf >> stages[s]->pop_grow.f_Linf_sd >> stages[s]->pop_grow.f_G_K >> stages[s]->pop_grow.Linf_mut_size >> stages[s]->pop_grow.GK_mut_size;
//		//	}
//		//}
//		//else {
//		if (pop_trans.pop_grow.sex_dep == 1) { //if it's sex dependent, do males first
//			growthfile >> pop_trans.pop_grow.m_m >> pop_trans.pop_grow.m_m_sd >> pop_trans.pop_grow.m_b >> pop_trans.pop_grow.m_mut_size >> pop_trans.pop_grow.b_mut_size
//				>> pop_trans.pop_grow.m_Linf >> pop_trans.pop_grow.m_Linf_sd >> pop_trans.pop_grow.m_G_K >> pop_trans.pop_grow.Linf_mut_size >> pop_trans.pop_grow.GK_mut_size;
//			for (int h = 0; h < 5; h++) { growthfile >> holder; } //get it to the next line for females
//		}
//		//then females
//		growthfile >> pop_trans.pop_grow.f_m >> pop_trans.pop_grow.f_m_sd >> pop_trans.pop_grow.f_b >> pop_trans.pop_grow.m_mut_size >> pop_trans.pop_grow.b_mut_size
//			>> pop_trans.pop_grow.f_Linf >> pop_trans.pop_grow.f_Linf_sd >> pop_trans.pop_grow.f_G_K >> pop_trans.pop_grow.Linf_mut_size >> pop_trans.pop_grow.GK_mut_size;
//		//}
//		growthfile.close();
//	}
//
//}

void Population::set_settle_info() {
	int sim, holder;
	string header;
	ifstream settfile(mod->wrkdir + "/Inputs/" + mod->setfile);
	int nheaders = 22;
	for (int h = 0; h < nheaders; h++) {
		settfile >> header;
	}
	settfile >> sim >> stages_b[0]->s_sett.stagedep >> stages_b[0]->s_sett.sexdep >> stages_b[0]->s_sett.densdep >> holder >> holder;
	for (int st = 0; st < mod->stages; st++) {
		if (st != 0) {
			stages_b[st]->s_sett.stagedep = stages_b[0]->s_sett.stagedep; stages_b[st]->s_sett.sexdep = stages_b[0]->s_sett.sexdep;
			stages_b[st]->s_sett.densdep = stages_b[0]->s_sett.stagedep;
		}
		if (stages_b[0]->s_sett.stagedep == 1) {
			if (stages_b[0]->s_sett.sexdep == 1) { //sex dependent
				settfile >> stages_b[st]->s_sett.m_findmate >> stages_b[st]->s_sett.m_pld >> stages_b[st]->s_sett.m_pld_sd >> stages_b[st]->s_sett.m_comp >> stages_b[st]->s_sett.m_comp_size >> stages_b[st]->s_sett.m_comp_sd >>
					stages_b[st]->s_sett.m_S0 >> stages_b[st]->s_sett.m_alphaS >> stages_b[st]->s_sett.m_betaS >> stages_b[st]->s_sett.m_sett_sd;
				stages_b[st]->s_sett.m_pld *= 24;
				if (stages_b[st]->s_sett.m_pld_sd != -9) { stages_b[st]->s_sett.m_pld_sd *= 24; } //don't overwrite 
				if (stages_b[st]->s_sett.m_comp != -9) { stages_b[st]->s_sett.m_comp *= 24; }//need to multiply by 24 to make it in hours because speed of currents are in hours
				for (int h = 0; h < 11; h++) { settfile >> header; } //get to next sex
			}
			settfile >> stages_b[st]->s_sett.f_findmate >> stages_b[st]->s_sett.f_pld >> stages_b[st]->s_sett.f_pld_sd >> stages_b[st]->s_sett.f_comp >> 
				stages_b[st]->s_sett.f_comp_size >> stages_b[st]->s_sett.f_comp_sd >>
				stages_b[st]->s_sett.f_S0 >> stages_b[st]->s_sett.f_alphaS >> stages_b[st]->s_sett.f_betaS >> stages_b[st]->s_sett.f_sett_sd;
			stages_b[st]->s_sett.f_pld *= 24;
			if (stages_b[st]->s_sett.f_pld_sd != -9) { stages_b[st]->s_sett.f_pld_sd *= 24; } //don't overwrite 
			if (stages_b[st]->s_sett.f_comp != -9) { stages_b[st]->s_sett.f_comp *= 24; }//need to multiply by 24 to make it in hours because speed of currents are in hours
			settfile >> stages_b[st]->s_sett.buffer >> stages_b[st]->s_sett.buffer_xycells >> stages_b[st]->s_sett.buffer_zlayers >> stages_b[st]->s_sett.buffer_cond >> stages_b[st]->s_sett.buffer_thresh >> stages_b[st]->s_sett.buffer_capture;
			if (stages_b[0]->s_sett.sexdep == 0 && mod->repro > 0) { //not sex dependent but sexual model, male values=female values
				stages_b[st]->s_sett.m_findmate = stages_b[st]->s_sett.f_findmate; stages_b[st]->s_sett.m_pld = stages_b[st]->s_sett.f_pld;
				stages_b[st]->s_sett.m_pld_sd = stages_b[st]->s_sett.f_pld_sd; stages_b[st]->s_sett.m_comp = stages_b[st]->s_sett.f_comp;
				stages_b[st]->s_sett.m_comp_size = stages_b[st]->s_sett.f_comp_size; stages_b[st]->s_sett.m_comp_sd = stages_b[st]->s_sett.f_comp_sd;
				stages_b[st]->s_sett.m_S0 = stages_b[st]->s_sett.f_S0; stages_b[st]->s_sett.m_alphaS = stages_b[st]->s_sett.f_alphaS;
				stages_b[st]->s_sett.m_betaS = stages_b[st]->s_sett.f_betaS; stages_b[st]->s_sett.m_sett_sd = stages_b[st]->s_sett.f_sett_sd;
			}
			for (int h = 0; h < 6; h++) { settfile >> header; } //get to next stage
		}
		else {
			if (st != stages_b[0]->s_emig.emigstage) { //if this is not the stage that emigrates,
				stages_b[st]->s_sett.m_findmate = stages_b[st]->s_sett.f_findmate=-9; stages_b[st]->s_sett.m_pld = stages_b[st]->s_sett.f_pld=-9;
				stages_b[st]->s_sett.m_pld_sd = stages_b[st]->s_sett.f_pld_sd=-9; stages_b[st]->s_sett.m_comp = stages_b[st]->s_sett.f_comp=-9;
				stages_b[st]->s_sett.m_comp_size = stages_b[st]->s_sett.f_comp_size=-9; stages_b[st]->s_sett.m_comp_sd = stages_b[st]->s_sett.f_comp_sd=-9;
				stages_b[st]->s_sett.m_S0 = stages_b[st]->s_sett.f_S0=-9; stages_b[st]->s_sett.m_alphaS = stages_b[st]->s_sett.f_alphaS=-9;
				stages_b[st]->s_sett.m_betaS = stages_b[st]->s_sett.f_betaS=-9; stages_b[st]->s_sett.m_sett_sd = stages_b[st]->s_sett.f_sett_sd=-9;
			}
			else {
				if (stages_b[0]->s_sett.sexdep == 1) { //sex dependent
					settfile >> stages_b[st]->s_sett.m_findmate >> stages_b[st]->s_sett.m_pld >> stages_b[st]->s_sett.m_pld_sd >> stages_b[st]->s_sett.m_comp >> stages_b[st]->s_sett.m_comp_size >> stages_b[st]->s_sett.m_comp_sd >>
						stages_b[st]->s_sett.m_S0 >> stages_b[st]->s_sett.m_alphaS >> stages_b[st]->s_sett.m_betaS >> stages_b[st]->s_sett.m_sett_sd;
					stages_b[st]->s_sett.m_pld *= 24;
					if (stages_b[st]->s_sett.m_pld_sd != -9) { stages_b[st]->s_sett.m_pld_sd *= 24; } //don't overwrite 
					if (stages_b[st]->s_sett.m_comp != -9) { stages_b[st]->s_sett.m_comp *= 24; }//need to multiply by 24 to make it in hours because speed of currents are in hours
					for (int h = 0; h < 11; h++) { settfile >> header; } //get to next sex
				}
				settfile >> stages_b[st]->s_sett.f_findmate >> stages_b[st]->s_sett.f_pld >> stages_b[st]->s_sett.f_pld_sd >> stages_b[st]->s_sett.f_comp >> stages_b[st]->s_sett.f_comp_size >> stages_b[st]->s_sett.f_comp_sd >>
					stages_b[st]->s_sett.f_S0 >> stages_b[st]->s_sett.f_alphaS >> stages_b[st]->s_sett.f_betaS >> stages_b[st]->s_sett.f_sett_sd;
				stages_b[st]->s_sett.f_pld *= 24;
				if (stages_b[st]->s_sett.f_pld_sd != -9) { stages_b[st]->s_sett.f_pld_sd *= 24; } //don't overwrite 
				if (stages_b[st]->s_sett.f_comp != -9) { stages_b[st]->s_sett.f_comp *= 24; }//need to multiply by 24 to make it in hours because speed of currents are in hours
				settfile >> stages_b[st]->s_sett.buffer >> stages_b[st]->s_sett.buffer_xycells >> stages_b[st]->s_sett.buffer_zlayers >> stages_b[st]->s_sett.buffer_cond >> stages_b[st]->s_sett.buffer_thresh >> stages_b[st]->s_sett.buffer_capture;
				if (stages_b[0]->s_sett.sexdep == 0 && mod->repro > 0) { //not sex dependent but sexual model, male values=female values
					stages_b[st]->s_sett.m_findmate = stages_b[st]->s_sett.f_findmate; stages_b[st]->s_sett.m_pld = stages_b[st]->s_sett.f_pld;
					stages_b[st]->s_sett.m_pld_sd = stages_b[st]->s_sett.f_pld_sd; stages_b[st]->s_sett.m_comp = stages_b[st]->s_sett.f_comp;
					stages_b[st]->s_sett.m_comp_size = stages_b[st]->s_sett.f_comp_size; stages_b[st]->s_sett.m_comp_sd = stages_b[st]->s_sett.f_comp_sd;
					stages_b[st]->s_sett.m_S0 = stages_b[st]->s_sett.f_S0; stages_b[st]->s_sett.m_alphaS = stages_b[st]->s_sett.f_alphaS;
					stages_b[st]->s_sett.m_betaS = stages_b[st]->s_sett.f_betaS; stages_b[st]->s_sett.m_sett_sd = stages_b[st]->s_sett.f_sett_sd;
				}
			}
		}
		
	}
	
	//cout << "pld is " << stages_b[0]->s_sett.m_pld << ", with sd " << stages_b[0]->s_sett.m_pld_sd << endl;
	//cout << "comp size is " << stages_b[0]->s_sett.m_comp_size << ", with sd " << stages_b[0]->s_sett.m_comp_sd << endl;
	//cout << "comp time is " << stages_b[0]->s_sett.f_comp << ", with sd " << stages_b[0]->s_sett.f_comp_sd << endl;
	//cout << "settlemenet density dependence is " << stages_b[0]->s_sett.densdep << endl;
	//cout << "S0 is " << stages_b[0]->s_sett.f_S0 << ", with alphaS " << stages_b[0]->s_sett.f_alphaS << ", with betaS " << stages_b[0]->s_sett.f_betaS << endl;
	settfile.close();
}

//void Population::set_settle_info() {
//	int sim, holder;
//	string header;
//	ifstream settfile(mod->wrkdir + "/Inputs/" + mod->setfile);
//	int nheaders = 29;
//	for (int h = 0; h < nheaders; h++) {
//		settfile >> header;
//	}
//	settfile >> sim >> pop_sett.stagedep >> pop_sett.sexdep >> pop_sett.densdep >> holder >> holder;
//	
//	if (pop_sett.stagedep == 1) { //stage-dependent settlement parameters
//		for (int s = 0; s < mod->stages; s++) {
//			if (pop_sett.sexdep == 1) { //sex dependent
//				settfile >> stages[s]->m_findmate >> stages[s]->m_pld >> stages[s]->m_pld_sd >> stages[s]->m_comp >> stages[s]->m_comp_size >> stages[s]->m_comp_sd >> stages[s]->comp_evo >> stages[s]->comp_mut >> stages[s]->comp_mut_size >>
//				stages[s]->m_S0 >> stages[s]->m_alphaS >> stages[s]->m_betaS >> stages[s]->S0_evo >> stages[s]->S0_mut >> stages[s]->S0_mut_size >> 
//					stages[s]->alphaS_mut_size >> stages[s]->betaS_mut_size >> stages[s]->m_sett_sd;
//				stages[s]->m_pld *= 24; stages[s]->m_pld_sd *= 24;  stages[s]->m_comp *= 24; //need to multiply by 24 to make it in hours because speed of currents are in hours
//				settfile>> stages[s]->buffer >> stages[s]->buffer_xycells >> stages[s]->buffer_zlayers >> stages[s]->buffer_cond >> stages[s]->buffer_thresh;
//				for (int h = 0; h < 6; h++) { settfile >> header; } //get to next sex
//			}
//			settfile >> stages[s]->f_findmate >> stages[s]->f_pld >> stages[s]->f_pld_sd >> stages[s]->f_comp >> stages[s]->f_comp_size >> stages[s]->f_comp_sd >> stages[s]->comp_evo >> stages[s]->comp_mut >> stages[s]->comp_mut_size >>
//				stages[s]->f_S0 >> stages[s]->f_alphaS >> stages[s]->f_betaS >> stages[s]->S0_evo >> stages[s]->S0_mut >> stages[s]->S0_mut_size >>
//				stages[s]->alphaS_mut_size >> stages[s]->betaS_mut_size >> stages[s]->f_sett_sd;
//			stages[s]->f_pld *= 24; stages[s]->f_pld_sd *= 24; stages[s]->f_comp *= 24; //need to multiply by 24 to make it in hours because speed of currents are in hours
//			settfile >> stages[s]->buffer >> stages[s]->buffer_xycells >> stages[s]->buffer_zlayers >> stages[s]->buffer_cond >> stages[s]->buffer_thresh;
//		}
//		for (int h = 0; h < 6; h++) { settfile >> header; } //get to next stage
//	}
//	else { //stage-independent
//		if (pop_sett.sexdep == 1) { //sex dependent
//			settfile >> pop_sett.m_findmate >> pop_sett.m_pld >> pop_sett.m_comp >> pop_sett.m_comp_size >> pop_sett.m_comp_sd >> pop_sett.comp_evo >> pop_sett.comp_mut >> pop_sett.comp_mut_size >>
//				pop_sett.m_S0 >> pop_sett.m_alphaS >> pop_sett.m_betaS >> pop_sett.S0_evo >> pop_sett.S0_mut >> pop_sett.S0_mut_size >>
//				pop_sett.alphaS_mut_size >> pop_sett.betaS_mut_size >> pop_sett.m_sett_sd;
//			pop_sett.m_pld *= 24; pop_sett.m_comp *= 24;
//			settfile >> pop_sett.buffer >> pop_sett.buffer_xycells >> pop_sett.buffer_zlayers >> pop_sett.buffer_cond >> pop_sett.buffer_thresh;
//			for (int h = 0; h < 6; h++) { settfile >> header; } //get to next sex
//		}
//		settfile >> pop_sett.f_findmate >> pop_sett.f_pld >> pop_sett.f_comp >> pop_sett.f_comp_size >> pop_sett.f_comp_sd >> pop_sett.comp_evo >> pop_sett.comp_mut >> pop_sett.comp_mut_size >>
//			pop_sett.f_S0 >> pop_sett.f_alphaS >> pop_sett.f_betaS >> pop_sett.S0_evo >> pop_sett.S0_mut >> pop_sett.S0_mut_size >>
//			pop_sett.alphaS_mut_size >> pop_sett.betaS_mut_size >> pop_sett.f_sett_sd >>
//			pop_sett.buffer >> pop_sett.buffer_xycells >> pop_sett.buffer_zlayers >> pop_sett.buffer_cond >> pop_sett.buffer_thresh;
//		pop_sett.f_pld *= 24; pop_sett.f_comp *= 24;
//	}
//
//	settfile.close();
//}

void Population::set_evo_info() {
	ifstream evo_file(mod->wrkdir + "/Inputs/" + mod->evofile);
	string header;

	int nheaders = 27; //when i start using simulations, change this to be dependent on sim#
	for (int h = 0; h < nheaders; h++) {
		evo_file >> header;
	}

	int sim; //not doing anything with this yet
	evo_file >> sim >> pop_dyn.pop_evo.ep_evo >> pop_dyn.pop_evo.ep_mut_size >> pop_dyn.pop_evo.ep_mut >> pop_dyn.pop_evo.alpha_mut_size >> pop_dyn.pop_evo.beta_mut_size >>
		pop_dyn.pop_evo.active_evo >> pop_dyn.pop_evo.active_mut >> pop_dyn.pop_evo.active_mut_size >> pop_dyn.pop_evo.dv_evo >> pop_dyn.pop_evo.dv_mut >> pop_dyn.pop_evo.dv_mut_size >>
		pop_dyn.pop_evo.growth_evo >> pop_dyn.pop_evo.growth_mut >> pop_dyn.pop_evo.m_mut_size >> pop_dyn.pop_evo.b_mut_size >> pop_dyn.pop_evo.Linf_mut_size >> pop_dyn.pop_evo.GK_mut_size >> pop_dyn.pop_evo.Ti_mut_size >>
		pop_dyn.pop_evo.comp_evo >> pop_dyn.pop_evo.comp_mut >> pop_dyn.pop_evo.comp_mut_size >>
		pop_dyn.pop_evo.S0_evo >> pop_dyn.pop_evo.S0_mut >> pop_dyn.pop_evo.S0_mut_size >> pop_dyn.pop_evo.alphaS_mut_size >> pop_dyn.pop_evo.betaS_mut_size;
	//cout << "dv evo is " << pop_dyn.pop_evo.dv_evo << ", with mut prob " << pop_dyn.pop_evo.dv_mut << endl;
}

void Population::set_pop() {
	ifstream pop_file(mod->wrkdir + "/Inputs/" + mod->popfile);
	string header;
	
	int nheaders = 20; //when i start using simulations, change this to be dependent on sim#
	for (int h = 0; h < nheaders; h++) {
		pop_file >> header;
	}
	
	int sim; //not doing anything with this yet
	string depthfile;
	pop_file >> sim >> pop_dyn.bc >> pop_dyn.R >> pop_dyn.max_age >> pop_dyn.fishing_mort >> pop_dyn.Nrepro >> pop_dyn.RepInt >> pop_dyn.PRep >> pop_dyn.juv_rel >> 
		pop_dyn.juv_size >> pop_dyn.iiv_juv_size >> pop_dyn.juv_size_sd >> 
		pop_dyn.survsched >> pop_dyn.Hsize >> depthfile >> pop_dyn.fecdensdep >> pop_dyn.devdensdep >> pop_dyn.devdenscoef >> pop_dyn.survdensdep >> pop_dyn.survdenscoef;
	cout << "juv release is " << pop_dyn.juv_rel << endl;
	//cout << "finish mort is " << pop_dyn.fishing_mort << endl;
	//pop_dyn.pld *= 24; pop_dyn.comp *= 24; //it's given in days, but turn it into hours so that it is the same unit as current angle/speed measurements

	//cout << "max age is " << pop_dyn.max_age << endl;
	pop_file.close();

	//read in depth preference data
	ifstream depth_file(mod->wrkdir + "/Inputs/" + depthfile);
	nheaders = 9;
	for (int h = 0; h < nheaders; h++) {
		depth_file >> header;
	}
	depth_file >> sim >> pop_dyn.DepthPref >> pop_dyn.Depthmin >> pop_dyn.Depthmax >> pop_dyn.DepthSex;

	if (pop_dyn.DepthSex == 1) { //if depth-preference is sex dependent 
		for (int s = 0; s < 2; s++) {
			depth_info* s_depth = new depth_info;
			depth_file >> s_depth->Depthmin >> s_depth->Depthmax;
			sex_depths.push_back(s_depth);
		}

	}
	depth_file.close();
	if (mod->stagestruct == 1) { set_stages(); } //if its a stage structured population, read in stage information (transition matrix)
	else { //not stagestructured, so there is only info for one stage and i don't need to fill in the transition matrix information
		stage_infob* stage = new stage_infob;
		stage->stage = 0;
		stages_b.push_back(stage);
	}
	set_emig_info(); //read emigration file and save info
	set_transfer_info();//read transfer file and save info
	set_settle_info(); //read settlement file and save info
	set_evo_info(); //set evolution information
	//cout << "linf_mut_size is " << pop_dyn.pop_evo.Linf_mut_size << " and Ti_mut_size is " << pop_dyn.pop_evo.Ti_mut_size << endl;
}

void Population::create_subpop(Cell * pcell ) {
	Subpopulation* sp = new Subpopulation(pcell, mod, lands, &pop_dyn, &pop_matrix, stages_b); //create a new subpopulation for each patch
	pcell->c_subpop = sp; //assign the pointer to the subpop to that cell
	subpops.push_back(sp); //add the subpop to the Population's vector of subpops
}

void Population::create_subpop(Patch* ppatch) {
	Subpopulation* sp = new Subpopulation(ppatch, mod, lands, &pop_dyn, &pop_matrix, stages_b); //create a new subpopulation for each patch
	ppatch->p_subpop = sp; //assign the pointer to the subpop to that patch
	subpops.push_back(sp); //add the subpop to the Population's vector of subpops
}
//void Population::create_subpop(Cell * pcell) {
//	Subpopulation* sp = new Subpopulation(pcell, mod, lands, &pop_dyn, &pop_matrix, &pop_emig, &pop_trans, &pop_sett, stages); //create a new subpopulation for each patch
//	pcell->c_subpop = sp; //assign the pointer to the subpop to that cell
//	subpops.push_back(sp); //add the subpop to the Population's vector of subpops
//}
//
//void Population::create_subpop(Patch* ppatch) {
//	Subpopulation* sp = new Subpopulation(ppatch, mod, lands, &pop_dyn, &pop_matrix, &pop_emig, &pop_trans, &pop_sett, stages); //create a new subpopulation for each patch
//	ppatch->p_subpop = sp; //assign the pointer to the subpop to that patch
//	subpops.push_back(sp); //add the subpop to the Population's vector of subpops
//}

///////getting functions /////////////
float Population::get_attribute(int attribute) {
	//attributes: 

	//1: bc
	//2: R
	//3: pld
	//4: comp
	//5: max_age
	//6: Nrepro
	//7: afr
	//8: RepInt
	//9: PRep
	//10:depthmin
	//11: depthmax
	//12: EP
	//13: D0
	//14: alpha
	//15: beta

	float value;
	switch (attribute) {

	case 1:
		value = pop_dyn.bc;
		break;
	case 2:
		value = pop_dyn.R;
		break;
	/*case 3:
		value = pop_dyn.pld;
		break;
	case 4:
		value = pop_dyn.comp;
		break;*/
	case 5:
		value = pop_dyn.max_age;
		break;
	case 6:
		value = pop_dyn.Nrepro;
		break;
	/*case 7:
		value = pop_dyn.afr;
		break;*/
	case 8:
		value = pop_dyn.RepInt;
		break;
	case 9:
		value = pop_dyn.PRep;
		break;
	case 10:
		value = pop_dyn.Depthmin;
		cout << "called 10: " << pop_dyn.Depthmin;
		break;
	case 11:
		value = pop_dyn.Depthmax;
		cout << "called 11: " << pop_dyn.Depthmax << endl;
		break;
	case 12:
		//value = pop_emig.f_emig_prob;
		value = stages_b[0]->s_emig.f_emig_prob;
		break;
	case 13:
		//value = pop_emig.f_D0;
		value = stages_b[0]->s_emig.f_D0;
		break;
	case 14:
		//value = pop_emig.f_alpha;
		value = stages_b[0]->s_emig.f_alpha;
		break;
	case 15:
		//value = pop_emig.f_beta;
		value = stages_b[0]->s_emig.f_beta;
		break;

	default:
		cout << "that is not a valid atttribute number" << endl;
		break;
	}
	return value;
}

//Population::stage_info* Population::get_stage_info(int s) {
//	return(stages[s]);
//}

//this is only useful if i have sex-dependent depth preferences!
Population::depth_info* Population::get_depth_info(int s) {
	return sex_depths[s];
}

//Population::emig_info* Population::get_pop_emig() {
//	return &pop_emig;
//}
//
//Population::trans_info* Population::get_pop_trans() {
//	return &pop_trans;
//}
//
//Population::sett_info* Population::get_pop_sett() {
//	return &pop_sett;
//}


void Population::initialise() {

	/*default_random_engine p_generator;
	uniform_real_distribution<double> uni_dist(0.0, 1.0);*/
	
	if (mod->seedtype == 0) { //free initialisation
		if (mod->patchmodel == 1) { //if it's a patch-based model
			for (int p = 0; p < lands->patch_vector.size(); p++) {

				create_subpop(lands->patch_vector[p]);//initialise subpopulation for all possible patches, as the index in this vector should match the patch ID
				if (lands->patch_vector[p]->included_cells.size() == 0) {
					continue;
				}
				//figure out indepth cells
				int mindepth = pop_dyn.Depthmin; int maxdepth = pop_dyn.Depthmax; //get the population's depth range
				for (int c = 0; c < lands->patch_vector[p]->included_cells.size(); c++) {
					if (mindepth <= lands->patch_vector[p]->included_cells[c]->min_depth && lands->patch_vector[p]->included_cells[c]->max_depth <= maxdepth) {
						lands->patch_vector[p]->indepth_cells.push_back(lands->patch_vector[p]->included_cells[c]);
					}
				}

				lands->patch_vector[p]->calc_K("indepth_info"); //recalculate maximum individuals for how many cells are suitable
				if (mod->minX == -9) { //which will only be the case if the whole map should be initialised
					if (mod->freetype == 0) { //a proportion of suitable patches
						double pgen = unif_dist(eng); //using dist and eng defined in classes.h
						if (pgen < mod->freeprop) { //take a proportion of the suitable patches in the range
							cout << "initialised patch " << p  << endl;
							lands->patch_vector[p]->p_subpop->init_indivs(); //initialise individuals
						}
					}
					else { //all suitable patches
						lands->patch_vector[p]->p_subpop->init_indivs(); //initialise individuals
					}
				}
				else { //specified xyz extent
					if (lands->patch_vector[p]->check_limits(mod->minX, mod->maxX, mod->minY, mod->maxY, mod->minZ, mod->maxZ)) { //if patch is within the range for initialisation
						double pgen = unif_dist(eng); //using dist and eng defined in classes.h
						if (mod->freetype == 0) {//if a proportion of suitable patches
							if (pgen < mod->freeprop) { //take a proportion of the suitable patches in the range

								cout << "initialised patch " << p << endl;
								lands->patch_vector[p]->p_subpop->init_indivs(); //initialise individuals
							}
						}
						else {
							lands->patch_vector[p]->p_subpop->init_indivs(); //initialise individuals
						}
					}
				}
			}
		}
		else { //cell-based model
			for (int c = 0; c < lands->suitable_cells.size(); c++) {
				//need to check whether the limits set by the user are >=< the depth limits of the species
				//int min, max;
				//if (mod->minZ >= Depthmin) { min = mod->minZ; } else { min = Depthmin; } //if user-set minimum is deeper than species' min range, use user set minimum, otherwiser use species' minimum
				//if (mod->maxZ < Depthmax) { max = mod->maxZ; } else { max = Depthmax; } //if user-set maximum is shallower than species' range, use user set maximum, otherwise use species' maximum
				if (mod->minX == -9) { //initialise the whole map
					if (mod->freetype == 0) {
						double pgen = unif_dist(eng);
						if (pgen < mod->freeprop) { //take a proportion of suitable cells
							create_subpop(lands->suitable_cells[c]); //initialise subpop for each cell
							lands->suitable_cells[c]->c_subpop->init_indivs(); //initialise individuals
						}
					}
					else { //if all suitable cells
						create_subpop(lands->suitable_cells[c]); //initialise subpop for each cell
						lands->suitable_cells[c]->c_subpop->init_indivs(); //initialise individuals

					}
				}
				else {
					if (lands->suitable_cells[c]->check_limits(mod->minX, mod->maxX, mod->minY, mod->maxY, mod->minZ, mod->maxZ)) {
						if (mod->freetype == 0) {
							double pgen = unif_dist(eng);
							if (pgen < mod->freeprop) { //take a proportion of suitable cells
								create_subpop(lands->suitable_cells[c]); //initialise subpop for each cell
								lands->suitable_cells[c]->c_subpop->init_indivs(); //initialise individuals
							}
						}
						else { //if all suitable cells
							create_subpop(lands->suitable_cells[c]); //initialise subpop for each cell
							lands->suitable_cells[c]->c_subpop->init_indivs(); //initialise individuals

						}
					}
				}
				
			}
		}
		cout << "there are now " << subpops.size() << " subpopulations" << endl;
	}

	else if (mod->seedtype == 1) { //from species distribution
		if (mod->patchmodel == 1) { //patch-based model
			for (int p = 0; p < lands->patch_vector.size(); p++) {
				create_subpop(lands->patch_vector[p]);//initialise subpopulation
				if (lands->patch_vector[p]->in_spdist && 
					lands->patch_vector[p]->check_limits(mod->minX, mod->maxX, mod->minY, mod->maxY, mod->minZ, mod->maxZ)) { //if patch is within the range for initialisation
					if (mod->sptype == 0) { //if a proportion of suitable patches
						double pgen = unif_dist(eng);
						if (pgen < mod->spprop) { //take a proportion of the suitable patches in the range
							
							lands->patch_vector[p]->p_subpop->init_indivs(); //initialise individuals
						}
					}
					else { //if all suitable patches
						lands->patch_vector[p]->p_subpop->init_indivs(); //initialise individuals
					}
				}
			}
		}
		else { //cell-based model
			for (int c = 0; c < lands->spdist_cells.size(); c++) {
				//need to check whether the limits set by the user are >=< the depth limits of the species
				int min, max;
				if (mod->minZ > pop_dyn.Depthmin) { min = mod->minZ; }
				else { min = pop_dyn.Depthmin; } //if user-set minimum is deeper than species' range, use user set minimum, otherwiser use species' minimum
				if (mod->maxZ < pop_dyn.Depthmax) { max = mod->maxZ; }
				else { max = pop_dyn.Depthmax; } //if user-set maximum is shallower than species' range, use user set maximum, otherwise use species' maximum

				if (mod->minX == -9) { //which will only be the case if the whole map should be initialised
					if (mod->sptype == 0) { //a proportion of suitable patches
						double pgen = unif_dist(eng); //using dist and eng defined in classes.h
						if (pgen < mod->spprop) { //take a proportion of the suitable patches in the range
							cout << "initialised cell " << c << endl;
							create_subpop(lands->spdist_cells[c]); //initialise subpop for each cell
							lands->spdist_cells[c]->c_subpop->init_indivs();
						}
					}
					else { //all suitable patches
						create_subpop(lands->spdist_cells[c]); //initialise subpop for each cell
						lands->spdist_cells[c]->c_subpop->init_indivs();

					}
				}
				else{
					if (lands->spdist_cells[c]->check_limits(mod->minX, mod->maxX, mod->minY, mod->maxY, min, max)) { //if the cell is within the limits
						if (mod->sptype == 0) { //if only a proportion of suitable cells
							double pgen = unif_dist(eng);
							if (pgen < mod->spprop) { //take a proportion of suitable cells
								create_subpop(lands->spdist_cells[c]); //initialise subpop for each cell
								lands->spdist_cells[c]->c_subpop->init_indivs(); //initialise individuals
							}
						}
						else { //if all suitable cells
							cout << "creating subpop" << endl;
							create_subpop(lands->spdist_cells[c]); //initialise subpop for each cell
							lands->spdist_cells[c]->c_subpop->init_indivs(); //initialise individuals
						}
					}
				}
			}
		}

	}
	else { //from a specific patch
		//still need to initialise all the patches so there is something for the individual to enter
		for (int p = 0; p < lands->patch_vector.size(); p++) {
			create_subpop(lands->patch_vector[p]);
			//figure out indepth cells
			int mindepth = pop_dyn.Depthmin; int maxdepth = pop_dyn.Depthmax; //get the population's depth range
			for (int c = 0; c < lands->patch_vector[p]->included_cells.size(); c++) {
				if (mindepth <= lands->patch_vector[p]->included_cells[c]->min_depth && lands->patch_vector[p]->included_cells[c]->max_depth <= maxdepth) {
					lands->patch_vector[p]->indepth_cells.push_back(lands->patch_vector[p]->included_cells[c]);
				}
			}
			lands->patch_vector[p]->calc_K("indepth_info"); //recalculate maximum individuals for how many cells are suitable

		}
		cout << "initialising patch " << mod->whichpatch;
		//figure out indepth cells
		lands->patch_vector[mod->whichpatch]->p_subpop->init_indivs();
		cout << " which has " << lands->patch_vector[mod->whichpatch]->p_subpop->size << " indivs" << endl;
	}
}

//output saving functions


//float Population::mean_evolution(string param) {
//	int counter = 0;
//
//	float value = 0;
//	for (int ind = 0; ind < pop_matrix.matrix_indivs.size(); ind++) {
//		if (param == "emigration") {
//			if (pop_emig.densdep == 1) { //density dependent
//				value += pop_matrix.matrix_indivs[ind]->dinfo.D0;
//			}
//			value += pop_matrix.matrix_indivs[ind]->dinfo.emig_prob;
//		}
//		else if (param == "alpha") {
//			value += pop_matrix.matrix_indivs[ind]->dinfo.alpha;
//		}
//		else if (param == "beta") {
//			value += pop_matrix.matrix_indivs[ind]->dinfo.beta;
//		}
//		else if (param == "min_active") {
//			if (mod->size_or_time == 1) { //time depedent
//				value += pop_matrix.matrix_indivs[ind]->dinfo.min_active_time;
//			}
//			else {
//				value += pop_matrix.matrix_indivs[ind]->dinfo.min_active_size;
//			}
//		}
//		else if (param == "growth") {
//			value += pop_matrix.matrix_indivs[ind]->dinfo.m;
//		}
//		else if (param == "comp") {
//			if (mod->size_or_time == 1) { //time depedent
//				value += pop_matrix.matrix_indivs[ind]->dinfo.comp_time;
//			}
//			else {
//				value += pop_matrix.matrix_indivs[ind]->dinfo.comp_size;
//			}
//		}
//		else if (param == "settlement") {
//			value += pop_matrix.matrix_indivs[ind]->dinfo.S0;
//		}
//		else if (param == "alphaS") {
//			value += pop_matrix.matrix_indivs[ind]->dinfo.alphaS;
//		}
//		else if (param == "betaS") {
//			value += pop_matrix.matrix_indivs[ind]->dinfo.betaS;
//		}
//	}
//	value /= pop_matrix.matrix_indivs.size();
//	counter++;
//	
//	for (int s = 0; s < subpops.size(); s++) {
//		if (subpops[s]->surv_indivs.size() > 0) { //if the subpop has individuals in it
//			value+=subpops[s]->mean_evolution(param);
//			counter++;
//		}
//	}
//	value /= counter;
//
//	return value;
//
//}

void Population::write_popfile(ofstream& popfile, int rep, int year, int repseason) {
	cout << "writing popfile" << endl;
	ofstream population_file;
	stringstream ss;
	ss << mod->wrkdir + "/Outputs/Population.txt";
	popfile.open(ss.str().c_str(), std::ios_base::app);
	
	//save the annual data to the population file
	for (int s = 0; s < subpops.size(); s++) { //go through each subpop
		vector<float> res;
		vector<float> val;
		if (subpops[s]->size > 0) { //if the subpop has been settled
			//enter a new line into the population file
			popfile << rep << "\t" << year << "\t" << repseason << "\t";
			res = subpops[s]->summarise();  //get summary results for each subpopulatio
			for (int v = 0; v < res.size(); v++) {
				popfile << res[v] << "\t"; //add these to the population file
			}

			if (pop_dyn.pop_evo.ep_evo == 1) { //if emigration evolves
				//parameters will be ep, alpha, beta;
				if (year == 0) { val = subpops[s]->mean_evolution("emigration", true); }
				else { val = subpops[s]->mean_evolution("emigration", false); }
				popfile << val[0] << "\t";
				if (stages_b[0]->s_emig.densdep == 1) {
					if (year == 0) {
						popfile << val[1] << "\t" << val[2] << "\t";
					}
				}
			}
			if (pop_dyn.pop_evo.active_evo == 1 || pop_dyn.pop_evo.dv_evo==1) { //if active time/size evolves and/or min dv size/time
				
				if (year == 0) { val = subpops[s]->mean_evolution("transfer", true); }
				else { val = subpops[s]->mean_evolution("transfer", false); }
				if (pop_dyn.pop_evo.active_evo == 1) {
					popfile << val[0] << "\t";
				}
				if (pop_dyn.pop_evo.dv_evo == 1) {
					popfile << val[1] << "\t";
				}
				
			}
			if (pop_dyn.pop_evo.growth_evo == 1) { //if growth rate evolves
				//if method is linear, size of vector will be 2: m, b
				//if method is gompertz, size of vector will be 3: Linf, K, Ti
				if (year == 0) { val = subpops[s]->mean_evolution("growth", true); }
				else { val = subpops[s]->mean_evolution("growth", false); }
				popfile << val[0] << "\t" ; //only want m for now
				if(val.size()==3) { //gompertz
					popfile << val[1] << "\t" << val[2] << "\t";
				}
			}
			if (pop_dyn.pop_evo.comp_evo == 1 || pop_dyn.pop_evo.S0_evo == 1) { //if minimum competency time/size evolves
				//comp time/size, sett prob
				// if density dependent, also alphaS, betaS
				if (year == 0) { val = subpops[s]->mean_evolution("settlement", true); }
				else { val = subpops[s]->mean_evolution("settlement", false); }
				if (pop_dyn.pop_evo.comp_evo == 1) {
					popfile << val[0] << "\t";
				}
				if (pop_dyn.pop_evo.S0_evo == 1) {
					popfile << val[1] << "\t" << val[2] << "\t";
					if (val.size() == 4) {
						popfile << val[3] << "\t" << val[4] << "\t";
					}
				}
			}
			popfile << endl;
		}
		
	}

	popfile.close();
}


void Population::write_rangefile(ofstream& rangefile, int rep, int year, int repseason) {
	ofstream range_file;
	stringstream qq;
	qq << mod->wrkdir + "/Outputs/Range.txt";
	rangefile.open(qq.str().c_str(), std::ios_base::app);	
	rangefile << rep << "\t" << year << "\t" << repseason << "\t";

	float nocc=0; //number of occupied patches
	float occupsuit=0;
	float minx=-9, maxx=-9, miny=-9, maxy=-9, minz=-9, maxz=-9;
	

	if (mod->patchmodel == 1) { //if it's a patch-based model
		for (int p = 0; p < subpops.size(); p++) {
			if (subpops[p]->size > 0) {  //if the patch is currently occupied,
				nocc++; //increase nocc

				int min_col = subpops[p]->patch_num->patchlimits.xmin;
				int min_row = subpops[p]->patch_num->patchlimits.ymin;
				int max_col = subpops[p]->patch_num->patchlimits.xmax;
				int max_row = subpops[p]->patch_num->patchlimits.ymax;
				int min_depth = subpops[p]->patch_num->patchlimits.zmin;
				int max_depth = subpops[p]->patch_num->patchlimits.zmax;

				//cell_corners struct holds the coordinates of the corners of the cell

				float patch_minx = lands->raster_3d[min_row][min_col][0]->cell_corners.llc[1];//minimum x doesnt matter what y it is, so i just need to make sure i am using min_col
				float patch_maxx = lands->raster_3d[min_row][max_col][0]->cell_corners.lrc[1];
				float patch_miny = lands->raster_3d[min_row][min_col][0]->cell_corners.ulc[0]; //minimum y doesnt matter what x it is, so i just need to use min_row
				float patch_maxy = lands->raster_3d[max_row][min_col][0]->cell_corners.llc[0];
				float patch_minz = subpops[p]->patch_num->patchlimits.zmin, patch_maxz = subpops[p]->patch_num->patchlimits.zmax;

				//havent figured this out yet, because right now, the cell doesn't know how many indivs it's got! only the subpop knows, which in this case is the patch-level
				//for (int u = 0; u < lands->patch_vector[p]->included_cells.size(); u++) { //going through each cell within the patch to check how many of them are actually occupied, get an idea of depth range of indivs
				//	if (lands->patch_vector[p]->included_cells[u]->current_density > 0) {
				//		if (u = 0) { //if it's the first cell, make min and max the same
				//			patch_minz = lands->patch_vector[p]->included_cells[u]->min_depth; patch_maxz = lands->patch_vector[p]->included_cells[u]->max_depth;
				//		}
				//		else {
				//			if (lands->patch_vector[p]->included_cells[u]->min_depth < patch_minz) { patch_minz = lands->patch_vector[p]->included_cells[u]->min_depth; }
				//			if (lands->patch_vector[p]->included_cells[u]->max_depth > patch_maxz) { patch_maxz = lands->patch_vector[p]->included_cells[u]->max_depth; }
				//		}
				//	}
				//}

				if (p == 0) {
					minx = patch_minx; maxx = patch_maxx; miny = patch_miny; maxy = patch_maxy; minz = patch_minz; maxz = patch_maxz;
				}
				else {
					if (patch_minx < minx) { minx = patch_minx; }
					if (patch_maxx > maxx) { maxx = patch_maxx; }
					if (patch_miny < miny) { miny = patch_miny; }
					if (patch_maxy > maxy) { maxy = patch_maxy; }
					if (patch_minz < minz) { minz = patch_minz; }
					if (patch_maxz > maxz) { maxz = patch_maxz; }
				}
			}
		}
		occupsuit = nocc / lands->get_land_att(4); //number of occupied patches / available patches
	}

	else { //cell-based model
		for (int c = 0; c < subpops.size(); c++) {
			if (subpops[c]->size > 0) {
				nocc++;

				float cell_minx = subpops[c]->cell_num->cell_corners.llc[0];
				float cell_maxx = subpops[c]->cell_num->cell_corners.lrc[0];
				float cell_miny = subpops[c]->cell_num->cell_corners.ulc[1];
				float cell_maxy = subpops[c]->cell_num->cell_corners.llc[1];
				float cell_minz = subpops[c]->cell_num->min_depth;
				float cell_maxz = subpops[c]->cell_num->max_depth;

				if (c == 0) {
					minx = cell_minx; maxx = cell_maxx; miny = cell_miny; maxy = cell_maxy; minz = cell_minz; maxz = cell_maxz;
				}
				else {
					if (cell_minx < minx) { minx = cell_minx; }
					if (cell_maxx > maxx) { maxx = cell_maxx; }
					if (cell_miny < miny) { miny = cell_miny; }
					if (cell_maxy > maxy) { maxy = cell_maxy; }
					if (cell_minz < minz) { minz = cell_minz; }
					if (cell_maxz > maxz) { maxz = cell_maxz; }
				}
			}
		}
		occupsuit = roundf((nocc / lands->suitable_cells.size())*100)/100; //need to round it to nearest 3 decimals because otherwise the number is too small and will be outputted as 0
	}
	
	rangefile << nocc << "\t" << occupsuit << "\t" << minx << "\t" << maxx << "\t" << miny << "\t" << maxy << "\t" << minz << "\t" << maxz << endl;
	rangefile.close();
}

void Population::write_indivfile(ofstream& indfile, int rep, int year, int repseason) {//ofstream& indfile, 
	//ofstream indfile;
	stringstream tt;
	tt << mod->wrkdir + "/Outputs/Individuals.txt";
	indfile.open(tt.str().c_str(), std::ios_base::app);
	
	
	vector<int> recorded;
	//cout << "there are " << subpops.size() << " subpops" << endl;
	for (int i = 0; i < subpops.size(); i++) { //go through indivs in each subpopulation
		if (subpops[i]->surv_indivs.size() > 0) { //if the subpop has individuals in it
			for (int ind = 0; ind < subpops[i]->surv_indivs.size(); ind++) {
				indfile << rep << "\t" << year << "\t" << repseason << "\t";
				indfile << subpops[i]->surv_indivs[ind]->ID << "\t" << subpops[i]->surv_indivs[ind]->sex  << "\t" << subpops[i]->surv_indivs[ind]->status << "\t";
				if (mod->patchmodel == 1) {
					indfile << subpops[i]->surv_indivs[ind]->natal_patch << "\t";
				}
				indfile << subpops[i]->surv_indivs[ind]->m_x_movements[0] << "\t" << subpops[i]->surv_indivs[ind]->m_y_movements[0] << "\t" << subpops[i]->surv_indivs[ind]->m_z_movements[0] << "\t";
				if (mod->patchmodel == 1) { indfile << subpops[i]->surv_indivs[ind]->current_patch->patch_ID << "\t"; }
				indfile << subpops[i]->surv_indivs[ind]->current_pos.x << "\t" << subpops[i]->surv_indivs[ind]->current_pos.y << "\t" << subpops[i]->surv_indivs[ind]->current_pos.z
					<< "\t" << subpops[i]->surv_indivs[ind]->dinfo.disp_distance << "\t";
				if (mod->size_or_time == 0) {
					indfile<< subpops[i]->surv_indivs[ind]->size << "\t";
				}
				if (pop_dyn.pop_evo.ep_evo == 1) { //if emigration evolves
					if (stages_b[0]->s_emig.densdep == 1) {
						//indfile << subpops[i]->surv_indivs[ind]->dinfo.D0 << "\t";
						//indfile << subpops[i]->surv_indivs[ind]->dinfo.alpha << "\t" << subpops[i]->surv_indivs[ind]->dinfo.beta << "\t";
						if (mod->repro == 0) { //if it's asexual
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->D0[0][0] << "\t";
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->alpha[0][0] << "\t" << subpops[i]->surv_indivs[ind]->dinfo.parent->beta[0][0] << "\t";
						}
						else if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_emig.sexdep==0) { //female or not sex dependent
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->D0[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->D0[0][1])/2 << "\t";
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->alpha[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->alpha[0][1]) / 2 << "\t"
								<< (subpops[i]->surv_indivs[ind]->dinfo.parent->beta[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->beta[0][1]) / 2 << "\t";
						}
						else {
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->D0[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->D0[1][1]) / 2 << "\t";
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->alpha[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->alpha[1][1]) / 2 << "\t"
								<< (subpops[i]->surv_indivs[ind]->dinfo.parent->beta[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->beta[1][1]) / 2 << "\t";
						}
					}
					else {
						if (mod->repro == 0) { //if it's asexual
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->emig_prob[0][0] << "\t";
						}
						else if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_emig.sexdep == 0) { //female or not sex dependent
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->emig_prob[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->emig_prob[0][1] )/2 << "\t";
						}
						else if(subpops[i]->surv_indivs[ind]->sex == 0 && stages_b[0]->s_emig.sexdep == 1){
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->emig_prob[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->emig_prob[1][1]) / 2 << "\t";
						}
					}
				}
				if (pop_dyn.pop_evo.active_evo == 1) { //if active time/size evolves
					
					if (mod->size_or_time == 1) { //time-dependent
						if (mod->repro == 0) { //asexual
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_time[0][0] << "\t";
						}
						else {
							if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_trans.sexdep == 0) { //female or not sex dependent
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_time[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_time[0][1] )/2 << "\t";
							}
							else if (subpops[i]->surv_indivs[ind]->sex == 0 && stages_b[0]->s_trans.sexdep == 1) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_time[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_time[1][1]) / 2 << "\t";
							}
						}
					}
					else { //size-dependent
						if (mod->repro == 0) { //asexual
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_size[0][0] << "\t";
						}
						else {
							if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_trans.sexdep == 0) { //female or not sex dependent
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_size[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_size[0][1]) / 2 << "\t";
							}
							else if (subpops[i]->surv_indivs[ind]->sex == 0 && stages_b[0]->s_trans.sexdep == 1) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_size[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->min_active_size[1][1]) / 2 << "\t";
							}
						}
					}
				}
				if (pop_dyn.pop_evo.dv_evo == 1) { //if active time/size evolves
					if (mod->size_or_time == 1) { //time-dependent
						if (mod->repro == 0) { //if it's asexual
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_time[0][0] << "\t";
						}
						else {
							if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_trans.sexdep == 0) { //female or not sex dependent
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_time[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_time[0][1]) / 2 << "\t";
							}
							else if (subpops[i]->surv_indivs[ind]->sex == 0 && stages_b[0]->s_trans.sexdep == 1) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_time[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_time[1][1]) / 2 << "\t";
							}
						}
					}
					else { //size-dependent
						if (mod->repro == 0) { //if it's asexual
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_size[0][0] << "\t";
						}
						else {
							if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_trans.sexdep == 0) { //female or not sex dependent
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_size[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_size[0][1]) / 2 << "\t";
							}
							else if (subpops[i]->surv_indivs[ind]->sex == 0 && stages_b[0]->s_trans.sexdep == 1) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_size[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->min_dv_size[1][1]) / 2 << "\t";
							}
						}
					}
				}
				if (pop_dyn.pop_evo.growth_evo == 1) { //if growth rate evolves
					if (mod->repro == 0) { //if it's asexual
						if (subpops[i]->surv_indivs[ind]->dinfo.grow_info.method == 0) {
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->m[0][0] << "\t";
						}
						else if (subpops[i]->surv_indivs[ind]->dinfo.grow_info.method == 1) {
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->Linf[0][0] ) << "\t";
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->G_K[0][0])<< "\t";
						}
						else {

						}
					}
					else {
						if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_trans.pop_grow.sex_dep == 0) {
							if (subpops[i]->surv_indivs[ind]->dinfo.grow_info.method == 0) {
								indfile << ( subpops[i]->surv_indivs[ind]->dinfo.parent->m[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->m[0][1])/2 << "\t";
							}
							else if (subpops[i]->surv_indivs[ind]->dinfo.grow_info.method == 1) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->Linf[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->Linf[0][1]) / 2 << "\t";
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->G_K[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->G_K[0][1]) / 2 << "\t";
						
							}
							else {

							}
						}
						else if (subpops[i]->surv_indivs[ind]->sex == 0 && stages_b[0]->s_trans.pop_grow.sex_dep == 1) {
							if (subpops[i]->surv_indivs[ind]->dinfo.grow_info.method == 0) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->m[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->m[1][1]) / 2 << "\t";
							}
							else if (subpops[i]->surv_indivs[ind]->dinfo.grow_info.method == 1) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->Linf[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->Linf[1][1]) / 2 << "\t";
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->G_K[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->G_K[1][1]) / 2 << "\t";

							}
							else {

							}
						}
						

					}
				}
				if (pop_dyn.pop_evo.comp_evo == 1) { //if minimum competency time/size evolves
					if (mod->size_or_time == 1) { //time dependent
						if (mod->repro == 0) { //asexual
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->comp_time[0][0] << "\t";
						}
						else {
							if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_sett.sexdep == 0) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->comp_time[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->comp_time[0][1]) / 2 << "\t";
							}
							else if (subpops[i]->surv_indivs[ind]->sex == 0 && stages_b[0]->s_sett.sexdep == 1) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->comp_time[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->comp_time[1][1]) / 2 << "\t";
							}
						}
					}
					else { //size dependent
						if (mod->repro == 0) { //asexual
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->comp_size[0][0] << "\t";
						}
						else {
							if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_sett.sexdep == 0) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->comp_size[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->comp_size[0][1]) / 2 << "\t";
							}
							else if (subpops[i]->surv_indivs[ind]->sex == 0 && stages_b[0]->s_sett.sexdep == 1) {
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->comp_size[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->comp_size[1][1]) / 2 << "\t";
							}
						}
					}
				}
				if (pop_dyn.pop_evo.S0_evo == 1) { //if settlement rate evolves
					if (mod->repro == 0) { //asexual
						indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->S0[0][0] << "\t";
						if (stages_b[0]->s_sett.densdep == 1) { //if settlement is density dependent
							indfile << subpops[i]->surv_indivs[ind]->dinfo.parent->alphaS[0][0] << "\t" << subpops[i]->surv_indivs[ind]->dinfo.parent->betaS[0][0] << "\t";
						}
					}
					else {
						if (subpops[i]->surv_indivs[ind]->sex == 1 || stages_b[0]->s_sett.sexdep == 0) {
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->S0[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->S0[0][1])/2<< "\t";
							if (stages_b[0]->s_sett.densdep == 1) { //if settlement is density dependent
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->alphaS[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->alphaS[0][1]) / 2 << "\t";
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->betaS[0][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->betaS[0][1]) / 2 << "\t";
							}
						}
						else if (subpops[i]->surv_indivs[ind]->sex == 0 || stages_b[0]->s_sett.sexdep == 1) {
							indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->S0[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->S0[1][1]) / 2 << "\t";
							if (stages_b[0]->s_sett.densdep == 1) { //if settlement is density dependent
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->alphaS[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->alphaS[1][1]) / 2 << "\t";
								indfile << (subpops[i]->surv_indivs[ind]->dinfo.parent->betaS[1][0] + subpops[i]->surv_indivs[ind]->dinfo.parent->betaS[1][1]) / 2 << "\t";
							}
						}
					}
				}
				indfile << endl;

				recorded.push_back(subpops[i]->surv_indivs[ind]->ID);
			}
		}
	}

	if (pop_matrix.matrix_indivs.size() > 0) {

		for (int i = 0; i < pop_matrix.matrix_indivs.size(); i++) { //go through all individuals still in the matrix
			if (count(recorded.begin(), recorded.end(), pop_matrix.matrix_indivs[i]->ID)) { continue; } //if i have already recorded information for this individual this year, don't record again
			//otherwise
			indfile << rep << "\t" << year << "\t" << repseason << "\t";
			indfile << pop_matrix.matrix_indivs[i]->ID << "\t" << pop_matrix.matrix_indivs[i]->sex << "\t" << pop_matrix.matrix_indivs[i]->status << "\t";
			if (mod->patchmodel == 1) { indfile << pop_matrix.matrix_indivs[i]->natal_patch << "\t"; }
			indfile << pop_matrix.matrix_indivs[i]->m_x_movements[0] << "\t" << pop_matrix.matrix_indivs[i]->m_y_movements[0] << "\t" << pop_matrix.matrix_indivs[i]->m_z_movements[0] << "\t";
			if (mod->patchmodel == 1 && pop_matrix.matrix_indivs[i]->current_patch !=nullptr) { indfile << pop_matrix.matrix_indivs[i]->current_patch->patch_ID << "\t"; }
			else if(mod->patchmodel == 1 && pop_matrix.matrix_indivs[i]->current_patch == nullptr){ indfile << "NA" << "\t"; }
			indfile << pop_matrix.matrix_indivs[i]->current_pos.x << "\t" << pop_matrix.matrix_indivs[i]->current_pos.y << "\t" << pop_matrix.matrix_indivs[i]->current_pos.z << "\t" << pop_matrix.matrix_indivs[i]->dinfo.disp_distance << "\t";
			if (mod->size_or_time == 0) {
				indfile << pop_matrix.matrix_indivs[i]->size << "\t";
			}
			if (pop_dyn.pop_evo.ep_evo == 1) { //if emigration evolves
				if (stages_b[0]->s_emig.densdep == 1) {
					indfile << pop_matrix.matrix_indivs[i]->dinfo.D0 << "\t";
					indfile << pop_matrix.matrix_indivs[i]->dinfo.alpha << "\t" << pop_matrix.matrix_indivs[i]->dinfo.beta << "\t";
				}
				else {
					indfile << pop_matrix.matrix_indivs[i]->dinfo.emig_prob << "\t";
				}
			}
			if (pop_dyn.pop_evo.active_evo == 1) { //if active time/size evolves
				if (mod->size_or_time == 1) { //time-dependent
					indfile << pop_matrix.matrix_indivs[i]->dinfo.min_active_time << "\t";
				}
				else { //size-dependent
					indfile << pop_matrix.matrix_indivs[i]->dinfo.min_active_size << "\t";
				}
			}
			if (pop_dyn.pop_evo.dv_evo == 1) { //if active time/size evolves
				if (mod->size_or_time == 1) { //time-dependent
					indfile << pop_matrix.matrix_indivs[i]->dinfo.min_dv_time << "\t";
				}
				else { //size-dependent
					indfile << pop_matrix.matrix_indivs[i]->dinfo.min_dv_size << "\t";
				}
			}
			if (pop_dyn.pop_evo.growth_evo == 1) { //if growth rate evolves
				if (pop_matrix.matrix_indivs[i]->dinfo.grow_info.method == 0) {
					indfile << pop_matrix.matrix_indivs[i]->dinfo.grow_info.m << "\t";
				}
				else if (pop_matrix.matrix_indivs[i]->dinfo.grow_info.method == 1) {
					indfile << pop_matrix.matrix_indivs[i]->dinfo.grow_info.Linf << "\t" << pop_matrix.matrix_indivs[i]->dinfo.grow_info.G_K;
				}
				else {

				}
			}
			if (pop_dyn.pop_evo.comp_evo == 1) { //if minimum competency time/size evolves
				if (mod->size_or_time == 1) { //time dependent
					indfile << pop_matrix.matrix_indivs[i]->dinfo.comp_time << "\t";
				}
				else { //size dependent
					indfile << pop_matrix.matrix_indivs[i]->dinfo.comp_size << "\t";
				}
			}
			if (pop_dyn.pop_evo.S0_evo == 1) { //if settlement rate evolves
				indfile << pop_matrix.matrix_indivs[i]->dinfo.S0 << "\t";
				if (stages_b[0]->s_sett.densdep == 1) { //if settlement is density dependent
					indfile << pop_matrix.matrix_indivs[i]->dinfo.alphaS << "\t" << pop_matrix.matrix_indivs[i]->dinfo.betaS << "\t";
				}
			}
			indfile << endl;


		}

	}
	indfile.close();
}

void Population::write_indmovfile(ofstream& indmovfile, int rep, int year, int repseason, bool lastyear) {
	cout << "writing indmovfile...";
	stringstream ii; ofstream indmov_file;
	ii << mod->wrkdir + "/Outputs/Indiv_Mov_rep" << rep << ".txt";
	indmovfile.open(ii.str().c_str(), std::ios_base::app);
	vector<int> recorded;
	//cout << "going through indiv file: ";
	for (int i = 0; i < subpops.size(); i++) { //go through indivs in each subpopulation
		if (subpops[i]->surv_indivs.size() > 0) { //if the subpop has individuals in it
			//cout << "surviving indivs" << endl;
			for (int ind = 0; ind < subpops[i]->surv_indivs.size(); ind++) {
				if (subpops[i]->surv_indivs[ind]->emig == false) { continue; } //i only want to know the paths of those who emigrated
				//if (lastyear==false) { //if it isnt the last year, only tell me the paths of dead indivs?
				//	if (subpops[i]->surv_indivs[ind]->dead==false) { cout << "writing indmov file, indiv still alive, skipped" << endl; continue; } //skip those that havent died
				//}
				//cout << "writing indmov file, indiv dead, writing data" << endl;
				indmovfile << rep << "\t" << year << "\t" << repseason << "\t";
				indmovfile << subpops[i]->surv_indivs[ind]->ID << "\t" << subpops[i]->surv_indivs[ind]->sex << "\t" << subpops[i]->surv_indivs[ind]->status << "\t";
				//cout << "ID and status "; 
				
				for (int x = 0; x < subpops[i]->surv_indivs[ind]->m_x_movements.size(); x++) {
					if (x == subpops[i]->surv_indivs[ind]->m_x_movements.size() - 1) { indmovfile << subpops[i]->surv_indivs[ind]->m_x_movements[x] << endl; } //if it's the last coordinate, end the line
					else { indmovfile << subpops[i]->surv_indivs[ind]->m_x_movements[x] << "\t"; } //otherwise just separate with tabs
				}
				
				indmovfile << rep << "\t" << year << "\t" << repseason << "\t" << subpops[i]->surv_indivs[ind]->ID << "\t" << subpops[i]->surv_indivs[ind]->sex << "\t" << subpops[i]->surv_indivs[ind]->status << "\t";

				for (int y = 0; y < subpops[i]->surv_indivs[ind]->m_y_movements.size(); y++) {

					if (y == subpops[i]->surv_indivs[ind]->m_y_movements.size() - 1) { indmovfile << subpops[i]->surv_indivs[ind]->m_y_movements[y] << endl; } //if it's the last coordinate, end the line
					else { indmovfile << subpops[i]->surv_indivs[ind]->m_y_movements[y] << "\t"; } //otherwise just separate with tabs
				}
				
				indmovfile << rep << "\t" << year << "\t" << repseason << "\t" << subpops[i]->surv_indivs[ind]->ID << "\t" << subpops[i]->surv_indivs[ind]->sex << "\t"  << subpops[i]->surv_indivs[ind]->status << "\t";

				for (int z = 0; z < subpops[i]->surv_indivs[ind]->m_z_movements.size(); z++) {

					if (z == subpops[i]->surv_indivs[ind]->m_z_movements.size() - 1) { indmovfile << subpops[i]->surv_indivs[ind]->m_z_movements[z] << endl; } //if it's the last coordinate, end the line
					else { indmovfile << subpops[i]->surv_indivs[ind]->m_z_movements[z] << "\t"; } //otherwise just separate with tabs
				}
				
				recorded.push_back(subpops[i]->surv_indivs[ind]->ID);
			}
		}
	}

	for (int i = 0; i < pop_matrix.matrix_indivs.size(); i++) { //go through all individuals still in the matrix
		//i don't check for last year here becaus ethe only indivs in the matrix will die because they've run out of time
		if (pop_matrix.matrix_indivs[i]->emig == false) { continue; } //i only want to know the paths of those who emigrated

		if (count(recorded.begin(), recorded.end(), pop_matrix.matrix_indivs[i]->ID)) { continue; } //if i have already recorded information for this individual this year, don't record again
		//otherwise
		indmovfile << rep << "\t" << year << "\t" << repseason << "\t";
		indmovfile << pop_matrix.matrix_indivs[i]->ID << "\t" << pop_matrix.matrix_indivs[i]->sex << "\t" << pop_matrix.matrix_indivs[i]->status << "\t";

		//cout << "ID and status ";
		for (int x = 0; x < pop_matrix.matrix_indivs[i]->m_x_movements.size(); x++) {
			
			if (x == pop_matrix.matrix_indivs[i]->m_x_movements.size() - 1) { indmovfile << pop_matrix.matrix_indivs[i]->m_x_movements[x] << endl; } //if it's the last coordinate, end the line
			else { indmovfile << pop_matrix.matrix_indivs[i]->m_x_movements[x] << "\t"; } //otherwise just separate with tabs
		}
		
		indmovfile << rep << "\t" << year << "\t" << repseason << "\t" << pop_matrix.matrix_indivs[i]->ID << "\t" << pop_matrix.matrix_indivs[i]->sex << "\t" << pop_matrix.matrix_indivs[i]->status << "\t";

		for (int y = 0; y < pop_matrix.matrix_indivs[i]->m_y_movements.size(); y++) {
			if (y == pop_matrix.matrix_indivs[i]->m_y_movements.size() - 1) { indmovfile << pop_matrix.matrix_indivs[i]->m_y_movements[y] << endl; } //if it's the last coordinate, end the line
			else { indmovfile << pop_matrix.matrix_indivs[i]->m_y_movements[y] << "\t"; } //otherwise just separate with tabs
		}
		
		indmovfile << rep << "\t" << year << "\t" << repseason << "\t" << pop_matrix.matrix_indivs[i]->ID << "\t" << pop_matrix.matrix_indivs[i]->sex << "\t" << pop_matrix.matrix_indivs[i]->status << "\t";

		for (int z = 0; z < pop_matrix.matrix_indivs[i]->m_z_movements.size(); z++) {
			if (z == pop_matrix.matrix_indivs[i]->m_z_movements.size() - 1) { indmovfile << pop_matrix.matrix_indivs[i]->m_z_movements[z] << endl; } //if it's the last coordinate, end the line
			else { indmovfile << pop_matrix.matrix_indivs[i]->m_z_movements[z] << "\t"; } //otherwise just separate with tabs
		}
		

	}
	cout << "done" << endl;
	indmovfile.close();
}

//run the model
void Population::one_year(ofstream& popfile, ofstream& indfile, ofstream& indmovfile, ofstream& rangefile, int rep, int year) {
	cout << "one year called" << endl;
	for (int r = 0; r < pop_dyn.Nrepro; r++) { //for every reproductive season in a year
		//if (year==0 && r==0) { // i want to know what the stats are before anything happens but after indivs have been initialised
		//	write_popfile(popfile, rep, 0, 0);
		//}
		if (year >= 800) {
			write_popfile(popfile, rep, year, pop_dyn.Nrepro - 1);//since it's the end of year, it's the last repro season

		}
		else if (mod->Outpop == 1 && year%mod->outpop_int == 0) { //if we want the population file output, and it's in the right interval
			//ofstream population_file;
			//stringstream ss;
			//ss << mod->wrkdir + "/Outputs/Population.txt";
			//population_file.open(ss.str().c_str());
			//write_popfile(population_file, rep, year, pop_dyn.Nrepro - 1);//since it's the end of year, it's the last repro season
			//population_file.close();

			write_popfile(popfile, rep, year, pop_dyn.Nrepro - 1);//since it's the end of year, it's the last repro season

		}
		cout << "wrote popfile ";
		if (mod->Outrange == 1 && year%mod->outrange_int == 0) { //if we want range output and the year is in the right interval
			//ofstream range_file;
			//stringstream qq;
			//qq << mod->wrkdir + "/Outputs/Range.txt";
			//range_file.open(qq.str().c_str());
			//write_rangefile(range_file, rep, year, pop_dyn.Nrepro - 1); //since it's the end of year, it's the last repro season
			//range_file.close();
			
			write_rangefile(rangefile, rep, year, pop_dyn.Nrepro - 1); //since it's the end of year, it's the last repro season

		}
		
		cout << "size of subpops is " << subpops.size() << endl;
		
		//reproduction : need to let all subpops reproduce before they emigrate
		for (int s = 0; s < subpops.size(); s++) {
			if (subpops[s]->size>0) { subpops[s]->reproduction(); }
		}
		cout << "did repro" << endl;

		//emigration
		for (int s = 0; s < subpops.size(); s++) {
			if (subpops[s]->size > 0) {
				subpops[s]->emigration();
			}
			
		}

		//transfer +settlement
		cout << "starting nindivs in matrix is " << pop_matrix.matrix_indivs.size() << endl;
		int pld_time = 0, max_pld = 0;
		while (pop_matrix.total > 0 ) { //as long as there are still dispersers in the matrix
			
			for (int i = 0; i < pop_matrix.matrix_indivs.size(); i++) { //go through each individual
				//if (pld_time == 0) {
				//	//cout << "leaving natal" << endl;
				//	pop_matrix.matrix_indivs[i]->leave_natal();
				//	continue;
				//}
				
				//else if (pld_time > pop_matrix.matrix_indivs[i]->dinfo.pld && pop_matrix.matrix_indivs[i]->status != 4) { pop_matrix.matrix_indivs[i]->dead = true; pop_matrix.matrix_indivs[i]->status = 4; cout << "indiv ran out of time" << endl; pop_matrix.total--; continue; } //if they've exceeded their pld time, they die
				//if (pld_time > pop_matrix.matrix_indivs[i]->dinfo.pld && pop_matrix.matrix_indivs[r]->dead == false && pop_matrix.matrix_indivs[r]->sett==false) {
				//	//if they've run out of pld time but they're not dead and they haven't settled
				//	pop_matrix.matrix_indivs[i]->dead = true; pop_matrix.matrix_indivs[i]->status = 4; cout << "indiv ran out of time" << endl; pop_matrix.total--; continue;
				//}
				if (max_pld == 0) { max_pld = pop_matrix.matrix_indivs[i]->dinfo.pld; cout << "max pld is " << max_pld << endl; }
				if (pld_time > max_pld && pop_matrix.matrix_indivs[i]->status==1) { //still dispersing
					//if they've run out of pld time but they're not dead and they haven't settled
					pop_matrix.matrix_indivs[i]->dead = true; pop_matrix.matrix_indivs[i]->status = 10; /*cout << "indiv ran out of time" << endl;*/ 
					pop_matrix.total--; continue;
				}
				//else if(pld_time> pop_matrix.matrix_indivs[i]->dinfo.comp){ pop_matrix.matrix_indivs[i]->dinfo.SL}
				
				else if(pop_matrix.matrix_indivs[i]->sett == true || pop_matrix.matrix_indivs[i]->dead == true) {  continue; } //if the individual has settled already, or has died, skip it. since i go through the whole list every time, this will happen
				
				else if (pop_matrix.matrix_indivs[i]->transfer_to_settle()==false) { //let it take a step, update current position and current cell. if it settles or dies,
					
					if(pop_matrix.matrix_indivs[i]->dead==true){
						pop_matrix.total--; //reduce the number of dispersers in the matrix by 1
						continue;
					}
					if (mod->patchmodel == 1) { //if it's a patch based model, i've already initialised all patches so dont need to create more subpops
						//cout << "preparing to enter subpop" << pop_matrix.matrix_indivs[i]->current_patch->patch_ID << endl;
						//cout << "preparing to enter subpop" << pop_matrix.matrix_indivs[i]->current_cell->patch->patch_ID << endl;
						//pop_matrix.matrix_indivs[i]->current_patch->p_subpop->enter(pop_matrix.matrix_indivs[i]); //enter new subpopulation if it hasnt died in transit. this will -- from matrix
						pop_matrix.matrix_indivs[i]->current_cell->patch->p_subpop->enter(pop_matrix.matrix_indivs[i]); //enter new subpopulation if it hasnt died in transit. this will -- from matrix

					}
					else { //cell-based model
						if (pop_matrix.matrix_indivs[i]->current_cell->c_subpop) { //if the subpop has already been created
							pop_matrix.matrix_indivs[i]->current_cell->c_subpop->enter(pop_matrix.matrix_indivs[i]); //enter new subpopulation , this wil -- from matrix
						}
						else { //need to create a new subpop!
							create_subpop(pop_matrix.matrix_indivs[i]->current_cell);
							pop_matrix.matrix_indivs[i]->current_cell->c_subpop->enter(pop_matrix.matrix_indivs[i]); //enter new subpopulation, this will --from matrix
						}
						
					}
					//entering a subpop does pop_matrix.total--
				}
				//cout << "there are now " << pop_matrix.total << " indivs in the matrix" << endl;
				//otherwise move on to next individual
				//if (pop_matrix.matrix_indivs[i]->sett || pop_matrix.matrix_indivs[i]->dead) { cout << "there are now " << pop_matrix.matrix_indivs.size() << "indivs in the vector and " << pop_matrix.total << " total" << endl; }
				pop_matrix.matrix_indivs[i]->dinfo.disp_time++; //increment the number of hours an indiv has been dispersing for
				if (mod->size_or_time == 1) { //if behaviour is time-dependent
					if (pop_matrix.matrix_indivs[i]->dinfo.diel_vert == 1 && pld_time >= pop_matrix.matrix_indivs[i]->dinfo.min_dv_time) { //if diel vertical migration happens before it becomes active, and enough time has passed
						pop_matrix.matrix_indivs[i]->dvm = true;
						//cout << "made dvm true" << endl;
					}
					if (pop_matrix.matrix_indivs[i]->dinfo.mode == 1 && //if indiv is a hybrid passive then active disperser
						pop_matrix.matrix_indivs[i]->dinfo.phase == 0 && //if it is currently still passive 
						pld_time >= pop_matrix.matrix_indivs[i]->dinfo.min_active_time) { //if the minimum time before active swimming has passed
						
						pop_matrix.matrix_indivs[i]->dinfo.SL = pop_matrix.matrix_indivs[i]->dinfo.active_SL; //give it its new step length
						pop_matrix.matrix_indivs[i]->dinfo.phase == 1; //change its phase to active
						if (pop_matrix.matrix_indivs[i]->dinfo.dv_active == 0) {
							pop_matrix.matrix_indivs[i]->dvm = false;
						}
						//cout << "past min active time, updating SL" << endl;
					}
					if (pop_matrix.matrix_indivs[i]->dinfo.comp_time !=-9 && pld_time >= pop_matrix.matrix_indivs[i]->dinfo.comp_time) { //if the minimum time before it is able to settle has passed
						pop_matrix.matrix_indivs[i]->dinfo.comp = true; //it is now competent for settlement
						//cout << "now competent" << endl;
					}

				}
				//size-dependent behaviour changes happen in the growth function
				else {
					pop_matrix.matrix_indivs[i]->larval_growth();
				}
			}
			pld_time++;
			/*cout << "elapsed pld time is " << pld_time << " hours " << endl;
			cout << "nindivs in matrix: " << pop_matrix.total << endl;*/
			if (pld_time % 24 == 0) { //if a whole day has passed 
				cout << "elapsed pld time is " << pld_time / 24 << "day "<< endl;
				cout << "nindivs in matrix: " << pop_matrix.total << endl;
				if ((pld_time / 24) > (max_pld/24)) {
					cout << "pld time is " << pld_time << " but there were remaining indivs in matrix, overwrite" << endl;
					pop_matrix.total = 0;
					//cout << "pld of indiv is " << pop_matrix.matrix_indivs[0]->dinfo.pld << " and their status us " << pop_matrix.matrix_indivs[0]->status << " have they settled? " << pop_matrix.matrix_indivs[0]->sett << endl;
				}
			}

		}
		//pop_matrix.matrix_indivs.clear(); //once there are no more dispersers as they have all settled or died, clear the matrix population
		//cout << "cleared matrix indivs" << endl;
		
		
		//save to file between reproductive events
		if (r < pop_dyn.Nrepro - 1 && pop_dyn.Nrepro > 1) {//if there is more than one reproductive season and this is not the last one,
			cout << "writing files between repro events" << endl;
			//write indiv file after settlement, before survival
			if (mod->Outind == 1 && year%mod->outind_int == 0) { //if we want individual output and the year is in the right interval
				cout << "calling write_indiv" << endl;
				/*ofstream indivs_file;
				stringstream tt;
				tt << mod->wrkdir + "/Outputs/Individuals.txt";
				indivs_file.open(tt.str().c_str());
				write_indivfile(indivs_file, rep, year , r);
				indivs_file.close();*/

				write_indivfile(indfile, rep, year, r);
			}
			if (mod->Outindmov == 1 && year%mod->outind_int == 0) { //if we want individual output and the year is in the right interval
				cout << "calling outindmove" << endl;
				//ofstream indmov_file;
				/*stringstream ii;
				ii << mod->wrkdir + "/Outputs/Indiv_Mov_rep" << r << ".txt";
				indmov_file.open(ii.str().c_str());*/
				if (year == mod->nyears - 1) {
					write_indmovfile(indmovfile, rep, year, r, true); //since it's the last year, write everyone's paths, not just the dead indivs
					//write_indmovfile(indmov_file, rep, year , r, true); //since it's the last year, write everyone's paths, not just the dead indivs
				}
				else { 
					write_indmovfile(indmovfile, rep, year, r, false);
					//write_indmovfile(indmov_file, rep, year , r, false); 
				} //otherwise write only the indivs that have died
				//indmov_file.close();
			}
			//if (mod->Outpop == 1 && year%mod->outpop_int == 0) { // if we want population output and the year is in the right interval
			//	write_popfile(popfile, rep, year , r);
			//}
			//if (mod->Outrange == 1 && year%mod->outrange_int == 0) { //if we want range output and the year is in the right interval
			//	write_rangefile(rangefile, rep, year , r);
			//}

			//survival between reproductive events
			if (r < pop_dyn.Nrepro - 1 && pop_dyn.survsched == 0) { //only do survival until the last reproductive season in the year
				for (int s = 0; s < subpops.size(); s++) {
					cout << "between rep survival" << endl;
					subpops[s]->survival();
				}
				//if (mod->stagestruct == 1) { //if there are stages
				//	for (int s = 0; s < subpops.size(); s++) { //development between reproduction in case some stages need less than a year to be reproductively active (ie minimum age for several stages is 0)
				//		subpops[s]->development(false); //because aging is false, it will not increase the age of individuals
				//	}
				//}
			}
			pop_matrix.matrix_indivs.clear(); //once there are no more dispersers as they have all settled or died, clear the matrix population
			pop_matrix.total = 0;

		}
		//indivs cant stay in the matrix longer than their dispersal period
		
	}
	

	//annual survival/survival after last reproductive season (so that i can include aging)
	cout << "onto annual survival" << endl;
	for (int s = 0; s < subpops.size(); s++) {
		subpops[s]->survival();
		if (mod->stagestruct == 1) { subpops[s]->development(true); } //this will age the individuals, if the population is stage-structured
	}
	cout << "finished annual survival" << endl;
	//write end of year outputs
	//if (mod->Outpop == 1 && year%mod->outpop_int == 0) { //if we want the population file output, and it's in the right interval
	//	write_popfile(popfile, rep, year+1, pop_dyn.Nrepro-1);//since it's the end of year, it's the last repro season
	//}
	//cout << "wrote popfile ";
	//if (mod->Outrange == 1 && year%mod->outrange_int == 0) { //if we want range output and the year is in the right interval
	//	write_rangefile(rangefile, rep, year+1, pop_dyn.Nrepro-1); //since it's the end of year, it's the last repro season
	//}
	//cout << "wrote rangefile ";

	//this is specifically for the disturbance testing in chapter 2 evolution
	/*if (year >= 800) {
		write_indivfile(indfile, rep, year, pop_dyn.Nrepro - 1); 

	}*/
	if (year >= 800) {
		write_indivfile(indfile, rep, year, pop_dyn.Nrepro - 1); //since it's the end of year, it's the last repro season

	}
	else if (mod->Outind == 1 && year%mod->outind_int == 0) { //if we want individual output and the year is in the right interval
		/*ofstream indivs_file;
		stringstream tt;
		tt << mod->wrkdir + "/Outputs/Individuals.txt";
		indivs_file.open(tt.str().c_str());*/
		//write_indivfile(indivs_file, rep, year, pop_dyn.Nrepro - 1); //since it's the end of year, it's the last repro season
		//indivs_file.close();
		write_indivfile(indfile, rep, year, pop_dyn.Nrepro - 1); //since it's the end of year, it's the last repro season
	}
	cout << "wrote indfile ";

	if (mod->Outindmov == 1 && year%mod->outind_int == 0) { //if we want individual output and the year is in the right interval
		/*ofstream indmov_file;
		stringstream ii;
		ii << mod->wrkdir + "/Outputs/Indiv_Mov_rep" << rep << ".txt";
		indmov_file.open(ii.str().c_str());*/
		if (year == mod->nyears - 1) {
			//write_indmovfile(indmov_file, rep, year, pop_dyn.Nrepro - 1, true); //since it's the last year, write everyone's paths, not just the dead indivs

			write_indmovfile(indmovfile, rep, year, pop_dyn.Nrepro - 1, true); //since it's the last year, write everyone's paths, not just the dead indivs
		}
		else{ 
			//write_indmovfile(indmov_file, rep, year, pop_dyn.Nrepro - 1, false);
			write_indmovfile(indmovfile, rep, year, pop_dyn.Nrepro - 1, false); 
		} //otherwise write only the indivs that have died
		
		//indmov_file.close();
	}
	cout << "wrote indmovfile ";
	if (lands->patch_extinct == 1 ) {
		if (year >= (lands->p_ext_burnin) && (year - lands->p_ext_burnin) % lands->p_ext_int == 0) {
			cout << " local patch extinction" << endl;
			float ran_i_ext;
			if (lands->p_ext_method == 0) { //random proportion
				
				for (int p = 0; p < subpops.size(); p++) {
					float ran_p_ext = unif_dist(eng);
					if (subpops[p]->size > 0 && ran_p_ext < lands->p_ext_p_prop) { //if the subpop isn't empty, and it is going to go extinct	
						cout << "clearing patch " << p;
						for (int i = 0; i < subpops[p]->surv_indivs.size(); i++) {
							ran_i_ext = unif_dist(eng);
							if (subpops[p]->surv_indivs[i]->dead == false && ran_i_ext <= lands->p_ext_i_prop) {
								subpops[p]->surv_indivs[i]->dead = true;
								subpops[p]->size--;
							}
						}
						cout << ", size of pop is " << subpops[p]->size << endl;
					}
				}
			}
			else { //specific patch					
				if (subpops[lands->p_ext_patch]->size > 0) { //if the subpop isn't empty, and it is going to go extinct
					cout << "clearing patch " << lands->p_ext_patch << "which had " << subpops[lands->p_ext_patch]->size << "indivs " << endl;
					for (int i = 0; i < subpops[lands->p_ext_patch]->surv_indivs.size(); i++) {
						ran_i_ext = unif_dist(eng);
						if (subpops[lands->p_ext_patch]->surv_indivs[i]->dead==false && ran_i_ext <= lands->p_ext_i_prop) {
							
							subpops[lands->p_ext_patch]->surv_indivs[i]->dead = true;
							subpops[lands->p_ext_patch]->size--;
						}

					}
					cout << "there should now be " << subpops[lands->p_ext_patch]->size << " indivs in this subpop" << endl;
				}
			}
		}
		
	}

	cout << "cleaning up matrix indivs for next year, there are " <<  pop_matrix.matrix_indivs.size() << " indivs" << endl;
	
	for (int p = 0; p < pop_matrix.matrix_indivs.size(); p++) {
		/*if (pop_matrix.matrix_indivs[p]->status == 1 || pop_matrix.matrix_indivs[p]->status == 4 || 
			pop_matrix.matrix_indivs[p]->status == 5 || pop_matrix.matrix_indivs[p]->status == 8 || pop_matrix.matrix_indivs[p]->status == 10) {
			delete(pop_matrix.matrix_indivs[p]);
		}*/
		
		if (pop_matrix.matrix_indivs[p] == nullptr) { //if it has already been deleted in subpop cleanup
			cout << "clearing matrix indivs, an indiv is nullptr so i'm continuing..." << endl;
			continue;
		}
		//cout << "indiv " << p << " has status " << pop_matrix.matrix_indivs[p]->status << endl;
		/*if (pop_matrix.matrix_indivs[p]->status == 3 || pop_matrix.matrix_indivs[p]->status == 9 || pop_matrix.matrix_indivs[p]->status == 11) {
			continue;
		}*/
		if (pop_matrix.matrix_indivs[p]->current_patch!=nullptr) {
			continue;
		}
		else {
			delete(pop_matrix.matrix_indivs[p]);
		}
		//if (pop_matrix.matrix_indivs[p]->dead == true || pop_matrix.matrix_indivs[p]->patch_num == nullptr || pop_matrix.matrix_indivs[p]->sett == false) {
		//	//if it has died, or hasn't settled, delete it from memory and deallocate.
		//	//the only individuals that should be safe are those that settled somewhere
		//	delete(pop_matrix.matrix_indivs[p]);
		//}
		//1 = disperser, 2 = local recruit, 3 = settled in suitable habitat, 4 = forced to settle on unsuitable habitat and died,
			// 5=died by dispersal-related mortality, 6= died by annual mortality, 7=died by exceeding max age, 8=died by being absorbed, 
			// 9=forced settlement(passive) on suitable habitat, 10= ran out of time,
	}
	cout << "cleaning up subpops ready for next year" << endl;
	for (int s = 0; s < subpops.size(); s++) {
		subpops[s]->clean_up();
	}
	cout << "clearing " << endl;
	pop_matrix.matrix_indivs.clear(); //once there are no more dispersers as they have all settled or died, clear the matrix population
	pop_matrix.total = 0;

}


//delete subpops
void Population::delete_subpops() {
	for (int s = 0; s < subpops.size(); s++) {
		subpops[s]->delete_indivs();
		subpops[s]->rep_males.clear(); subpops[s]->rep_females.clear();
		if (mod->patchmodel == 1) { subpops[s]->patch_num->p_subpop = nullptr; }
		else { subpops[s]->cell_num->c_subpop = nullptr; }
		//subpops[s]->size = subpops[s]->indivs.size();
		delete subpops[s];
	}
	subpops.clear(); //clear the vector 

	//i think i did this already but there's no harm in doing it again
	/*for (int m = 0; m < pop_matrix.matrix_indivs.size(); m++) {
		delete pop_matrix.matrix_indivs[m];
	}
	pop_matrix.matrix_indivs.clear();
	pop_matrix.total = 0;*/

}


////NOT DONE YET!
//reset population for next rep
void Population::reset() {
	for (int p = 0; p < subpops.size(); p++) {
		if (mod->patchmodel == 1) {
			subpops[p]->patch_num->indepth_cells.clear(); //because i am initialising them again 

		}
	}
	
	delete_subpops(); //this gets rid of all individuals in each subpop, 

	if (pop_matrix.matrix_indivs.size() != 0) { cout << "there are still indivs in the matrix" << endl; }
	
	pop_dyn.total_births = 0; //reset indiv IDs
}