#include "classes.h"
////////////////CONSTRUCTORS ////////////////
//default constructor
Individual::Individual() {

}

Individual::~Individual()
{
}

//non- stage structured individual constructor
//Individual::Individual(Model* mo, pop_info* p, Landscape* la, emig_info* e, trans_info* t, sett_info* s, int x, int IDnum) {
//	emig = false; //has it emigrated yet?
//	dvm = false; //is it old/big enough to undergo diel vertical migration?
//	sett = false; //has it settled yet?
//	//rep = false; //reproduced yet?
//	dead = false; //has it died?
//	age = 0;
//	pstage = 0;
//	sex = x;
//	ID = IDnum;
//	rep_int = 0; //start at 0, will increment if population RepInt>0
//	skip_breeding = false;
//
//	mod = mo;
//	lands = la;
//	pop_dyn = p;
//	/*pop_einfo = e;
//	pop_tinfo = t;
//	pop_sinfo = s;*/
//
//	set_emig_info();
//	set_trans_info();
//	set_settle_info();
//}

//stage structured individual constructor
Individual::Individual(Model* mo, pop_info* p, Landscape* la, /*emig_info* e, trans_info* t, sett_info* s, */stage_infob* st, int x, int IDnum, bool invent, vector<stage_infob*>* emigst) {
	status = 0;
	emig = false; //has it emigrated yet?
	dvm = false; //is it old/big enough to undergo diel vertical migration?
	sett = false; //has it settled yet?
	//rep = false; //reproduced yet?
	dead = false; //has it died?
	age = 0;
	sex = x; //by default, each indiv is female
	ID = IDnum;
	rep_int = 0; //start at 0, will increment if population RepInt>0
	skip_breeding = false;

	mod = mo;
	lands = la;
	pop_dyn = p;
	//pop_einfo = e;
	//pop_tinfo = t;
	//pop_sinfo = s;
	pstage = st;
	subpop_stages = emigst;
	dinfo.assign_parent();
	//if (mod->repro != 0) { //if it's a sexual model
	//	dinfo.assign_parent();
	//}

	//initial size
	if (invent == true && mod->size_or_time==0 && mod->repro > 0) { //if parental information needs to be invented and model is size-dependent and sexual reproduction 
		//this can't be sex-dependent from input
		dinfo.parent->size_at_birth[0][0] = pop_dyn->juv_size; //assign mother's allele
		dinfo.parent->size_at_birth[0][1] = pop_dyn->juv_size; //assign father's allele
		size_at_birth = (dinfo.parent->size_at_birth[0][0] + dinfo.parent->size_at_birth[0][1]) / 2; //take the average of those two values  (redundant because they're the same value)
		size = size_at_birth;
	}
	
	set_emig_info(invent);
	set_trans_info(invent);
	set_settle_info(invent);

}


//general use functions
float Individual::get_localK() {
	float max_inds;
	if (mod->patchmodel == 1) { //patch model
		max_inds = current_patch->indepth_info.max_inds; 
	}
	else { //cell based
		max_inds = current_cell->max_inds;
	}
	return max_inds;
}

int Individual::get_localN() {
	int N;
	if (mod->patchmodel == 1) { N = current_patch->p_subpop->size; }
	else { N = current_cell->c_subpop->size; }
	return N;
}

//**************reproduction************//

//give new individuals the information from their parents
void Individual::disp_gene_expression() { //this function is called when individual is in the right stage to disperse so can use pstage.
	if (mod->repro > 0) { //if it's sexual reproduction, therefore inheritance came from both parents
		int n_sexes;
		if (pstage->s_emig.sexdep == 1) { n_sexes = 2; } //need separate values for both sexes
		else { n_sexes = 1; } //same value for both sexes, so no need to keep it separate
	

		//here it is also important whether it is sex dependent
		//if sex INDEPENDENT, then there is 1 locus with 2 alleles. the offspring expresses the AVERAGE of the parents' values
		//if sex dependent then there are 2 loci with 2 alleles and indiv randomy inherits one from each parent and expresses the average
		//remember the full matrix looks like this:
		//			Mother		Father
		// female	
		// male

		//emigration parameters
		for (int s = 0; s < n_sexes; s++) {
			if (pstage->s_emig.densdep == 0) {

				if (s == 1 && sex == 0) { //for a male
					dinfo.emig_prob = (dinfo.parent->emig_prob[s][0] + dinfo.parent->emig_prob[s][1]) / 2; //take the average of those two values
				}
				else {
					//calculate mean
					dinfo.emig_prob = (dinfo.parent->emig_prob[0][0] + dinfo.parent->emig_prob[0][1]) / 2; //take the average of those two values
				}
			}
			else {

				if (s == 1 && sex == 0) { //if it's for a male
					dinfo.D0 = (dinfo.parent->D0[s][0] + dinfo.parent->D0[s][1]) / 2;
					dinfo.alpha = (dinfo.parent->alpha[s][0] + dinfo.parent->alpha[s][1]) / 2;
					dinfo.beta = (dinfo.parent->beta[s][0] + dinfo.parent->beta[s][1]) / 2;
				}
				else {
					//calculate mean
					dinfo.D0 = (dinfo.parent->D0[0][0] + dinfo.parent->D0[0][1]) / 2;
					dinfo.alpha = (dinfo.parent->alpha[0][0] + dinfo.parent->alpha[0][1]) / 2;
					dinfo.beta = (dinfo.parent->beta[0][0] + dinfo.parent->beta[0][1]) / 2;
				}

			}

			//transfer
			if (mod->size_or_time == 0) { //size dependent, so need the slope

				if (s == 1 && sex == 0) {
					dinfo.SL_m = (dinfo.parent->SL_m[s][0] + dinfo.parent->SL_m[s][1]) / 2; //take the average of those two values
					dinfo.rho_m = (dinfo.parent->rho_m[s][0] + dinfo.parent->rho_m[s][1]) / 2;
					dinfo.min_active_size = (dinfo.parent->min_active_size[s][0] + dinfo.parent->min_active_size[s][1]) / 2;
					dinfo.min_dv_size = (dinfo.parent->min_dv_size[s][0] + dinfo.parent->min_dv_size[s][1]) / 2;
					dinfo.grow_info.m = (dinfo.parent->m[s][0] + dinfo.parent->m[s][1]) / 2;
					dinfo.grow_info.Linf = (dinfo.parent->Linf[s][0] + dinfo.parent->Linf[s][1]) / 2;
				}
				else {
					// calculate mean
					dinfo.SL_m = (dinfo.parent->SL_m[0][0] + dinfo.parent->SL_m[0][1]) / 2; //take the average of those two values
					dinfo.rho_m = (dinfo.parent->rho_m[0][0] + dinfo.parent->rho_m[0][1]) / 2;
					dinfo.min_active_size = (dinfo.parent->min_active_size[0][0] + dinfo.parent->min_active_size[0][1]) / 2;
					dinfo.min_dv_size = (dinfo.parent->min_dv_size[0][0] + dinfo.parent->min_dv_size[0][1]) / 2;
					dinfo.grow_info.m = (dinfo.parent->m[0][0] + dinfo.parent->m[0][1]) / 2;
					dinfo.grow_info.Linf = (dinfo.parent->Linf[0][0] + dinfo.parent->Linf[0][1]) / 2;
				}
			}
			else {

				
				if (s == 1 && sex == 0) { //for a male
					dinfo.SL = (dinfo.parent->SL[s][0] + dinfo.parent->SL[s][1]) / 2; //take the average of those two values
					dinfo.active_SL = (dinfo.parent->active_SL[s][0] + dinfo.parent->active_SL[s][1]) / 2;
					dinfo.rho = (dinfo.parent->rho[s][0] + dinfo.parent->rho[s][1]) / 2;
					dinfo.min_active_time = (dinfo.parent->min_active_time[s][0] + dinfo.parent->min_active_time[s][1]) / 2;
					dinfo.min_dv_time = (dinfo.parent->min_dv_time[s][0] + dinfo.parent->min_dv_time[s][1]) / 2;

				}
				else {
					//calculate mean
					dinfo.SL = (dinfo.parent->SL[0][0] + dinfo.parent->SL[0][1]) / 2; //take the average of those two values
					dinfo.active_SL = (dinfo.parent->active_SL[0][0] + dinfo.parent->active_SL[0][1]) / 2;
					dinfo.rho = (dinfo.parent->rho[0][0] + dinfo.parent->rho[0][1]) / 2;
					dinfo.min_active_time = (dinfo.parent->min_active_time[0][0] + dinfo.parent->min_active_time[0][1]) / 2;
					dinfo.min_dv_time = (dinfo.parent->min_dv_time[0][0] + dinfo.parent->min_dv_time[0][1]) / 2;
				}
			}

			//settlement
			
			if (s == 1 && sex == 0) { //for a male
				dinfo.S0 = (dinfo.parent->S0[s][0] + dinfo.parent->S0[s][1]) / 2; //take the average of those two values
				dinfo.alphaS = (dinfo.parent->alphaS[s][0] + dinfo.parent->alphaS[s][1]) / 2;
				dinfo.betaS = (dinfo.parent->betaS[s][0] + dinfo.parent->betaS[s][1]) / 2;
				if (mod->size_or_time == 0) { //size dependent
					dinfo.comp_size = (dinfo.parent->comp_size[s][0] + dinfo.parent->comp_size[s][1]) / 2;
				}
				else {
					dinfo.comp_time = (dinfo.parent->comp_time[s][0] + dinfo.parent->comp_time[s][1]) / 2;

				}
			}
			else {
				//calculate mean
				dinfo.S0 = (dinfo.parent->S0[0][0] + dinfo.parent->S0[0][1]) / 2; //take the average of those two values
				dinfo.alphaS = (dinfo.parent->alphaS[0][0] + dinfo.parent->alphaS[0][1]) / 2;
				dinfo.betaS = (dinfo.parent->betaS[0][0] + dinfo.parent->betaS[0][1]) / 2;
				if (mod->size_or_time == 0) { //size dependent
					dinfo.comp_size = (dinfo.parent->comp_size[0][0] + dinfo.parent->comp_size[0][1]) / 2;
				}
				else {
					dinfo.comp_time = (dinfo.parent->comp_time[0][0] + dinfo.parent->comp_time[0][1]) / 2;

				}
			}
			
		}
	}
	else { //asexual, female only, 
		if (pstage->s_emig.densdep == 0) {
			dinfo.emig_prob = dinfo.parent->emig_prob[0][0];
		}
		else {
			dinfo.D0 = dinfo.parent->D0[0][0];
			dinfo.alpha = dinfo.parent->alpha[0][0];
			dinfo.beta = dinfo.parent->beta[0][0];
		}

		//transfer
		if (mod->size_or_time == 0) { //size dependent, so need the slope
			dinfo.SL_m = dinfo.parent->SL_m[0][0];
			dinfo.rho_m = dinfo.parent->rho_m[0][0];
			dinfo.min_active_size = dinfo.parent->min_active_size[0][0]; 
			dinfo.grow_info.m = dinfo.parent->m[0][0];
			dinfo.grow_info.Linf = dinfo.parent->Linf[0][0];
		}
		else {
			dinfo.SL = dinfo.parent->SL[0][0];
			dinfo.active_SL = dinfo.parent->active_SL[0][0];
			dinfo.rho = dinfo.parent->rho[0][0];
			dinfo.min_active_time = dinfo.parent->min_active_time[0][0];
		}

		//settlement
		dinfo.S0 = dinfo.parent->S0[0][0];
		dinfo.alphaS = dinfo.parent->alphaS[0][0];
		dinfo.betaS = dinfo.parent->betaS[0][0];
		if (mod->size_or_time == 0) { //size dependent
			dinfo.comp_size = dinfo.parent->comp_size[0][0];
		}
		else {
			dinfo.comp_time = dinfo.parent->comp_time[0][0];

		}
	}

	//apply any evolution
	//if (dinfo.ep_evo == 1) { //if emigration rate evolves
	//	uniform_real_distribution<float> mut_dist;
	//	float evo_prob = unif_dist(eng);
	//	if (evo_prob < dinfo.ep_mut) { //if random probability is less than mutation rate, then mutation occurs
	//		if (dinfo.D0 != -9) { //if emigration is density dependent
	//			//do D0
	//			mut_dist = uniform_real_distribution<float>(-dinfo.ep_mut_size, dinfo.ep_mut_size);
	//			dinfo.D0 += mut_dist(eng);
	//			//first do alpha
	//			mut_dist = uniform_real_distribution<float>(-dinfo.alpha_mut_size, dinfo.alpha_mut_size);
	//			dinfo.alpha += mut_dist(eng);
	//			//then do beta
	//			mut_dist = uniform_real_distribution<float>(-dinfo.beta_mut_size, dinfo.beta_mut_size);
	//			dinfo.beta += mut_dist(eng);
	//		}
	//		else {
	//			mut_dist = uniform_real_distribution<float>(-dinfo.ep_mut_size, dinfo.ep_mut_size);
	//			dinfo.emig_prob += mut_dist(eng);
	//		}
	//	}
	//}

	//if (dinfo.active_evo == 1) { //if transfer rate evolves
	//	uniform_real_distribution<float> mut_dist;
	//	float evo_prob = unif_dist(eng);
	//	if (evo_prob < dinfo.S0_mut) { //if random probability is less than mutation rate, then mutation occurs
	//		mut_dist = uniform_real_distribution<float>(-dinfo.active_mut_size, dinfo.active_mut_size);
	//		if (mod->size_or_time == 1) { //time dependent 
	//			dinfo.min_active_time += mut_dist(eng);
	//		}
	//		else { //size_dependent
	//			dinfo.min_active_size += mut_dist(eng);
	//		}
	//	}
	//}

	//if (dinfo.S0_evo == 1) { //if settlement rate evolves
	//	uniform_real_distribution<float> mut_dist;
	//	float evo_prob = unif_dist(eng);
	//	if (evo_prob < dinfo.S0_mut) { //if random probability is less than mutation rate, then mutation occurs
	//		if (dinfo.alphaS != -9) { //if emigration is density dependent
	//			//do S0
	//			mut_dist = uniform_real_distribution<float>(-dinfo.S0_mut_size, dinfo.S0_mut_size);
	//			dinfo.S0 += mut_dist(eng);
	//			//first do alphaS
	//			mut_dist = uniform_real_distribution<float>(-dinfo.alphaS_mut_size, dinfo.alphaS_mut_size);
	//			dinfo.alphaS += mut_dist(eng);
	//			//then do betaS
	//			mut_dist = uniform_real_distribution<float>(-dinfo.betaS_mut_size, dinfo.betaS_mut_size);
	//			dinfo.betaS += mut_dist(eng);
	//		}
	//		else {
	//			mut_dist = uniform_real_distribution<float>(-dinfo.S0_mut_size, dinfo.S0_mut_size);
	//			dinfo.S0 += mut_dist(eng);
	//		}
	//	}
	//}
	//if (dinfo.comp_evo == 1) { //if competency time/size evolves
	//	uniform_real_distribution<float> mut_dist(-dinfo.comp_mut_size, dinfo.comp_mut_size);
	//	float evo_prob = unif_dist(eng);
	//	if (mod->size_or_time == 0) {//size
	//		if (evo_prob < dinfo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
	//			dinfo.comp_size += mut_dist(eng);
	//		}
	//	}
	//	else {//time
	//		if (evo_prob < dinfo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
	//			dinfo.comp_time += mut_dist(eng);
	//		}
	//	}
	//}

}

void Individual:: parameter_iiv() {

}
//void Individual::inheritance(Individual* parent) { //this one is used for asexual reproduction, when only genetic material from the mother is inherited
//	//give new individual the starting coordinate from parent's settlement spot
//	
//	current_pos.x = parent->current_pos.x; current_pos.y = parent->current_pos.y; current_pos.z = parent->current_pos.z;
//	m_x_movements.push_back(current_pos.x); m_y_movements.push_back(current_pos.y); m_z_movements.push_back(current_pos.z);
//	memory.push(current_pos);
//
//	//give it the current patch number from the parent
//	if (mod->patchmodel==1) { current_patch = parent->current_patch; natal_patch = parent->current_patch->patch_ID; dinfo.tested_sett.push_back(natal_patch); } //if the patch pointer points to anything, give that pointer to individual too
//	current_cell = parent->current_cell;
//	natal_cell = current_cell->cell_ID;
//	if (mod->patchmodel == 0) { dinfo.tested_sett.push_back(natal_cell); } //only if it's a cell based model will settlement be only looking at individual cells
//	
//	if (mod->size_or_time == 0) {
//		size_at_birth = parent->size_at_birth;
//		size = parent->size_at_birth;
//	}
//
//
//	//I check for stagedependency because:
//	//if a phase is stagedependent, then no inter individual variability or evolution is allowed, so asexual inheritance matters very little and new dispersal information can be updated without issues
//		//when individual is initialised (before inheritance) it will already have the base values for its stage (if stage structured), so this is just about which values to update
//	//if a phase is NOT stage dependent, iiv and evo are possible but ONLY ONE STAGE CAN EMIGRATE, so individual can inherit varied values from parent
//
//	if (parent->pstage->s_emig.stagedep == 0) { //iiv and evo are possible but only one stage can emigrate, so can just take parent's value
//		//emigration parameters
//		dinfo.emig_prob = parent->dinfo.emig_prob;
//		dinfo.D0 = parent->dinfo.D0;
//		dinfo.alpha = parent->dinfo.alpha;
//		dinfo.beta = parent->dinfo.beta;
//		if (dinfo.ep_evo == 1) { //if emigration rate evolves
//			uniform_real_distribution<float> mut_dist;
//			float evo_prob = unif_dist(eng);
//			if (evo_prob < dinfo.ep_mut) { //if random probability is less than mutation rate, then mutation occurs
//				if (dinfo.D0 != -9) { //if emigration is density dependent
//					//do D0
//					mut_dist = uniform_real_distribution<float>(-dinfo.ep_mut_size, dinfo.ep_mut_size);
//					dinfo.D0 += mut_dist(eng);
//					//first do alpha
//					mut_dist = uniform_real_distribution<float>(-dinfo.alpha_mut_size, dinfo.alpha_mut_size);
//					dinfo.alpha += mut_dist(eng);
//					//then do beta
//					mut_dist = uniform_real_distribution<float>(-dinfo.beta_mut_size, dinfo.beta_mut_size);
//					dinfo.beta += mut_dist(eng);
//				}
//				else {
//					mut_dist = uniform_real_distribution<float>(-dinfo.ep_mut_size, dinfo.ep_mut_size);
//					dinfo.emig_prob += mut_dist(eng);
//				}
//			}
//		}
//	}
//	//transfer parameters
//	if (parent->pstage->s_trans.stagedep == 0) {
//		if (mod->size_or_time == 0) { //size dependent, so need the slope
//			dinfo.SL_m = parent->dinfo.SL_m;
//			dinfo.rho_m = parent->dinfo.rho_m;
//		}
//		else {
//			dinfo.SL = parent->dinfo.SL;
//			dinfo.active_SL = parent->dinfo.active_SL;
//			dinfo.rho = parent->dinfo.rho;
//		}
//		dinfo.min_active_size = parent->dinfo.min_active_size;
//		dinfo.min_active_time = parent->dinfo.min_active_time;
//		dinfo.grow_info.m = parent->dinfo.grow_info.m; //linear method
//		dinfo.grow_info.Linf = parent->dinfo.grow_info.Linf; //gompertz 
//		if (dinfo.active_evo == 1) { //if emigration rate evolves
//			uniform_real_distribution<float> mut_dist;
//			float evo_prob = unif_dist(eng);
//			if (evo_prob < dinfo.S0_mut) { //if random probability is less than mutation rate, then mutation occurs
//				mut_dist = uniform_real_distribution<float>(-dinfo.active_mut_size, dinfo.active_mut_size);
//				if (mod->size_or_time == 1) { //time dependent 
//					dinfo.min_active_time += mut_dist(eng);
//				}
//				else { //size_dependent
//					dinfo.min_active_size += mut_dist(eng);
//				}
//			}
//		}
//
//	}
//	
//	//settlement parameters
//	if (parent->pstage->s_sett.stagedep == 0) { 
//		dinfo.S0 = parent->dinfo.S0; //settlement prob
//		dinfo.alphaS = parent->dinfo.alphaS; //for density dependent sett
//		dinfo.betaS = parent->dinfo.betaS;
//		dinfo.comp_size = parent->dinfo.comp_size; //minimum size for settlement
//		dinfo.comp_time = parent->dinfo.comp_time; //minimum time before settlement
//
//		if (dinfo.S0_evo == 1) { //if emigration rate evolves
//			uniform_real_distribution<float> mut_dist;
//			float evo_prob = unif_dist(eng);
//			if (evo_prob < dinfo.S0_mut) { //if random probability is less than mutation rate, then mutation occurs
//				if (dinfo.alphaS != -9) { //if emigration is density dependent
//					//do S0
//					mut_dist = uniform_real_distribution<float>(-dinfo.S0_mut_size, dinfo.S0_mut_size);
//					dinfo.S0 += mut_dist(eng);
//					//first do alphaS
//					mut_dist = uniform_real_distribution<float>(-dinfo.alphaS_mut_size, dinfo.alphaS_mut_size);
//					dinfo.alphaS += mut_dist(eng);
//					//then do betaS
//					mut_dist = uniform_real_distribution<float>(-dinfo.betaS_mut_size, dinfo.betaS_mut_size);
//					dinfo.betaS += mut_dist(eng);
//				}
//				else {
//					mut_dist = uniform_real_distribution<float>(-dinfo.S0_mut_size, dinfo.S0_mut_size);
//					dinfo.S0 += mut_dist(eng);
//				}
//			}
//		}
//		if (dinfo.comp_evo == 1) { //if competency time/size evolves
//			uniform_real_distribution<float> mut_dist(-dinfo.comp_mut_size, dinfo.comp_mut_size);
//			float evo_prob = unif_dist(eng);
//			if (mod->size_or_time == 0) {//size
//				if (evo_prob < dinfo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
//					dinfo.comp_size += mut_dist(eng);
//				}
//			}
//			else {//time
//				if (evo_prob < dinfo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
//					dinfo.comp_time += mut_dist(eng);
//				}
//			}
//		}
//	}
//
//}

void Individual::inheritance(Individual* mother, Individual* father) {
	//give new individual the starting coordinate from mother's settlement spot
	vector<float>::iterator it_x = mother->m_x_movements.end() - 1; //end of vector storing movement data
	vector<float>::iterator it_y = mother->m_y_movements.end() - 1;
	vector<float>::iterator it_z = mother->m_z_movements.end() - 1;
	m_x_movements.push_back(*it_x);
	m_y_movements.push_back(*it_y);
	m_z_movements.push_back(*it_z);
	current_pos.x = *it_x; current_pos.y = *it_y; current_pos.z= * it_z;
	memory.push(current_pos);

	//give it the current patch number from the parent
	if (mother->current_patch) { current_patch = mother->current_patch; natal_patch = mother->current_patch->patch_ID; }//if the patch pointer points to anything, give that pointer to individual too
	current_cell = mother->current_cell;
	natal_cell = current_cell->cell_ID;
	
	if (mod->patchmodel == 1) {
		dinfo.tested_sett.push_back(natal_patch);
	}
	else {
		dinfo.tested_sett.push_back(natal_cell);
	}
	if (current_cell == nullptr) { cout << "didn't assign current cell properly in inheritance" << endl; }

	//filling parental info
	if (father == nullptr) { //so asexual, inheritance only from mother
			//initial size, because this can't be sexdep from input, i only need female alleles but from both parents
		if (mod->size_or_time == 0) { //if disp is size dependent
			dinfo.parent->size_at_birth[0][0] = mother->dinfo.parent->size_at_birth[0][0];
			size_at_birth = dinfo.parent->size_at_birth[0][0]; 
			size = size_at_birth;
		}
		if (pstage->s_emig.stagedep == 1) { return; } //since there is no iiv or evo, taking values from population is fine
		//emigration
		if (pstage->s_emig.densdep == 0) {
			//fill parent details
			dinfo.parent->emig_prob[0][0] = mother->dinfo.parent->emig_prob[0][0];
			//apply any evolution
			if (pop_dyn->pop_evo.ep_evo == 1) { //if emigration rate evolves
				uniform_real_distribution<float> mut_dist;
				float evo_prob = unif_dist(eng);
				if (evo_prob < pop_dyn->pop_evo.ep_mut) { //if random probability is less than mutation rate, then mutation occurs
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.ep_mut_size, pop_dyn->pop_evo.ep_mut_size);
					if ((dinfo.parent->emig_prob[0][0] + mut_dist(eng)) < 0) { dinfo.parent->emig_prob[0][0] = 0; }
					else {
						dinfo.parent->emig_prob[0][0] += mut_dist(eng);
					}
					
					
				}
			}
		}
		else {
			//fill parent details
			dinfo.parent->D0[0][0] = mother->dinfo.parent->D0[0][0];
			dinfo.parent->alpha[0][0] = mother->dinfo.parent->alpha[0][0];
			dinfo.parent->beta[0][0] = mother->dinfo.parent->beta[0][0];
			//apply any evolution
			if (pop_dyn->pop_evo.ep_evo == 1) { //if emigration rate evolves
				uniform_real_distribution<float> mut_dist;
				float evo_prob = unif_dist(eng);
				if (evo_prob < pop_dyn->pop_evo.ep_mut) { //if random probability is less than mutation rate, then mutation occurs
					//do D0
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.ep_mut_size,pop_dyn->pop_evo.ep_mut_size);
					if ((dinfo.parent->D0[0][0] + mut_dist(eng)) < 0) { dinfo.parent->D0[0][0] = 0; }
					else { dinfo.parent->D0[0][0] += mut_dist(eng); }
					//first do alpha
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.alpha_mut_size, pop_dyn->pop_evo.alpha_mut_size);
					if ((dinfo.parent->alpha[0][0] + mut_dist(eng)) < 0) { dinfo.parent->alpha[0][0] = 0; }
					else { dinfo.parent->alpha[0][0] += mut_dist(eng); }
					//then do beta
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.beta_mut_size, pop_dyn->pop_evo.beta_mut_size);
					if ((dinfo.parent->beta[0][0] + mut_dist(eng)) < 0) { dinfo.parent->beta[0][0] = 0; }
					else { dinfo.parent->beta[0][0] += mut_dist(eng); }
				}
			}
			
		}

		//transfer
		if (mod->size_or_time == 0) { //size dependent, so need the slope
			dinfo.parent->SL_m[0][0] = mother->dinfo.parent->SL_m[0][0];
			dinfo.parent->rho_m[0][0] = mother->dinfo.parent->rho_m[0][0];
			dinfo.parent->min_active_size[0][0] = mother->dinfo.parent->min_active_size[0][0];
			dinfo.parent->min_dv_size[0][0] = mother->dinfo.parent->min_dv_size[0][0];
			dinfo.parent->m[0][0] = mother->dinfo.parent->m[0][0]; //linear method
			dinfo.parent->Linf[0][0] = mother->dinfo.parent->Linf[0][0]; //gompertz
			dinfo.parent->G_K[0][0] = mother->dinfo.parent->G_K[0][0];
			dinfo.parent->Ti[0][0] = mother->dinfo.parent->Ti[0][0];
		}
		else {
			dinfo.parent->SL[0][0] = mother->dinfo.parent->SL[0][0];
			dinfo.parent->active_SL[0][0] = mother->dinfo.parent->active_SL[0][0];
			dinfo.parent->rho[0][0] = mother->dinfo.parent->rho[0][0];
			dinfo.parent->min_active_time[0][0] = mother->dinfo.parent->min_active_time[0][0];
			dinfo.parent->min_active_time[0][0] = mother->dinfo.parent->min_dv_time[0][0];

		}
		if (pop_dyn->pop_evo.active_evo == 1) { //if transfer rate evolves
			uniform_real_distribution<float> mut_dist;
			float evo_prob = unif_dist(eng);
			if (evo_prob < pop_dyn->pop_evo.active_mut) { //if random probability is less than mutation rate, then mutation occurs
				mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.active_mut_size, pop_dyn->pop_evo.active_mut_size);
				if (mod->size_or_time == 1) { //time dependent 
					if ((dinfo.parent->min_active_time[0][0] + mut_dist(eng)) < 0) { dinfo.parent->min_active_time[0][0] = 0; }
					else { dinfo.parent->min_active_time[0][0] += mut_dist(eng); }
				}
				else { //size_dependent
					if ((dinfo.parent->min_active_size[0][0] + mut_dist(eng)) < 0) { dinfo.parent->min_active_size[0][0] = 0; }
					else { dinfo.parent->min_active_size[0][0] += mut_dist(eng); }
				}
			}
		}
		if (pop_dyn->pop_evo.dv_evo == 1) { //if min dv size/time evolves
			uniform_real_distribution<float> mut_dist;
			float evo_prob = unif_dist(eng);
			if (evo_prob < pop_dyn->pop_evo.dv_mut) { //if random probability is less than mutation rate, then mutation occurs
				mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.dv_mut_size, pop_dyn->pop_evo.dv_mut_size);
				if (mod->size_or_time == 1) { //time dependent 
					if ((dinfo.parent->min_dv_time[0][0] + mut_dist(eng)) < 0) { dinfo.parent->min_dv_time[0][0] = 0; }
					else { dinfo.parent->min_dv_time[0][0] += mut_dist(eng); }
				}
				else { //size_dependent
					if ((dinfo.parent->min_dv_size[0][0] + mut_dist(eng)) < 0) { dinfo.parent->min_dv_size[0][0] = 0; }
					else { dinfo.parent->min_dv_size[0][0] += mut_dist(eng); }
				}
			}
		}
		if (pop_dyn->pop_evo.growth_evo == 1) { //if growth evolves
			uniform_real_distribution<float> mut_dist;
			float evo_prob = unif_dist(eng);
			if (evo_prob < pop_dyn->pop_evo.growth_mut) { //if random probability is less than mutation rate, then mutation occurs
				if (mod->grow_method == "linear") { //linear
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.m_mut_size, pop_dyn->pop_evo.m_mut_size);
					if ((dinfo.parent->m[0][0] + mut_dist(eng)) < 0) { dinfo.parent->m[0][0] = 0; }
					else { dinfo.parent->m[0][0] += mut_dist(eng); }
				}
				else if (mod->grow_method == "gompertz") {
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.Linf_mut_size, pop_dyn->pop_evo.Linf_mut_size);
					if ((dinfo.parent->Linf[0][0] + mut_dist(eng)) < 0) { dinfo.parent->Linf[0][0] = 0; }
					else { dinfo.parent->Linf[0][0] += mut_dist(eng); }
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.GK_mut_size, pop_dyn->pop_evo.GK_mut_size);
					if ((dinfo.parent->G_K[0][0] + mut_dist(eng)) < 0) { dinfo.parent->G_K[0][0] = 0; }
					else { dinfo.parent->G_K[0][0] += mut_dist(eng); }
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.Ti_mut_size, pop_dyn->pop_evo.Ti_mut_size);
					if ((dinfo.parent->Ti[0][0] + mut_dist(eng)) < 0) { dinfo.parent->Ti[0][0] = 0; }
					else { dinfo.parent->Ti[0][0] += mut_dist(eng); }
				}
				else {

				}

			}
		}

		//settlement
		dinfo.parent->S0[0][0] = mother->dinfo.parent->S0[0][0];
		dinfo.parent->alphaS[0][0] = mother->dinfo.parent->alphaS[0][0];
		dinfo.parent->betaS[0][0] = mother->dinfo.parent->betaS[0][0];
		if (mod->size_or_time == 0) { //size dependent
			dinfo.parent->comp_size[0][0] = mother->dinfo.parent->comp_size[0][0];
		}
		else {
			dinfo.parent->comp_time[0][0] = mother->dinfo.parent->comp_time[0][0];
		}
		if (pop_dyn->pop_evo.S0_evo == 1) { //if settlement rate evolves
			uniform_real_distribution<float> mut_dist;
			float evo_prob = unif_dist(eng);
			if (evo_prob < pop_dyn->pop_evo.S0_mut) { //if random probability is less than mutation rate, then mutation occurs
				//do S0
				mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.S0_mut_size, pop_dyn->pop_evo.S0_mut_size);
				if ((dinfo.parent->S0[0][0] + mut_dist(eng)) < 0) { dinfo.parent->S0[0][0] = 0; }
				else { dinfo.parent->S0[0][0] += mut_dist(eng); }
				if (pstage->s_sett.densdep == 1) { //if emigration is density dependent
					//first do alphaS
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.alphaS_mut_size, pop_dyn->pop_evo.alphaS_mut_size);
					if ((dinfo.parent->alphaS[0][0] + mut_dist(eng)) < 0) { dinfo.parent->alphaS[0][0] = 0; }
					else { dinfo.parent->alphaS[0][0] += mut_dist(eng); }
					//then do betaS
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.betaS_mut_size, pop_dyn->pop_evo.betaS_mut_size);
					if ((dinfo.parent->betaS[0][0] + mut_dist(eng)) < 0) { dinfo.parent->betaS[0][0] = 0; }
					else { dinfo.parent->betaS[0][0] += mut_dist(eng); }
				}
			}
		}
		if (pop_dyn->pop_evo.comp_evo == 1) { //if competency time/size evolves
			uniform_real_distribution<float> mut_dist(-pop_dyn->pop_evo.comp_mut_size, pop_dyn->pop_evo.comp_mut_size);
			float evo_prob = unif_dist(eng);
			if (mod->size_or_time == 0) {//size
				if (evo_prob < pop_dyn->pop_evo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
					if ((dinfo.parent->comp_size[0][0] + mut_dist(eng)) < 0) { dinfo.parent->comp_size[0][0] = 0; }
					else { dinfo.parent->comp_size[0][0] += mut_dist(eng); }
				}
			}
			else {//time
				if (evo_prob < pop_dyn->pop_evo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
					if ((dinfo.parent->comp_time[0][0] + mut_dist(eng)) < 0) { dinfo.parent->comp_time[0][0] = 0; }
					else { dinfo.parent->comp_time[0][0] += mut_dist(eng); }
				}
			}
		}

		if (pstage->s_emig.stagedep == 0 && pstage->s_emig.emigstage == 0) { //if emigration stage is stage 0, so there is no stage development before indiv needs its dispersal information, call gene expression
			disp_gene_expression();
		}
		else if (sex == 0 && pstage->s_emig.f_emig_prob != 9 || pstage->s_emig.f_D0 != -9) { //if female emig prob is not -9, therefore it should emigrate now
			disp_gene_expression();
		}
		else if (sex == 1 && pstage->s_emig.m_emig_prob != -9 || pstage->s_emig.m_D0 != -9) {
			disp_gene_expression();
		}
		return;
	}

	//I check for stagedependency because:
	//if a phase is stagedependent, then no inter individual variability or evolution is allowed, so inheritance matters very little and new dispersal information can be updated without issues
	//when individual is initialised (before inheritance) it will already have the base values for its stage (if stage structured), so this is just about which values to update
	//if a phase is NOT stage dependent, iiv and evo are possible but ONLY ONE STAGE CAN EMIGRATE, so individual can inherit varied values from parent
	bernoulli_distribution pd(0.5);	
	int n_sexes;
	if (pstage->s_emig.sexdep == 1) { n_sexes = 2; } //need separate values for both sexes
	else { n_sexes = 1; } //same value for both sexes, so no need to keep it separate

	//initial size, because this can't be sexdep from input, i only need female alleles but from both parents
	if (mod->size_or_time == 0) { //if disp is size dependent
		dinfo.parent->size_at_birth[0][0] = mother->dinfo.parent->size_at_birth[0][pd(eng)];
		dinfo.parent->size_at_birth[0][1] = father->dinfo.parent->size_at_birth[0][pd(eng)];

		size_at_birth = (dinfo.parent->size_at_birth[0][0] + dinfo.parent->size_at_birth[0][1]) / 2; //take the average of those two values
		size = size_at_birth;
		if (size == -9) { cout << "size at birth isn't set right in inheritance" << endl; }
	}

	if (pstage->s_emig.stagedep == 1) { return; } //since there is no iiv or evo, taking values from population is fine
	else {

		//here it is also important whether it is sex dependent
		//if sex INDEPENDENT, then there is 1 locus with 2 alleles. the offspring expresses the AVERAGE of the parents' values
		//if sex dependent then there are 2 loci with 2 alleles and indiv randomy inherits one from each parent and expresses the average
		//remember the full matrix looks like this:
		//			Mother		Father
		// female	
		// male

		for (int s = 0; s < n_sexes; s++) {
			//emigration
			if (pstage->s_emig.densdep == 0) {
				//fill parent details
				dinfo.parent->emig_prob[s][0] = mother->dinfo.parent->emig_prob[s][pd(eng)];
				dinfo.parent->emig_prob[s][1] = father->dinfo.parent->emig_prob[s][pd(eng)];
				if (pop_dyn->pop_evo.ep_evo == 1) { //if emigration rate evolves
					uniform_real_distribution<float> mut_dist;
					float evo_prob = unif_dist(eng);
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.ep_mut_size, pop_dyn->pop_evo.ep_mut_size);
					if (evo_prob < pop_dyn->pop_evo.ep_mut) { //if random probability is less than mutation rate, then mutation occurs
						if ((dinfo.parent->emig_prob[s][0] + mut_dist(eng)) < 0) { dinfo.parent->emig_prob[s][0] = 0; }
						else { dinfo.parent->emig_prob[s][0] += mut_dist(eng); }
					}
					evo_prob = unif_dist(eng); //draw another random number 
					if (evo_prob < pop_dyn->pop_evo.ep_mut) { //if random probability is less than mutation rate, then mutation occurs
						if ((dinfo.parent->emig_prob[s][1] + mut_dist(eng)) < 0) { dinfo.parent->emig_prob[s][1] = 0; }
						else { dinfo.parent->emig_prob[s][1] += mut_dist(eng); }
					}
					//need to draw them separately to allow either or both to mutate if the probability happens

				}
			}
			else {
				//fill parent details
				dinfo.parent->D0[s][0] = mother->dinfo.parent->D0[s][pd(eng)];
				dinfo.parent->D0[s][1] = father->dinfo.parent->D0[s][pd(eng)];
				dinfo.parent->alpha[s][0] = mother->dinfo.parent->alpha[s][pd(eng)];
				dinfo.parent->alpha[s][1] = father->dinfo.parent->alpha[s][pd(eng)];
				dinfo.parent->beta[s][0] = mother->dinfo.parent->beta[s][pd(eng)];
				dinfo.parent->beta[s][1] = father->dinfo.parent->beta[s][pd(eng)];		
				if (pop_dyn->pop_evo.ep_evo == 1) { //if emigration rate evolves
					uniform_real_distribution<float> mut_dist;
					float evo_prob = unif_dist(eng);
					float evo_prob2 = unif_dist(eng);
					if (evo_prob < pop_dyn->pop_evo.ep_mut) { //if random probability is less than mutation rate, then mutation occurs
						//do D0
						mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.ep_mut_size, pop_dyn->pop_evo.ep_mut_size);
						if ((dinfo.parent->D0[s][0] + mut_dist(eng)) < 0) { dinfo.parent->D0[s][0] = 0; }
						else { dinfo.parent->D0[s][0] += mut_dist(eng); }
						//first do alpha
						if (pop_dyn->pop_evo.alpha_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.alpha_mut_size, pop_dyn->pop_evo.alpha_mut_size);
							if ((dinfo.parent->alpha[s][0] + mut_dist(eng)) < 0) { dinfo.parent->alpha[s][0] = 0; }
							else { dinfo.parent->alpha[s][0] += mut_dist(eng); }
						}
						//then do beta
						if (pop_dyn->pop_evo.beta_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.beta_mut_size, pop_dyn->pop_evo.beta_mut_size);
							dinfo.parent->beta[s][0] += mut_dist(eng);
							if (dinfo.parent->beta[s][0] < 0) { dinfo.parent->beta[s][0] = 0; }
						}
					}
					if (evo_prob2 < pop_dyn->pop_evo.ep_mut) { //if random probability is less than mutation rate, then mutation occurs
						//do D0
						mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.ep_mut_size, pop_dyn->pop_evo.ep_mut_size);
						dinfo.parent->D0[s][1] += mut_dist(eng);
						if (dinfo.parent->D0[s][1] < 0) { dinfo.parent->D0[s][1] = 0; }
						//first do alpha
						if (pop_dyn->pop_evo.alpha_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.alpha_mut_size, pop_dyn->pop_evo.alpha_mut_size);
							dinfo.parent->alpha[s][1] += mut_dist(eng);
							if (dinfo.parent->alpha[s][1] < 0) { dinfo.parent->alpha[s][1] = 0; }
						}
						//then do beta
						if (pop_dyn->pop_evo.beta_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.beta_mut_size, pop_dyn->pop_evo.beta_mut_size);
							dinfo.parent->beta[s][1] += mut_dist(eng);
							if (dinfo.parent->beta[s][1] < 0) { dinfo.parent->beta[s][1] = 0; }
						}
					}
				}
			}

			//transfer
			if (mod->size_or_time == 0) { //size dependent, so need the slope
				dinfo.parent->SL_m[s][0] = mother->dinfo.parent->SL_m[s][pd(eng)];
				dinfo.parent->SL_m[s][1] = father->dinfo.parent->SL_m[s][pd(eng)];
				dinfo.parent->rho_m[s][0] = mother->dinfo.parent->rho_m[s][pd(eng)];
				dinfo.parent->rho_m[s][1] = father->dinfo.parent->rho_m[s][pd(eng)];
				dinfo.parent->min_active_size[s][0] = mother->dinfo.parent->min_active_size[s][pd(eng)];
				dinfo.parent->min_active_size[s][1] = father->dinfo.parent->min_active_size[s][pd(eng)];
				dinfo.parent->min_dv_size[s][0] = mother->dinfo.parent->min_dv_size[s][pd(eng)];
				dinfo.parent->min_dv_size[s][1] = father->dinfo.parent->min_dv_size[s][pd(eng)];
				dinfo.parent->m[s][0] = mother->dinfo.parent->m[s][pd(eng)]; //linear method
				dinfo.parent->m[s][1] = father->dinfo.parent->m[s][pd(eng)]; //linear method
				dinfo.parent->Linf[s][0] = mother->dinfo.parent->Linf[s][pd(eng)]; //gompertz
				dinfo.parent->Linf[s][1] = father->dinfo.parent->Linf[s][pd(eng)]; //gompertz
				//if (dinfo.min_active_size == -9) { cout << "min active size not set right at inheritance" << endl; }
			}
			else {
				dinfo.parent->SL[s][0] = mother->dinfo.parent->SL[s][pd(eng)];
				dinfo.parent->SL[s][1] = father->dinfo.parent->SL[s][pd(eng)];
				dinfo.parent->active_SL[s][0] = mother->dinfo.parent->active_SL[s][pd(eng)];
				dinfo.parent->active_SL[s][1] = father->dinfo.parent->active_SL[s][pd(eng)];
				dinfo.parent->rho[s][0] = mother->dinfo.parent->rho[s][pd(eng)];
				dinfo.parent->rho[s][1] = father->dinfo.parent->rho[s][pd(eng)];
				dinfo.parent->min_active_time[s][0] = mother->dinfo.parent->min_active_time[s][pd(eng)];
				dinfo.parent->min_active_time[s][1] = father->dinfo.parent->min_active_time[s][pd(eng)];
				dinfo.parent->min_dv_time[s][0] = mother->dinfo.parent->min_dv_time[s][pd(eng)];
				dinfo.parent->min_dv_time[s][1] = father->dinfo.parent->min_dv_time[s][pd(eng)];

			}
			if (pop_dyn->pop_evo.active_evo == 1) { //if transfer rate evolves
				uniform_real_distribution<float> mut_dist;
				float evo_prob = unif_dist(eng);
				float evo_prob2 = unif_dist(eng);
				if (evo_prob < pop_dyn->pop_evo.active_mut) { //if random probability is less than mutation rate, then mutation occurs
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.active_mut_size, pop_dyn->pop_evo.active_mut_size);
					if (mod->size_or_time == 1) { //time dependent 
						dinfo.parent->min_active_time[s][0] += mut_dist(eng);
						if (dinfo.parent->min_active_time[s][0] < 0) { dinfo.parent->min_active_time[s][0] = 0; }
					}
					else { //size_dependent
						dinfo.parent->min_active_size[s][0] += mut_dist(eng);
						if (dinfo.parent->min_active_size[s][0] < 0) { dinfo.parent->min_active_size[s][0] = 0; }
					}
				}
				if (evo_prob2 < pop_dyn->pop_evo.active_mut) { //if random probability is less than mutation rate, then mutation occurs
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.active_mut_size, pop_dyn->pop_evo.active_mut_size);
					if (mod->size_or_time == 1) { //time dependent 
						dinfo.parent->min_active_time[s][1] += mut_dist(eng);
						if (dinfo.parent->min_active_time[s][1] < 0) { dinfo.parent->min_active_time[s][1] = 0; }
					}
					else { //size_dependent
						dinfo.parent->min_active_size[s][1] += mut_dist(eng);
						if (dinfo.parent->min_active_size[s][1] < 0) { dinfo.parent->min_active_size[s][1] = 0; }
					}
				}
			}
			if (pop_dyn->pop_evo.dv_evo == 1) { //if transfer rate evolves
				uniform_real_distribution<float> mut_dist;
				float evo_prob = unif_dist(eng);
				float evo_prob2 = unif_dist(eng);
				if (evo_prob < pop_dyn->pop_evo.dv_mut) { //if random probability is less than mutation rate, then mutation occurs
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.dv_mut_size, pop_dyn->pop_evo.dv_mut_size);
					if (mod->size_or_time == 1) { //time dependent 
						dinfo.parent->min_dv_time[s][0] += mut_dist(eng);
						if (dinfo.parent->min_dv_time[s][0] < 0) { dinfo.parent->min_dv_time[s][0] = 0; }
					}
					else { //size_dependent
						dinfo.parent->min_dv_size[s][0] += mut_dist(eng);
						if (dinfo.parent->min_dv_size[s][0] < 0) { dinfo.parent->min_dv_size[s][0] = 0; }
					}
				}
				if (evo_prob2 < pop_dyn->pop_evo.active_mut) { //if random probability is less than mutation rate, then mutation occurs
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.dv_mut_size, pop_dyn->pop_evo.dv_mut_size);
					if (mod->size_or_time == 1) { //time dependent 
						dinfo.parent->min_dv_time[s][1] += mut_dist(eng);
						if (dinfo.parent->min_dv_time[s][1] < 0) { dinfo.parent->min_dv_time[s][1] = 0; }
					}
					else { //size_dependent
						dinfo.parent->min_dv_size[s][1] += mut_dist(eng);
						if (dinfo.parent->min_dv_size[s][0] < 0) { dinfo.parent->min_dv_size[s][0] = 0; }
					}
				}
			}
			if (pop_dyn->pop_evo.growth_evo == 1) { //if growth evolves
				uniform_real_distribution<float> mut_dist;
				float evo_prob = unif_dist(eng);
				float evo_prob2 = unif_dist(eng);
				if (evo_prob < pop_dyn->pop_evo.growth_mut) { //if random probability is less than mutation rate, then mutation occurs
					if (mod->grow_method == "linear") { //linear
						mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.m_mut_size, pop_dyn->pop_evo.m_mut_size);
						dinfo.parent->m[s][0] += mut_dist(eng);
						if (dinfo.parent->m[s][0] < 0) { dinfo.parent->m[s][0] = 0; }
					}
					else if (mod->grow_method == "gompertz") {
						if (pop_dyn->pop_evo.Linf_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.Linf_mut_size, pop_dyn->pop_evo.Linf_mut_size);
							dinfo.parent->Linf[s][0] += mut_dist(eng);
							if (dinfo.parent->Linf[s][0] < 0) { dinfo.parent->Linf[s][0] = 0; }
						}
						if (pop_dyn->pop_evo.GK_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.GK_mut_size, pop_dyn->pop_evo.GK_mut_size);
							dinfo.parent->G_K[s][0] += mut_dist(eng);
							if (dinfo.parent->G_K[s][0] < 0) { dinfo.parent->G_K[s][0] = 0; }
						}
						if (pop_dyn->pop_evo.Ti_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.Ti_mut_size, pop_dyn->pop_evo.Ti_mut_size);
							dinfo.parent->Ti[s][0] += mut_dist(eng);
							if (dinfo.parent->Ti[s][0] < 0) { dinfo.parent->Ti[s][0] = 0; }
						}
					}
					else {

					}

				}
				if (evo_prob2 < pop_dyn->pop_evo.growth_mut) { //if random probability is less than mutation rate, then mutation occurs
					if (mod->grow_method == "linear") { //linear
						mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.m_mut_size, pop_dyn->pop_evo.m_mut_size);
						dinfo.parent->m[s][1] += mut_dist(eng);
						if (dinfo.parent->m[s][1] < 0) { dinfo.parent->m[s][1] = 0; }
					}
					else if (mod->grow_method == "gompertz") {
						if (pop_dyn->pop_evo.Linf_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.Linf_mut_size, pop_dyn->pop_evo.Linf_mut_size);
							dinfo.parent->Linf[s][1] += mut_dist(eng);
							if (dinfo.parent->Linf[s][1] < 0) { dinfo.parent->Linf[s][1] = 0; }
						}
						if (pop_dyn->pop_evo.GK_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.GK_mut_size, pop_dyn->pop_evo.GK_mut_size);
							dinfo.parent->G_K[s][1] += mut_dist(eng);
							if (dinfo.parent->G_K[s][1] < 0) { dinfo.parent->G_K[s][1] = 0; }
						}
						if (pop_dyn->pop_evo.Ti_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.Ti_mut_size, pop_dyn->pop_evo.Ti_mut_size);
							dinfo.parent->Ti[s][1] += mut_dist(eng);
							if (dinfo.parent->Ti[s][1] < 0) { dinfo.parent->Ti[s][1] = 0; }
						}
					}
					else {

					}

				}
			}

			//settlement
			dinfo.parent->S0[s][0] = mother->dinfo.parent->S0[s][pd(eng)];
			dinfo.parent->S0[s][1] = father->dinfo.parent->S0[s][pd(eng)];
			dinfo.parent->alphaS[s][0] = mother->dinfo.parent->alphaS[s][pd(eng)];
			dinfo.parent->alphaS[s][1] = father->dinfo.parent->alphaS[s][pd(eng)];
			dinfo.parent->betaS[s][0] = mother->dinfo.parent->betaS[s][pd(eng)];
			dinfo.parent->betaS[s][1] = father->dinfo.parent->betaS[s][pd(eng)];
			if (mod->size_or_time == 0) { //size dependent
				dinfo.parent->comp_size[s][0] = mother->dinfo.parent->comp_size[s][pd(eng)];
				dinfo.parent->comp_size[s][1] = father->dinfo.parent->comp_size[s][pd(eng)];
			}
			else {
				dinfo.parent->comp_time[s][0] = mother->dinfo.parent->comp_time[s][pd(eng)];
				dinfo.parent->comp_time[s][1] = father->dinfo.parent->comp_time[s][pd(eng)];
			}
			if (pop_dyn->pop_evo.S0_evo == 1) { //if settlement rate evolves
				uniform_real_distribution<float> mut_dist;
				float evo_prob = unif_dist(eng);
				float evo_prob2 = unif_dist(eng);
				if (evo_prob < pop_dyn->pop_evo.S0_mut) { //if random probability is less than mutation rate, then mutation occurs
					//do S0
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.S0_mut_size, pop_dyn->pop_evo.S0_mut_size);
					dinfo.parent->S0[s][0] += mut_dist(eng);
					if (dinfo.parent->S0[s][0] < 0) { dinfo.parent->S0[s][0] = 0; }
					if (pstage->s_sett.densdep == 1) { //if emigration is density dependent
						//first do alphaS
						if (pop_dyn->pop_evo.alphaS_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.alphaS_mut_size, pop_dyn->pop_evo.alphaS_mut_size);
							dinfo.parent->alphaS[s][0] += mut_dist(eng);
							if (dinfo.parent->alphaS[s][0] < 0) { dinfo.parent->alphaS[s][0] = 0; }
						}
						//then do betaS
						if (pop_dyn->pop_evo.betaS_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.betaS_mut_size, pop_dyn->pop_evo.betaS_mut_size);
							dinfo.parent->betaS[s][0] += mut_dist(eng);
							if (dinfo.parent->betaS[s][0] < 0) { dinfo.parent->betaS[s][0] = 0; }
						}
					}
				}
				if (evo_prob2 < pop_dyn->pop_evo.S0_mut) { //if random probability is less than mutation rate, then mutation occurs
					//do S0
					mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.S0_mut_size, pop_dyn->pop_evo.S0_mut_size);
					dinfo.parent->S0[s][1] += mut_dist(eng);
					if (dinfo.parent->S0[s][1] < 0) { dinfo.parent->S0[s][1] = 0; }

					if (pstage->s_sett.densdep == 1) { //if emigration is density dependent
						//first do alphaS
						if (pop_dyn->pop_evo.alphaS_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.alphaS_mut_size, pop_dyn->pop_evo.alphaS_mut_size);
							dinfo.parent->alphaS[s][1] += mut_dist(eng);
							if (dinfo.parent->alphaS[s][1] < 0) { dinfo.parent->alphaS[s][1] = 0; }
						}
						//then do betaS
						if (pop_dyn->pop_evo.betaS_mut_size != -9) {
							mut_dist = uniform_real_distribution<float>(-pop_dyn->pop_evo.betaS_mut_size, pop_dyn->pop_evo.betaS_mut_size);
							dinfo.parent->betaS[s][1] += mut_dist(eng);
							if (dinfo.parent->betaS[s][1] < 0) { dinfo.parent->betaS[s][1] = 0; }
						}
					}
				}
			}
			if (pop_dyn->pop_evo.comp_evo == 1) { //if competency time/size evolves
				uniform_real_distribution<float> mut_dist(-pop_dyn->pop_evo.comp_mut_size, pop_dyn->pop_evo.comp_mut_size);
				float evo_prob = unif_dist(eng);
				float evo_prob2 = unif_dist(eng);
				if (mod->size_or_time == 0) {//size
					if (evo_prob < pop_dyn->pop_evo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
						dinfo.parent->comp_size[s][0] += mut_dist(eng);
						if (dinfo.parent->comp_size[s][0] < 0) { dinfo.parent->comp_size[s][0] = 0; }
					}
					if (evo_prob2 < pop_dyn->pop_evo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
						dinfo.parent->comp_size[s][1] += mut_dist(eng);
						if (dinfo.parent->comp_size[s][1] < 0) { dinfo.parent->comp_size[s][1] = 0; }
					}
				}
				else {//time
					if (evo_prob < pop_dyn->pop_evo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
						dinfo.parent->comp_time[s][0] += mut_dist(eng);
						if (dinfo.parent->comp_time[s][0] < 0) { dinfo.parent->comp_time[s][0] = 0; }
					}
					if (evo_prob2 < pop_dyn->pop_evo.comp_mut) { //if random probability is less than mutation rate, then mutation occurs
						dinfo.parent->comp_time[s][1] += mut_dist(eng);
						if (dinfo.parent->comp_time[s][1] < 0) { dinfo.parent->comp_time[s][1] = 0; }
					}
				}
			}
		}

		
	}
	//cout << " for indiv " << ID << ", before disp_expression, S0 is " << dinfo.S0 << ", with alphaS " << dinfo.alphaS << ", with betaS " << dinfo.betaS;
	//cout << " and parent values are " << dinfo.parent->S0[0][0] << ", " << dinfo.parent->S0[0][1] << ", " << dinfo.parent->alphaS[0][0] << ", " << dinfo.parent->alphaS[0][1]
	//	<< ", " << dinfo.parent->betaS[0][0] << ", " << dinfo.parent->betaS[0][1] << endl;
	//dispersal gene expression
	if (pstage->s_emig.stagedep == 0 && pstage->s_emig.emigstage == 0) { //if emigration stage is stage 0, so there is no stage development before indiv needs its dispersal information, call gene expression
		disp_gene_expression();
	}
	else if (sex == 0 && pstage->s_emig.f_emig_prob != 9 || pstage->s_emig.f_D0!=-9) { //if female emig prob is not -9, therefore it should emigrate now
		disp_gene_expression();
	}
	else if (sex == 1 && pstage->s_emig.m_emig_prob != -9 || pstage->s_emig.m_D0 !=-9) {
		disp_gene_expression();
	}
	//cout << " for indiv " << ID << ", after inheritance S0 is " << dinfo.S0 << ", with alphaS " << dinfo.alphaS << ", with betaS " << dinfo.betaS << endl;

}

bool Individual::check_repro() {//check whether individual reproduces
	//default_random_engine unif_generator;
	if (sex == 0) {return false; } //if individual is male, they dont produce offspring (i am assuming that i dont have a harem situation, where male reproductive behaviour matters)
	
	
	//reproduction interval
	if (skip_breeding == true) { //if it needs more time before reproducing (based on Population RepInt)
		rep_int++; //increment rep_int, this individual has skipped this reproductive event
		if (rep_int == pop_dyn->RepInt) { skip_breeding = false; } //if after this reproductive event, the individual is allowed to reproduce, update skip_breeding
		return false;
	}

	//stage
	if (pstage) { //if pstage points to something, it is a stagestructured model. look at sex-specific fecundity rates
		if (pstage->f_fec == 0) { return false; } //if this stage doesn't reproduce yet
	}

	//if it is all clear to reproduce, what is the probability that it will do so?
	/*uniform_real_distribution<float> prob_dist(0, 1);*/
	float rep_prob = unif_dist(eng);
	if (rep_prob < pop_dyn->PRep) { //if the individual reproduces (based on population PRep)
		 return true; //let the individual reproduce
	}
	else {return false; }

}

int Individual::reproduction() {
	//default_random_engine pois_generator;

	//check whether individual reproduces
	if (!check_repro()) { return -9; } //if the individual doesn't reproduce, exit

	//otherwise,
	int n_offspring;

	if (!pstage) { //if it's not stage structured
		float R = pop_dyn->R;
		float max_inds;
		float mean_offspring;
		max_inds = get_localK();
		
		if (mod->repro > 0) { R *= 2; } //if it's a sexual model multiply R by 2, to allow male offspring
		if(mod->patchmodel==1){
			mean_offspring = R / (1 + abs(R - 1)*pow(current_patch->p_subpop->size / max_inds, pop_dyn->bc)); //Maynard Smith and Slatkin model
		}
		else {
			mean_offspring = R / (1 + abs(R - 1)*pow(current_cell->c_subpop->size / max_inds, pop_dyn->bc)); //Maynard Smith and Slatkin model
		}
		poisson_distribution<int> poisson(mean_offspring); //poisson distribution
		n_offspring = poisson(eng);
	}
	else { //stage-structured
		float fec;
		if (pop_dyn->fecdensdep == 1) { //if fecundity is density dependent
			float k = get_localK();
			fec = pstage->f_fec * exp(-1*(1/k) * get_localN());
		}
		else {
			fec = pstage->f_fec; //this way, whether it's a sexual or female-only model, it doesn't matter
		}
		
		poisson_distribution<int> poisson(fec); //poisson distribution
		n_offspring = poisson(eng);
	}

	//rep = true;
	if (pop_dyn->RepInt > 0) { rep_int = 0; skip_breeding = true; } //if they need to skip the next breeding season, reset their interval counter and make skip_breeding true
	

	return n_offspring;
}

//emigration
void Individual::set_emig_info(bool invent) {
	if (pstage->s_emig.sexdep == 1 && sex == 0) { //sex dependent and male
		dinfo.emig_prob = pstage->s_emig.m_emig_prob;
		dinfo.D0 = pstage->s_emig.m_D0;
		dinfo.alpha = pstage->s_emig.m_alpha;
		dinfo.beta = pstage->s_emig.m_beta;
	}
	else {
		dinfo.emig_prob = pstage->s_emig.f_emig_prob;
		dinfo.D0 = pstage->s_emig.f_D0;
		dinfo.alpha = pstage->s_emig.f_alpha;
		dinfo.beta = pstage->s_emig.f_beta;
	}
	//cout << "indiv is stage "<< pstage->stage << ", setting emig prob: " << dinfo.emig_prob << endl;
	/*dinfo.ep_evo = pstage->s_emig.ep_evo;
	dinfo.ep_mut = pstage->s_emig.ep_mut;
	dinfo.ep_mut_size = pstage->s_emig.ep_mut_size;
	dinfo.alpha_mut_size = pstage->s_emig.alpha_mut_size;
	dinfo.beta_mut_size = pstage->s_emig.beta_mut_size;*/
	//if (first == true) {
	//	if (pstage->s_emig.indvar == 1) { //if there is inter-individual variation in emigration parameters
	//		if (pstage->s_emig.densdep == 1) { //if emigration is density dependent
	//			normal_distribution<float> D0_dist(dinfo.D0, pstage->s_emig.iiv_sd);
	//			normal_distribution<float> alpha_dist(dinfo.alpha, pstage->s_emig.iiv_sd);
	//			normal_distribution<float> beta_dist(dinfo.beta, pstage->s_emig.iiv_sd);

	//			dinfo.D0 = D0_dist(eng);
	//			dinfo.alpha = alpha_dist(eng);
	//			dinfo.beta = beta_dist(eng);
	//		}
	//		else {
	//			normal_distribution<float> ep_dist(dinfo.emig_prob, pstage->s_emig.iiv_sd);
	//			dinfo.emig_prob = ep_dist(eng);
	//		}

	//	}

	//}

	if (invent == true && pstage->s_emig.stagedep==0) { //if this individual is part of the first generation of the model, then it doesn't have any parent data so i need to invent some
		//also dispersal needs to be stage independent for parental data to matter. if stage dependent, then can't be iiv or evo so no recombination etc
		stage_infob* temp_si= (*subpop_stages)[pstage->s_emig.emigstage];
		//cout << "checking temp_si, emigrating stage is " << temp_si->stage << endl;
		//cout << "emig stage is " << pstage->s_emig.emigstage << endl;
		float rsamp;
		//emigration parameters
		if (temp_si->s_emig.densdep == 0) {
			dinfo.parent->emig_prob[0][0] = dinfo.parent->emig_prob[0][1]= temp_si->s_emig.f_emig_prob;
		}
		else {
			dinfo.parent->D0[0][0] = dinfo.parent->D0[0][1] = temp_si->s_emig.f_D0;
			dinfo.parent->alpha[0][0] = dinfo.parent->alpha[0][1] = temp_si->s_emig.f_alpha;
			dinfo.parent->beta[0][0] = dinfo.parent->beta[0][1] = temp_si->s_emig.f_beta;
		}
			
		if (temp_si->s_emig.sexdep == 1) { //sex dependent, so need male values in a 2x2 matrix
			if (pstage->s_emig.densdep == 0) {
				dinfo.parent->emig_prob[1][0] = dinfo.parent->emig_prob[1][1] = temp_si->s_emig.m_emig_prob;
			}
			else {
				dinfo.parent->D0[1][0] = dinfo.parent->D0[1][1] = temp_si->s_emig.m_D0;
				dinfo.parent->alpha[1][0] = dinfo.parent->alpha[1][1] = temp_si->s_emig.m_alpha;
				dinfo.parent->beta[1][0] = dinfo.parent->beta[1][1] = temp_si->s_emig.m_beta;
			}
		}
		
		if (temp_si->s_emig.indvar == 1) {
			if (temp_si->s_emig.densdep == 0) {
				normal_distribution<float> ep_dist(temp_si->s_emig.f_emig_prob, pstage->s_emig.iiv_sd);
				rsamp = ep_dist(eng);
				if (rsamp < 0) { rsamp = 0; }
				if (rsamp > 1) { rsamp = 1; }
				dinfo.parent->emig_prob[0][0] = dinfo.parent->emig_prob[0][1] = rsamp;
				if (temp_si->s_emig.sexdep == 1) {
					normal_distribution<float> ep_dist(temp_si->s_emig.m_emig_prob, pstage->s_emig.iiv_sd);
					rsamp = ep_dist(eng);
					if (rsamp < 0) { rsamp = 0; }
					if (rsamp > 1) { rsamp = 1; }
					dinfo.parent->emig_prob[1][0] = dinfo.parent->emig_prob[1][1] = rsamp;
				}				
			}
			else {
				normal_distribution<float> D0_dist(temp_si->s_emig.f_D0, temp_si->s_emig.iiv_sd);
				normal_distribution<float> alpha_dist(temp_si->s_emig.f_alpha, temp_si->s_emig.iiv_sd);
				normal_distribution<float> beta_dist(temp_si->s_emig.f_beta, temp_si->s_emig.iiv_sd);
				rsamp = D0_dist(eng);
				if (rsamp < 0) { rsamp = 0; }
				if (rsamp > 1) { rsamp = 1; }
				dinfo.parent->D0[0][0] = dinfo.parent->D0[0][1] = rsamp;
				rsamp = alpha_dist(eng);
				dinfo.parent->alpha[0][0] = dinfo.parent->alpha[0][1] = rsamp;
				dinfo.parent->beta[0][0] = dinfo.parent->beta[0][1] = rsamp;
				if (temp_si->s_emig.sexdep == 1) {
					normal_distribution<float> D0_dist(temp_si->s_emig.m_D0, temp_si->s_emig.iiv_sd);
					normal_distribution<float> alpha_dist(temp_si->s_emig.m_alpha, temp_si->s_emig.iiv_sd);
					normal_distribution<float> beta_dist(temp_si->s_emig.m_beta, temp_si->s_emig.iiv_sd);
					rsamp = D0_dist(eng);
					if (rsamp < 0) { rsamp = 0; }
					if (rsamp > 1) { rsamp = 1; }
					dinfo.parent->D0[1][0] = dinfo.parent->D0[1][1] = rsamp;
					rsamp = alpha_dist(eng);
					dinfo.parent->alpha[1][0] = dinfo.parent->alpha[1][1] = rsamp;
					rsamp = beta_dist(eng);
					dinfo.parent->beta[1][0] = dinfo.parent->beta[1][1] = rsamp;
				}
			}
		}		
	}
}

//void Individual::set_emig_info() {
//	if (pop_einfo->stagedep == 1) { //if dispersal is stage dependent
//		if (pop_einfo->sexdep == 1 && sex==0) { //and sex dependent, assign male info
//
//			dinfo.emig_prob = pstage->m_ep;
//			dinfo.D0 = pstage->m_D0;
//			dinfo.alpha = pstage->m_alpha;
//			dinfo.beta = pstage->m_beta;
//		}
//		else { //else assign female info
//			dinfo.emig_prob = pstage->f_ep;
//			dinfo.D0 = pstage->f_D0;
//			dinfo.alpha = pstage->f_alpha;
//			dinfo.beta = pstage->f_beta;
//		}
//		dinfo.ep_evo = pstage->ep_evo;
//		dinfo.ep_mut = pstage->ep_mut;
//		dinfo.ep_mut_size = pstage->ep_mut_size;
//		dinfo.alpha_mut_size = pstage->alpha_mut_size;
//		dinfo.beta_mut_size = pstage->beta_mut_size;
//
//	}
//	else { //if emigration is not stage dependent
//		if (pop_einfo->sexdep==1 && sex ==0) { //but it is sex dependent, assign male info
//			dinfo.emig_prob = pop_einfo->m_emig_prob;
//			dinfo.D0 = pop_einfo->m_D0;
//			dinfo.alpha = pop_einfo->m_alpha;
//			dinfo.beta = pop_einfo->m_beta;
//		}
//		else { //else assign female info
//			dinfo.emig_prob = pop_einfo->f_emig_prob;
//			dinfo.D0 = pop_einfo->f_D0;
//			dinfo.alpha = pop_einfo->f_alpha;
//			dinfo.beta = pop_einfo->f_beta;
//
//		}
//		dinfo.ep_evo = pop_einfo->ep_evo;
//		dinfo.ep_mut = pop_einfo->ep_mut;
//		dinfo.ep_mut_size = pop_einfo->ep_mut_size;
//		dinfo.alpha_mut_size = pop_einfo->alpha_mut_size;
//		dinfo.beta_mut_size = pop_einfo->beta_mut_size;
//	}
//	if (pop_einfo->indvar == 1) { //if there is inter-individual variation in emigration parameters
//		if (pop_einfo->densdep == 1) { //if emigration is density dependent
//			normal_distribution<float> D0_dist(dinfo.D0, pop_einfo->iiv_sd);
//			normal_distribution<float> alpha_dist(dinfo.alpha, pop_einfo->iiv_sd);
//			normal_distribution<float> beta_dist(dinfo.beta, pop_einfo->iiv_sd);
//			
//			dinfo.D0 = D0_dist(eng);
//			dinfo.alpha= alpha_dist(eng);
//			dinfo.beta = beta_dist(eng);
//		}
//		else {
//			normal_distribution<float> ep_dist(dinfo.emig_prob, pop_einfo->iiv_sd);
//			dinfo.emig_prob = ep_dist(eng);
//		}
//
//	}
//
//}

bool Individual::check_emig() {
	/*default_random_engine unif_generator;
	uniform_real_distribution<float> emig_dist(0, 1);*/
	float ran_em = unif_dist(eng); //sample randomly from a uniform distribution
	if (sett == true) { return false; } //if the individual has already emigrated and settled, it can't emigrate again. UNLESS I ALLOW RESUSPENSION

	if (pstage->s_emig.densdep == 0) { //if emigration is density independent
		if (ran_em < dinfo.emig_prob) { status = 1; return true; } //if the random sample is less than individual's emigration prob, then it can emigrate, status=disperser
		else { status = 2; return false; } //status=local recruit
	}
	else { //if emigration is density dependent
		float d, n, k;
		if (mod->patchmodel == 1) { n = current_patch->p_subpop->size; }
		else { n = current_cell->c_subpop->size; }

		if (pstage) { //if it's a stage structured model
			d = dinfo.D0 / (1 + exp(-1 * (pop_dyn->bc*n - dinfo.beta)*dinfo.alpha));
		}
		else { //if it's non overlapping
			k = get_localK();
			d = dinfo.D0 / (1 + exp(-1 * (n / k - dinfo.beta)*dinfo.alpha));
		}
		if (ran_em < d) {
			status = 1; return true;
		}
		else { status = 2; return false; }
	}
}

void Individual::calc_SL_rho() { 
	//this function assumes that the model is using size-dependent dispersal, therefore dinfo.SL_m is a MULTIPLIER, how many body lengths
	//and rho_m is the slope of the negative linear relationship of body size to rho
	if (dinfo.mode == 0) { //if the individual never disperses actively, it goes from passive to competent (still passive)
		dinfo.SL = 1;
		dinfo.rho = 0.9;
		return;
	}

	//need to check whether the individual is ready for SL and rho to be updated, ie if it is >minimum size for next stage 
	else if (size < dinfo.min_active_size) { //i dont want it to do anything if it's not big enough to be active yet
		//a lot of the transfer phase conditions are if(dinfo.SL==1) for passive movement, so to update this early will mess things up
		return;
	}
	else {
		//cout << "indiv " << ID << " has SL_m " << dinfo.SL_m << " and size " << size << " mm" ;
		dinfo.SL = dinfo.SL_m*size / 1000; //calculate the step length. /1000 because size is in mm so needs to covnert to m
		dinfo.rho = dinfo.rho_m*size + 1; //the intercept is 1 because the smaller an indiv is, the closer to 1 i am assuming it gets
		//cout << "calculated SL is " << dinfo.SL << " and rho " << dinfo.rho << endl;
	}
}

//transfer
void Individual::set_growth(bool invent) {
	dinfo.grow_info.method = pstage->s_trans.pop_grow.method;
	/*dinfo.grow_info.growth_evo = pstage->s_trans.pop_grow.growth_evo;
	dinfo.grow_info.growth_mut = pstage->s_trans.pop_grow.growth_mut;*/

	if (dinfo.grow_info.method == 0) { //linear
		if (pstage->s_trans.sexdep == 1 && sex == 0) { //and sex dependent, if sex is male then get male info
			dinfo.grow_info.m = pstage->s_trans.pop_grow.m_m;
			dinfo.grow_info.m_sd = pstage->s_trans.pop_grow.m_m_sd;
			dinfo.grow_info.b = pstage->s_trans.pop_grow.m_b;
		}
		else {
			dinfo.grow_info.m = pstage->s_trans.pop_grow.f_m;
			dinfo.grow_info.m_sd = pstage->s_trans.pop_grow.f_m_sd;
			dinfo.grow_info.b = pstage->s_trans.pop_grow.f_b;

		}
		//if (first == true) {
		//	if (dinfo.grow_info.m_sd != -9) { //iiv in the slope of the linear growth function, ie faster/slower growers
		//		normal_distribution<float>m_dist(dinfo.grow_info.m, dinfo.grow_info.m_sd);
		//		dinfo.grow_info.m = m_dist(eng);
		//	}

		//}
		/*dinfo.grow_info.m_mut_size = pstage->s_trans.pop_grow.m_mut_size;
		dinfo.grow_info.b_mut_size = pstage->s_trans.pop_grow.b_mut_size;*/
	}
	else if (dinfo.grow_info.method == 1) { //modified gompertz
		if (pstage->s_trans.sexdep == 1 && sex == 0) { //and sex dependent, if sex is male then get male info
			dinfo.grow_info.Linf = pstage->s_trans.pop_grow.m_Linf; //in mm
			dinfo.grow_info.G_K = pstage->s_trans.pop_grow.m_G_K;
			dinfo.grow_info.Linf_sd = pstage->s_trans.pop_grow.m_Linf_sd;
		}
		else { //female
			dinfo.grow_info.Linf = pstage->s_trans.pop_grow.f_Linf; //in mm
			dinfo.grow_info.G_K = pstage->s_trans.pop_grow.f_G_K;
			dinfo.grow_info.Linf_sd = pstage->s_trans.pop_grow.f_Linf_sd;
		}
		dinfo.grow_info.Ti = pstage->s_trans.pop_grow.Ti; //inflection point, earliest settlement date
		dinfo.grow_info.Ti_sd = pstage->s_trans.pop_grow.Ti_sd;

		
		/*dinfo.grow_info.Linf_mut_size = pstage->s_trans.pop_grow.Linf_mut_size;
		dinfo.grow_info.GK_mut_size = pstage->s_trans.pop_grow.GK_mut_size;*/
	}
	else if (dinfo.grow_info.method == 2) {
		if (pstage->s_trans.sexdep == 1 && sex == 0) { //and sex dependent, if sex is male then get male info

		}
		else { //female

		}

	}

	if (invent == true && pstage->s_trans.pop_grow.stage_dep == 0) { //if this individual is part of the first generation of the model, then it doesn't have any parent data so i need to invent some
		//also dispersal needs to be stage independent for parental data to matter. if stage dependent, then can't be iiv or evo so no recombination etc
		stage_infob* temp_si = (*subpop_stages)[pstage->s_emig.emigstage];
		dinfo.parent->m[0][0] = dinfo.parent->m[0][1] = temp_si->s_trans.pop_grow.f_m;
		dinfo.parent->Linf[0][0] = dinfo.parent->Linf[0][1] = temp_si->s_trans.pop_grow.f_Linf;
		dinfo.parent->Ti[0][0] = dinfo.parent->Ti[0][1] = temp_si->s_trans.pop_grow.Ti;
		if (temp_si->s_trans.pop_grow.f_m_sd != -9) { //iiv in the slope of the linear growth function, ie faster/slower growers
			normal_distribution<float>m_dist(temp_si->s_trans.pop_grow.f_m, temp_si->s_trans.pop_grow.f_m_sd);
			dinfo.parent->m[0][0] = dinfo.parent->m[0][1] = m_dist(eng);
		}
		if (temp_si->s_trans.pop_grow.f_Linf_sd != -9) { //iiv in the slope of the linear growth function, ie faster/slower growers
			normal_distribution<float>L_dist(temp_si->s_trans.pop_grow.f_Linf, temp_si->s_trans.pop_grow.f_Linf_sd);
			dinfo.parent->Linf[0][0] = dinfo.parent->Linf[0][1] = L_dist(eng);
		}
		if (dinfo.grow_info.Ti_sd != -9) {
			normal_distribution<float>Ti_dist(dinfo.grow_info.Ti, dinfo.grow_info.Ti_sd);
			dinfo.parent->Ti[0][0] = dinfo.parent->Ti[0][1] = Ti_dist(eng);
		}
			
		if(temp_si->s_trans.sexdep==1) { //do the male values
			
			dinfo.parent->m[1][0] = dinfo.parent->m[1][1] = temp_si->s_trans.pop_grow.m_m;
			dinfo.parent->Linf[1][0] = dinfo.parent->Linf[1][1] = temp_si->s_trans.pop_grow.m_Linf;
			dinfo.parent->Ti[1][0] = dinfo.parent->Ti[1][1] = temp_si->s_trans.pop_grow.Ti;
			if (temp_si->s_trans.pop_grow.m_m_sd != -9) { //iiv in the slope of the linear growth function, ie faster/slower growers
				normal_distribution<float>m_dist(temp_si->s_trans.pop_grow.m_m, temp_si->s_trans.pop_grow.m_m_sd);
				dinfo.parent->m[1][0] = dinfo.parent->m[1][1] = m_dist(eng);
			}
			if (temp_si->s_trans.pop_grow.m_Linf_sd != -9) { //iiv in the slope of the linear growth function, ie faster/slower growers
				normal_distribution<float>L_dist(temp_si->s_trans.pop_grow.m_Linf, temp_si->s_trans.pop_grow.m_Linf_sd);
				dinfo.parent->Linf[1][0] = dinfo.parent->Linf[1][1] = L_dist(eng);
			}
			if (temp_si->s_trans.pop_grow.Ti_sd != -9) {
				normal_distribution<float>Ti_dist(temp_si->s_trans.pop_grow.Ti, temp_si->s_trans.pop_grow.Ti_sd);
				dinfo.parent->Ti[1][0] = dinfo.parent->Ti[1][1] = Ti_dist(eng);
			}
		}
		
	}

	if (pstage->stage==0 && dinfo.grow_info.method == -9) { cout << "didnt set grow params right at init" << endl; }
	if (pstage->stage == 0 && dinfo.grow_info.Ti==0) { cout << "didnt set Ti right at init" << endl; }

}

void Individual::set_trans_info(bool invent) {
	if (pstage->s_trans.sexdep == 1 && sex == 0) { //sex dependent and male
		dinfo.rho = pstage->s_trans.m_rho;
		dinfo.rho_sd = pstage->s_trans.m_rho_sd;
		dinfo.mortprob = pstage->s_trans.m_mortprob;
		dinfo.slope = pstage->s_trans.m_slope;
		dinfo.inflpt = pstage->s_trans.m_inflpoint;
		dinfo.SL = pstage->s_trans.m_SL;
		dinfo.active_SL = pstage->s_trans.m_comp_SL;
		dinfo.SL_sd = pstage->s_trans.m_SL_sd;
		dinfo.min_active_size = pstage->s_trans.m_min_active_size;
		dinfo.min_active_time = pstage->s_trans.m_min_active_time;
		dinfo.active_sd = pstage->s_trans.m_active_sd;
	}
	else {
		dinfo.rho = pstage->s_trans.f_rho;
		dinfo.rho_sd = pstage->s_trans.f_rho_sd;
		dinfo.mortprob = pstage->s_trans.f_mortprob;
		dinfo.slope = pstage->s_trans.f_slope;
		dinfo.inflpt = pstage->s_trans.f_inflpoint;
		dinfo.SL = pstage->s_trans.f_SL;
		dinfo.active_SL = pstage->s_trans.f_comp_SL;
		dinfo.SL_sd = pstage->s_trans.f_SL_sd;
		dinfo.min_active_size = pstage->s_trans.f_min_active_size;
		dinfo.min_active_time = pstage->s_trans.f_min_active_time;
		dinfo.active_sd = pstage->s_trans.f_active_sd;

	}
	/*dinfo.active_evo = pstage->s_trans.active_evo;
	dinfo.active_mut = pstage->s_trans.active_mut;
	dinfo.active_mut_size = pstage->s_trans.active_mut_size;*/
	dinfo.buoy_min = pstage->s_trans.buoy_min;
	dinfo.buoy_max = pstage->s_trans.buoy_max;
	dinfo.diel_vert = pstage->s_trans.diel_vert;
	dinfo.dv_range = pstage->s_trans.dv_range;
	dinfo.min_dv_size = pstage->s_trans.min_dv_size;
	dinfo.min_dv_time = pstage->s_trans.min_dv_time;
	dinfo.min_dv_sd = pstage->s_trans.min_dv_sd;
	dinfo.down_bias = pstage->s_trans.down_bias;
	dinfo.memory = pstage->s_trans.memory;

	//if (pstage->stage==0 && dinfo.min_active_size == -9) { cout << "min active size wasn't set right and initialisation" << endl; }
	//if (pstage->stage==0 && dinfo.min_dv_size == -9) { cout << "min dv size wasn't set right and initialisation" << endl; }
	//mortprob is in probability/day so i need to make it per hour
	//dinfo.mortprob /= 24;

	if (mod->size_or_time == 0) { //size dependent behaviour, then SL_m is used to calculate SL
		//then dinfo.SL and dinfo.active_SL are in BL/s, make them BL/hr
		dinfo.rho_m = dinfo.rho;
		//passive disperser
		if (dinfo.min_active_size == -9) { //they dont get active. if they're hybrid, this will have a value and if they're active it will be the juvenile release size
			dinfo.mode = 0; //passive
			dinfo.phase = 0;
		}
		//hybrid disperser
		else if (dinfo.active_SL != -9) { //first passive then has a new SL
			dinfo.active_SL = dinfo.active_SL * 60 * 60; //this is now BL/hr
			dinfo.SL_m = dinfo.active_SL;

			dinfo.mode = 1; //hybrid
			dinfo.phase = 0; //passive
		}

		//active disperser
		else {
			dinfo.SL_m = dinfo.SL * 60 * 60; //BL/hr

			dinfo.mode = 2; //active
			dinfo.phase = 1; //active
		}

		//call the function that sets growth parameters 
		if (invent == true) {
			set_growth(true);
		}
		else {
			set_growth(false);
		}
		
	}
	else { //time dependent behaviour, the SL is the distance in cm/s
		//passive
		if (dinfo.min_active_time == -9) { //they dont get active. if they're hybrid, this will have a value and if they're active it will be the juvenile release size
			dinfo.mode = 0; //passive
			dinfo.phase = 0; //passive
		}
		//hybrid disperser
		else if (dinfo.active_SL != -9) { //first passive then has a new SL
			dinfo.active_SL = dinfo.active_SL * 60 * 60 / 100; //this is now m/hr
			dinfo.mode = 1; //hybrid
			dinfo.phase = 0; //starts out passive;
		}

		//active disperser
		else {
			
			dinfo.SL = dinfo.SL * 60 * 60 / 100; //input is given in cm/s so make it m/hr

			dinfo.mode = 2; //active
			dinfo.phase = 1; //active
			
		}
	}


	if (mod->size_or_time == 0) { //if behaviour is size dependent, now that i have sampled the iiv, calculation SL and rho from the slope values
		calc_SL_rho();
	}

	//cout << "indiv's SL is " << dinfo.SL;
	//cout << "indiv's rho is " << dinfo.rho; 
	//cout << "indiv's active_SL is " << dinfo.active_SL; 
	//cout << "indiv's buoy min is" << dinfo.buoy_min;
	//cout << "indiv's diel vert is " << dinfo.diel_vert;

	if (invent == true && pstage->s_emig.stagedep == 0) { //if this individual is part of the first generation of the model, then it doesn't have any parent data so i need to invent some
		//also dispersal needs to be stage independent for parental data to matter. if stage dependent, then can't be iiv or evo so no recombination etc
		stage_infob* temp_si = (*subpop_stages)[pstage->s_emig.emigstage];
		//transfer parameters
		float rsamp;
		if (mod->size_or_time == 0) { //size dep
			
			dinfo.parent->rho_m[0][0] = dinfo.parent->rho_m[0][1] = temp_si->s_trans.f_rho;
			//passive
			if (temp_si->s_trans.f_min_active_size == -9) {
				dinfo.parent->SL_m[0][0] = dinfo.parent->SL_m[0][1] = -9;
			}
			//hybrid
			else if (temp_si->s_trans.f_comp_SL != -9) { 
				dinfo.parent->SL_m[0][0] = dinfo.parent->SL_m[0][1] = (temp_si->s_trans.f_comp_SL * 60 * 60); //this is now in BL/hr, this is slope of size dependence equation
				//dinfo.parent->active_SL[0][0] = dinfo.parent->active_SL[0][1] = (temp_si->s_trans.f_comp_SL * 60 * 60); 
			}
			//active
			else { //they don't switch so use SL
				dinfo.parent->SL_m[0][0] = dinfo.parent->SL_m[0][1] = (temp_si->s_trans.f_SL * 60 * 60);
				//dinfo.parent->active_SL[0][0] = dinfo.parent->active_SL[0][1] = (temp_si->s_trans.f_comp_SL * 60 * 60 / 100);
			}
			dinfo.parent->min_active_size[0][0] = dinfo.parent->min_active_size[0][1] = temp_si->s_trans.f_min_active_size;
			dinfo.parent->min_dv_size[0][0] = dinfo.parent->min_dv_size[0][1] = temp_si->s_trans.min_dv_size;
			if (temp_si->s_trans.f_rho_sd != -9) { //if it is applied in rho (ie more or less likely to choose to go against the current
				normal_distribution<float>rho_dist(temp_si->s_trans.f_rho, temp_si->s_trans.f_rho_sd);
				rsamp = rho_dist(eng);
				if (rsamp < 0) { rsamp = 0; }
				if (rsamp > 1) { rsamp = 1; }
				dinfo.parent->rho_m[0][0] = dinfo.parent->rho_m[0][1] = rsamp;
			}
			if (temp_si->s_trans.f_SL_sd != -9) { //if it is applied to step length (ie stronger and weaker swimmers)
				normal_distribution<float>SL_dist;

				//if comp_SL has a value, then SL =1 until indiv becomes active, there is no iiv there, so can just sample for comp_SL
				if (temp_si->s_trans.f_comp_SL != -9) {
					SL_dist = normal_distribution<float>(temp_si->s_trans.f_comp_SL, temp_si->s_trans.f_SL_sd);
				}
				else { //otherwise, SL can change now, indiv is always active (if indiv were always passive, there is no iiv)
					SL_dist = normal_distribution<float>(temp_si->s_trans.f_SL, temp_si->s_trans.f_SL_sd);
				}	
				//dont overwrite dinfo.SL because it needs to remain ==1 until indiv is the right size/age
				dinfo.parent->SL_m[0][0] = dinfo.parent->SL_m[0][1] = SL_dist(eng); //this will then be used to calculate SL based on body size	
			}
			if (temp_si->s_trans.f_active_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
				normal_distribution<float>active_dist;					
				active_dist = normal_distribution<float>(temp_si->s_trans.f_min_active_size, temp_si->s_trans.f_active_sd);
				dinfo.parent->min_active_size[0][0]=dinfo.parent->min_active_size[0][1] = active_dist(eng);
			}
			if (temp_si->s_trans.min_dv_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
				normal_distribution<float>dv_dist;
				dv_dist = normal_distribution<float>(temp_si->s_trans.min_dv_size, temp_si->s_trans.min_dv_sd);
				dinfo.parent->min_dv_size[0][0] = dinfo.parent->min_dv_size[0][1] = dv_dist(eng);
			}
				
			/*dinfo.parent->m[0][0] = dinfo.parent->m[0][1] = temp_si->s_trans.pop_grow.f_m;
			dinfo.parent->Linf[0][0] = dinfo.parent->Linf[0][1] = temp_si->s_trans.pop_grow.f_Linf;*/
			if(temp_si->s_trans.sexdep==1) { //do the male values
				dinfo.parent->rho_m[1][0] = dinfo.parent->rho_m[1][1] = temp_si->s_trans.m_rho;
				//passive
				if (temp_si->s_trans.m_min_active_size == -9) {
					dinfo.parent->SL_m[1][0] = dinfo.parent->SL_m[1][1] = -9;
				}
				//hybrid
				if (temp_si->s_trans.m_comp_SL != -9) { //if they switch to active behaviour after a certain size, which means dinfo.SL==1
					dinfo.parent->SL_m[1][0] = dinfo.parent->SL_m[1][1] = (temp_si->s_trans.m_comp_SL * 60 * 60); //this is now in BL/hr
					//dinfo.parent->active_SL[1][0] = dinfo.parent->active_SL[1][1] = (temp_si->s_trans.m_comp_SL * 60 * 60);
				}
				//active
				else { //they don't switch so use SL
					dinfo.parent->SL_m[1][0] = dinfo.parent->SL_m[1][1] = (temp_si->s_trans.m_SL * 60 * 60);
					//dinfo.parent->active_SL[1][0] = dinfo.parent->SL[1][1] = (temp_si->s_trans.m_comp_SL * 60 * 60 / 100);
				}
				dinfo.parent->min_active_size[1][0] = dinfo.parent->min_active_size[1][1] = temp_si->s_trans.m_min_active_size;
				dinfo.parent->min_dv_size[1][0] = dinfo.parent->min_dv_size[1][1] = temp_si->s_trans.min_dv_size;
				/*dinfo.parent->m[1][0] = dinfo.parent->m[1][1] = temp_si->s_trans.pop_grow.m_m;
				dinfo.parent->Linf[1][0] = dinfo.parent->Linf[1][1] = temp_si->s_trans.pop_grow.m_Linf;*/
				if (temp_si->s_trans.m_rho_sd != -9) { //if it is applied in rho (ie more or less likely to choose to go against the current
					normal_distribution<float>rho_dist(temp_si->s_trans.m_rho, temp_si->s_trans.m_rho_sd);
					rsamp = rho_dist(eng);
					if (rsamp < 0) { rsamp = 0; }
					if (rsamp > 1) { rsamp = 1; }
					dinfo.parent->rho_m[1][0] = dinfo.parent->rho_m[1][1] = rsamp;
				}
				if (temp_si->s_trans.m_SL_sd != -9) { //if it is applied to step length (ie stronger and weaker swimmers)
					normal_distribution<float>SL_dist;

					//if comp_SL has a value, then SL =1 until indiv becomes active, there is no iiv there, so can just sample for comp_SL
					if (temp_si->s_trans.m_comp_SL != -9) {
						SL_dist = normal_distribution<float>(temp_si->s_trans.m_comp_SL, temp_si->s_trans.m_SL_sd);
									
					}
					else { //otherwise, SL can change now, indiv is always active (if indiv were always passive, there is no iiv)
						SL_dist = normal_distribution<float>(temp_si->s_trans.m_SL, temp_si->s_trans.m_SL_sd);
					}
					//dont overwrite dinfo.SL because it needs to remain ==1 until indiv is the right size/age
					dinfo.parent->SL_m[1][0] = dinfo.parent->SL_m[1][1] = SL_dist(eng); //this will then be used to calculate SL based on body size
				}
				if (temp_si->s_trans.m_active_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
					normal_distribution<float>active_dist;
					active_dist = normal_distribution<float>(temp_si->s_trans.m_min_active_size, temp_si->s_trans.m_active_sd);
					dinfo.parent->min_active_size[1][0] = dinfo.parent->min_active_size[1][1] = active_dist(eng);
				}
				if (temp_si->s_trans.min_dv_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
					normal_distribution<float>dv_dist;
					dv_dist = normal_distribution<float>(temp_si->s_trans.min_dv_size, temp_si->s_trans.min_dv_sd);
					dinfo.parent->min_dv_size[1][0] = dinfo.parent->min_dv_size[1][1] = dv_dist(eng);
				}
			}

			

		}
		else {//time dependent
			//passive
			if (temp_si->s_trans.f_min_active_time == -9) {
				dinfo.parent->SL[0][0] = dinfo.parent->SL[0][1] = 1;
				dinfo.parent->active_SL[0][0] = dinfo.parent->active_SL[0][1] = -9;
			}
			//hybrid
			else if (temp_si->s_trans.f_comp_SL != -9) { //if they switch to active behaviour after a certain size, which means dinfo.SL==1
				dinfo.parent->SL[0][0] = dinfo.parent->SL[0][1] = 1;
				dinfo.parent->active_SL[0][0] = dinfo.parent->active_SL[0][1] = (temp_si->s_trans.f_comp_SL * 60 * 60 / 100);
			}
			else {
				dinfo.parent->SL[0][0] = dinfo.parent->SL[0][1] = (temp_si->s_trans.f_SL * 60 * 60 / 100);
				dinfo.parent->active_SL[0][0] = dinfo.parent->active_SL[0][1] = -9;
			}
			dinfo.parent->rho[0][0] = dinfo.parent->rho[0][1] = temp_si->s_trans.f_rho;
			dinfo.parent->min_active_time[0][0] = dinfo.parent->min_active_time[0][1] = temp_si->s_trans.f_min_active_time;
			dinfo.parent->min_dv_time[0][0] = dinfo.parent->min_dv_time[0][1] = temp_si->s_trans.min_dv_time;
			if (temp_si->s_trans.f_rho_sd != -9) { //if it is applied in rho (ie more or less likely to choose to go against the current
				normal_distribution<float>rho_dist(temp_si->s_trans.f_rho, temp_si->s_trans.f_rho_sd);
				rsamp = rho_dist(eng);
				if (rsamp < 0) { rsamp = 0; }
				if (rsamp > 1) { rsamp = 1; }
				dinfo.parent->rho[0][0] = dinfo.parent->rho[0][1] = rsamp;
			}
			if (temp_si->s_trans.f_SL_sd != -9) { //if it is applied to step length (ie stronger and weaker swimmers)
				normal_distribution<float>SL_dist;

				//if comp_SL has a value, then SL =1 until indiv becomes active, there is no iiv there, so can just sample for comp_SL
				if (temp_si->s_trans.f_comp_SL != -9) {
					SL_dist = normal_distribution<float>(temp_si->s_trans.f_comp_SL, temp_si->s_trans.f_SL_sd);
					//dont overwrite dinfo.SL because it needs to remain ==1 until indiv is the right size/age
					dinfo.parent->active_SL[0][0] = dinfo.parent->active_SL[0][1] = SL_dist(eng); //this will then be used to calculate SL based on body size	
				}
				else { //otherwise, SL can change now, indiv is always active (if indiv were always passive, there is no iiv)
					SL_dist = normal_distribution<float>(temp_si->s_trans.f_SL, temp_si->s_trans.f_SL_sd);
					dinfo.parent->SL[0][0] = dinfo.parent->SL[0][1] = SL_dist(eng); //this will then be used to calculate SL based on body size
				}
				
			}
			if (temp_si->s_trans.f_active_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
				normal_distribution<float>active_dist;
				active_dist = normal_distribution<float>(temp_si->s_trans.f_min_active_time, temp_si->s_trans.f_active_sd);
				dinfo.parent->min_active_time[0][0] = dinfo.parent->min_active_time[0][1] = active_dist(eng);
			}
			if (temp_si->s_trans.min_dv_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
				normal_distribution<float>dv_dist;
				dv_dist = normal_distribution<float>(temp_si->s_trans.min_dv_time, temp_si->s_trans.min_dv_sd);
				dinfo.parent->min_dv_time[1][0] = dinfo.parent->min_dv_size[0][1] = dv_dist(eng);
			}
			
			if (temp_si->s_trans.sexdep == 1) { //do the male values
				//passive
				if (temp_si->s_trans.m_min_active_time == -9) {
					dinfo.parent->SL[1][0] = dinfo.parent->SL[1][1] = 1;
					dinfo.parent->active_SL[1][0] = dinfo.parent->active_SL[1][1] = -9;
				}

				//hybrid
				if (temp_si->s_trans.m_comp_SL != -9) { //if they switch to active behaviour after a certain size, which means dinfo.SL==1
					dinfo.parent->SL[1][0] = dinfo.parent->SL[1][1] = 1;
					dinfo.parent->active_SL[1][0] = dinfo.parent->active_SL[1][1] = (temp_si->s_trans.m_comp_SL * 60 * 60 / 100);
				}
				//active
				else {
					dinfo.parent->SL[1][0] = dinfo.parent->SL[1][1] = (temp_si->s_trans.m_SL * 60 * 60 / 100);
					dinfo.parent->active_SL[1][0] = dinfo.parent->active_SL[1][1] = -9;
				}
				dinfo.parent->rho[1][0] = dinfo.parent->rho[1][1] = temp_si->s_trans.m_rho;
				dinfo.parent->min_active_time[1][0] = dinfo.parent->min_active_time[1][1] = temp_si->s_trans.m_min_active_time;
				dinfo.parent->min_dv_time[1][0] = dinfo.parent->min_dv_time[1][1] = temp_si->s_trans.min_dv_time;
				if (temp_si->s_trans.f_rho_sd != -9) { //if it is applied in rho (ie more or less likely to choose to go against the current
					normal_distribution<float>rho_dist(temp_si->s_trans.m_rho, temp_si->s_trans.m_rho_sd);
					rsamp = rho_dist(eng);
					if (rsamp < 0) { rsamp = 0; }
					if (rsamp > 1) { rsamp = 1; }
					dinfo.parent->rho[1][0] = dinfo.parent->rho[1][1] = rsamp;
				}
				if (temp_si->s_trans.m_SL_sd != -9) { //if it is applied to step length (ie stronger and weaker swimmers)
					normal_distribution<float>SL_dist;

					//if comp_SL has a value, then SL =1 until indiv becomes active, there is no iiv there, so can just sample for comp_SL
					if (temp_si->s_trans.m_comp_SL != -9) {
						SL_dist = normal_distribution<float>(temp_si->s_trans.m_comp_SL, temp_si->s_trans.m_SL_sd);
						//dont overwrite dinfo.SL because it needs to remain ==1 until indiv is the right size/age
						dinfo.parent->active_SL[1][0] = dinfo.parent->active_SL[1][1] = SL_dist(eng); //this will then be used to calculate SL based on body size	
					}
					else { //otherwise, SL can change now, indiv is always active (if indiv were always passive, there is no iiv)
						SL_dist = normal_distribution<float>(temp_si->s_trans.m_SL, temp_si->s_trans.m_SL_sd);
						dinfo.parent->SL[1][0] = dinfo.parent->SL[1][1] = SL_dist(eng); //this will then be used to calculate SL based on body size
					}

				}
				if (temp_si->s_trans.m_active_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
					normal_distribution<float>active_dist;
					active_dist = normal_distribution<float>(temp_si->s_trans.m_min_active_size, temp_si->s_trans.m_active_sd);
					dinfo.parent->min_active_time[1][0] = dinfo.parent->min_active_time[1][1] = active_dist(eng);
				}
				if (temp_si->s_trans.min_dv_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
					normal_distribution<float>dv_dist;
					dv_dist = normal_distribution<float>(temp_si->s_trans.min_dv_time, temp_si->s_trans.min_dv_sd);
					dinfo.parent->min_dv_time[1][0] = dinfo.parent->min_dv_size[1][1] = dv_dist(eng);
				}
			}
		}
	}
	//if (dinfo.mortprob!=-9 && dinfo.mortprob != (0.042/24)) { cout << "mort prob is " << dinfo.mortprob << " and indiv is a " << sex << endl; }
}

void Individual::dispersal_mort() {
	/*default_random_engine unif_generator;
	uniform_real_distribution<float> mort_dist(0, 1);*/
	float ran_mort = unif_dist(eng); //sample randomly from a uniform distribution
	//if (pop_tinfo->distmort == 1) { //distance dependent mortality
	//	float dispmort = 1.0 / (1.0 + exp(-(dinfo.disp_distance - dinfo.inflpt)*dinfo.slope));
	//	if (ran_mort < dispmort) {
	//		dead = true;
	//	}
	//}
	//else {
	if (ran_mort < dinfo.mortprob) { //if the random mortality rate is less than the constant mortality probability
		//cout << "indiv died from mortprob " << dinfo.mortprob << endl;
		status = 5; dead = true; //individual dies
		//cout << "indiv " << ID << " has died due to disp mort" << endl;
	}
	//}
}
//void Individual::set_trans_info() {
//	
//	if (pop_tinfo->stagedep == 1) { //if transfer is stage dependent
//		if (pop_tinfo->sexdep == 1 && sex == 0) { //and sex dependent, assign male info
//			dinfo.rho = pstage->m_rho;
//			dinfo.rho_sd = pstage->m_rho_sd;
//			dinfo.mortprob = pstage->m_mortprob;
//			dinfo.slope = pstage->m_slope;
//			dinfo.inflpt = pstage->m_inflpoint;
//			dinfo.SL = pstage->m_SL;
//			dinfo.active_SL = pstage->m_comp_SL;
//			dinfo.SL_sd = pstage->m_SL_sd;
//			dinfo.min_active_size = pstage->m_min_active_size;
//			dinfo.min_active_time = pstage->m_min_active_time;
//			dinfo.active_sd = pstage->m_active_sd;
//		}
//		else { //else assign female info
//			dinfo.rho = pstage->f_rho;
//			dinfo.rho_sd = pstage->f_rho_sd;
//			dinfo.mortprob = pstage->f_mortprob;
//			dinfo.slope = pstage->f_slope;
//			dinfo.inflpt = pstage->f_inflpoint;
//			dinfo.SL = pstage->f_SL;
//			dinfo.active_SL = pstage->f_comp_SL;
//			dinfo.SL_sd = pstage->f_SL_sd;
//			dinfo.min_active_size = pstage->f_min_active_size;
//			dinfo.min_active_time = pstage->f_min_active_time;
//			dinfo.active_sd = pstage->f_active_sd;
//		}
//		dinfo.memory = pstage->memory;
//		dinfo.dp = pstage->dp;
//		dinfo.diel_vert = pstage->diel_vert;
//		dinfo.dv_range = pstage->dv_range;
//		dinfo.dv_active = pstage->dv_active;
//		dinfo.min_dv_time = pstage->min_dv_time * 24; //input is in days, need it in hours
//		dinfo.min_dv_size = pstage->min_dv_size;
//		dinfo.down_bias = pstage->down_bias;
//		dinfo.buoy_min = pstage->buoy_min;
//		dinfo.buoy_max = pstage->buoy_max;
//		if (dinfo.min_active_time != -9) { dinfo.min_active_time *= 24; } //input is in days, need it in hours
//		dinfo.active_evo = pstage->active_evo;
//		dinfo.active_mut = pstage->active_mut;
//		dinfo.active_mut_size = pstage->active_mut_size;
//	}
//	else { //if transfer is not stage dependent
//		if (pop_tinfo->sexdep == 1 && sex == 0) { //but it is sex dependent, assign male info
//			dinfo.rho = pop_tinfo->m_rho;
//			dinfo.rho_sd = pop_tinfo->m_rho_sd;
//			dinfo.mortprob = pop_tinfo->m_mortprob;
//			dinfo.slope = pop_tinfo->m_slope;
//			dinfo.inflpt = pop_tinfo->m_inflpoint;
//			dinfo.SL = pop_tinfo->m_SL;
//			dinfo.active_SL = pop_tinfo->m_comp_SL;
//			dinfo.SL_sd = pop_tinfo->m_SL_sd;
//			dinfo.min_active_size = pop_tinfo->m_min_active_size;
//			dinfo.min_active_time = pop_tinfo->m_min_active_time;
//			dinfo.active_sd = pop_tinfo->m_active_sd;
//		}
//		else { //else assign female info
//			dinfo.rho = pop_tinfo->f_rho;
//			dinfo.rho_sd = pop_tinfo->f_rho_sd;
//			dinfo.mortprob = pop_tinfo->f_mortprob;
//			dinfo.slope = pop_tinfo->f_slope;
//			dinfo.inflpt = pop_tinfo->f_inflpoint;
//			dinfo.SL = pop_tinfo->f_SL;
//			dinfo.active_SL = pop_tinfo->f_comp_SL;
//			dinfo.SL_sd = pop_tinfo->f_SL_sd;
//			dinfo.min_active_size = pop_tinfo->f_min_active_size;
//			dinfo.min_active_time = pop_tinfo->f_min_active_time;
//			dinfo.active_sd = pop_tinfo->f_active_sd;
//		}
//		dinfo.memory = pop_tinfo->memory;
//		dinfo.dp = pop_tinfo->dp;
//		dinfo.diel_vert = pop_tinfo->diel_vert;
//		dinfo.dv_range = pop_tinfo->dv_range;
//		dinfo.dv_active = pop_tinfo->dv_active;
//		dinfo.min_dv_time = pop_tinfo->min_dv_time * 24; //input is in days, need it in hours
//		dinfo.min_dv_size = pop_tinfo->min_dv_size;
//		dinfo.down_bias = pop_tinfo->down_bias;
//		dinfo.buoy_min = pop_tinfo->buoy_min;
//		dinfo.buoy_max = pop_tinfo->buoy_max;
//		if (dinfo.min_active_time != -9) { dinfo.min_active_time *= 24; } //input is in days, need it in hours
//		dinfo.active_evo = pop_tinfo->active_evo;
//		dinfo.active_mut = pop_tinfo->active_mut;
//		dinfo.active_mut_size = pop_tinfo->active_mut_size;
//	}
//	//mortprob is in probability/day so i need to make it per hour
//	dinfo.mortprob /= 24;
//
//	if (mod->size_or_time == 0) { //size dependent behaviour, then SL_m is used to calculate SL
//		//then dinfo.SL and dinfo.active_SL are in BL/s, make them BL/hr
//		dinfo.rho_m = dinfo.rho;
//		if (dinfo.active_SL != -9) { //if they switch to active behaviour after a certain size, which means dinfo.SL==1
//			dinfo.active_SL = dinfo.active_SL * 60 * 60; //this is now BL/hr
//			dinfo.SL_m = dinfo.active_SL;
//		}
//		else { //otherwise, use SL
//			dinfo.SL_m = dinfo.SL * 60 * 60; //BL/hr
//		}
//		//call the function that sets growth parameters 
//		set_growth();
//	}
//	else { //time dependent behaviour, the SL is the distance in cm/s
//		if (dinfo.active_SL != -9) { //if they switch to active behaviour after a certain time, which means dinfo.SL==1
//			dinfo.active_SL = dinfo.active_SL * 60 * 60 / 100; //this is now m/hr
//		}
//		else { //otherwise, use SL
//			if(dinfo.SL != 1){ //if SL==1, that means it's passive and needs to remain 1 until active_SL (if applicable), or throughout dispersal
//				dinfo.SL = dinfo.SL * 60 * 60 /100 ; //input is given in cm/s so make it m/hr
//			}
//		}
//	}
//	//all step length parameters should now at least be /hr, BL get covnerted to m/hr when SL is calculated
//
//	//since all parameters can be set separately for inter individual variation, check whether each sd ==-9
//	if (dinfo.rho_sd != -9) { //if it is applied in rho (ie more or less likely to choose to go against the current
//		normal_distribution<float>rho_dist(dinfo.rho, dinfo.rho_sd);
//		if (mod->size_or_time == 0) { //if the model is using size-dependent dispersal, therefore growth rate, dinfo.SL is a MULTIPLIER, how many body lengths
//			dinfo.rho_m = rho_dist(eng);
//		}
//		else {
//			dinfo.rho = rho_dist(eng);
//		}
//	}
//	if (dinfo.SL_sd != -9) { //if it is applied to step length (ie stronger and weaker swimmers)
//		normal_distribution<float>SL_dist;
//
//		//if comp_SL has a value, then SL =1 until indiv becomes active, there is no iiv there, so can just sample for comp_SL
//		if (dinfo.active_SL != -9) {
//			SL_dist = normal_distribution<float>(dinfo.active_SL, dinfo.SL_sd);
//			//dont overwrite dinfo.SL because it needs to remain ==1 until indiv is the right size/age
//			if (mod->size_or_time == 0) { //size dependent
//				dinfo.SL_m = SL_dist(eng); //this will then be used to calculate SL based on body size
//			}
//			else { //time dependent
//				dinfo.active_SL = SL_dist(eng);
//			}
//		}
//		else { //otherwise, SL can change now, indiv is always active (if indiv were always passive, there is no iiv)
//			SL_dist = normal_distribution<float>(dinfo.SL, dinfo.SL_sd);
//			if (mod->size_or_time == 0) { //size dependent
//				dinfo.SL_m = SL_dist(eng); //this will then be used to calculate SL based on body size
//			}
//			else { //time dependent
//				dinfo.SL = SL_dist(eng);
//			}
//		}
//	}
//	if (dinfo.active_sd != -9) { //iiv in size/time to active movement during dispersal, ie developed at smaller size
//		normal_distribution<float>active_dist;
//		if (mod->size_or_time == 0) { //if behaviour is size dependent
//			active_dist = normal_distribution<float>(dinfo.min_active_size, dinfo.active_sd);
//			dinfo.min_active_size = active_dist(eng);
//		}
//		else { //if it is time dependent
//			active_dist = normal_distribution<float>(dinfo.min_active_time, dinfo.active_sd);
//			dinfo.min_active_time = active_dist(eng);
//		}
//	}
//
//	if (mod->size_or_time == 0) { //if behaviour is size dependent, now that i have sampled the iiv, calculation SL and rho from the slope values
//		calc_SL_rho();
//	}
//	
//	//cout << "indiv's SL is " << dinfo.SL;
//	//cout << "indiv's rho is " << dinfo.rho; 
//	//cout << "indiv's active_SL is " << dinfo.active_SL; 
//	//cout << "indiv's buoy min is" << dinfo.buoy_min;
//	//cout << "indiv's diel vert is " << dinfo.diel_vert;
//}
//
//void Individual::dispersal_mort() {
//	/*default_random_engine unif_generator;
//	uniform_real_distribution<float> mort_dist(0, 1);*/
//	float ran_mort = unif_dist(eng); //sample randomly from a uniform distribution
//	//if (pop_tinfo->distmort == 1) { //distance dependent mortality
//	//	float dispmort = 1.0 / (1.0 + exp(-(dinfo.disp_distance - dinfo.inflpt)*dinfo.slope));
//	//	if (ran_mort < dispmort) {
//	//		dead = true;
//	//	}
//	//}
//	//else {
//		if (ran_mort < dinfo.mortprob) { //if the random mortality rate is less than the constant mortality probability
//			status = 5; dead = true; //individual dies
//			cout << "indiv " << ID << " has died due to disp mort" << endl;
//		}
//	//}
//}

//settlement
void Individual::set_settle_info(bool invent) {
	if (pstage->s_sett.sexdep == 1 && sex == 0) { //and sex dependent, assign male info
		dinfo.S0 = pstage->s_sett.m_S0;
		dinfo.alphaS = pstage->s_sett.m_alphaS;
		dinfo.betaS = pstage->s_sett.m_betaS;
		dinfo.sett_sd = pstage->s_sett.m_sett_sd;
		dinfo.pld = pstage->s_sett.m_pld;
		dinfo.pld_sd = pstage->s_sett.m_pld_sd;
		dinfo.comp_time = pstage->s_sett.m_comp;
		dinfo.comp_size = pstage->s_sett.m_comp_size;
		dinfo.comp_sd = pstage->s_sett.m_comp_sd;
		dinfo.findmate = pstage->s_sett.m_findmate;
	}
	else {
		dinfo.S0 = pstage->s_sett.f_S0;
		dinfo.alphaS = pstage->s_sett.f_alphaS;
		dinfo.betaS = pstage->s_sett.f_betaS;
		dinfo.sett_sd = pstage->s_sett.f_sett_sd;
		dinfo.pld = pstage->s_sett.f_pld;
		dinfo.pld_sd = pstage->s_sett.f_pld_sd;
		dinfo.comp_time = pstage->s_sett.f_comp;
		dinfo.comp_size = pstage->s_sett.f_comp_size;
		dinfo.comp_sd = pstage->s_sett.f_comp_sd;
		dinfo.findmate = pstage->s_sett.f_findmate;
	}

	/*dinfo.S0_evo = pstage->s_sett.S0_evo;
	dinfo.S0_mut = pstage->s_sett.S0_mut;
	dinfo.S0_mut_size = pstage->s_sett.S0_mut_size;
	dinfo.alphaS_mut_size = pstage->s_sett.alphaS_mut_size;
	dinfo.betaS_mut_size = pstage->s_sett.betaS_mut_size;
	dinfo.comp_evo = pstage->s_sett.comp_evo;
	dinfo.comp_mut = pstage->s_sett.comp_mut;
	dinfo.comp_mut_size = pstage->s_sett.comp_mut_size;*/



	if (invent == true && pstage->s_emig.stagedep == 0) { //if this individual is part of the first generation of the model, then it doesn't have any parent data so i need to invent some
		//also dispersal needs to be stage independent for parental data to matter. if stage dependent, then can't be iiv or evo so no recombination etc
		stage_infob* temp_si = (*subpop_stages)[pstage->s_emig.emigstage];
		float rsamp;
		dinfo.parent->S0[0][0] = dinfo.parent->S0[0][1] = temp_si->s_sett.f_S0;
		
		if (temp_si->s_sett.densdep == 1) {
			dinfo.parent->alphaS[0][0] = dinfo.parent->alphaS[0][1] = temp_si->s_sett.f_alphaS;
			dinfo.parent->betaS[0][0] = dinfo.parent->betaS[0][1] = temp_si->s_sett.f_betaS;
			//cout << "assigned dinfo.parent.alphaS: " << dinfo.parent->alphaS[0][0] << ", " << dinfo.parent->alphaS[0][1];
			//cout << " assigned dinfo.parent.betaS: " << dinfo.parent->betaS[0][0] << ", " << dinfo.parent->betaS[0][1] << endl;
		}
		if (temp_si->s_sett.f_sett_sd != -9) {
			normal_distribution<float> S0_dist(temp_si->s_sett.f_S0, temp_si->s_sett.f_sett_sd);
			rsamp = S0_dist(eng);
			if (rsamp < 0) { rsamp = 0; }
			if (rsamp > 1) { rsamp = 1; }
			dinfo.parent->S0[0][0] = dinfo.parent->S0[0][1] = rsamp;
			if (temp_si->s_sett.densdep == 1) {
				normal_distribution<float> alpha_dist(temp_si->s_sett.f_alphaS, temp_si->s_sett.f_sett_sd);
				normal_distribution<float> beta_dist(temp_si->s_sett.f_betaS, temp_si->s_sett.f_sett_sd);
				dinfo.parent->alphaS[0][0] = dinfo.parent->alphaS[0][1] = alpha_dist(eng);
				dinfo.parent->betaS[0][0] = dinfo.parent->betaS[0][1] = beta_dist(eng);
			}
		}
		/*if (temp_si->s_sett.f_pld_sd != -9) {
			normal_distribution<float> pld_dist(temp_si->s_sett.f_pld, temp_si->s_sett.f_pld_sd);
			dinfo.parent->pld[0][0]=dinfo.parent->pld[0][1] = pld_dist(eng);
		}*/
		if (mod->size_or_time == 0) {
			dinfo.parent->comp_size[0][0] = dinfo.parent->comp_size[0][1] = temp_si->s_sett.f_comp_size;
			if (temp_si->s_sett.f_comp_sd != -9) {
				normal_distribution<float> comp_dist;
				comp_dist = normal_distribution<float>(temp_si->s_sett.f_comp_size, temp_si->s_sett.f_comp_sd);
				dinfo.parent->comp_size[0][0] = dinfo.parent->comp_size[0][1] = comp_dist(eng);
			}
		}
		else {
			dinfo.parent->comp_time[0][0] = dinfo.parent->comp_time[0][1] = temp_si->s_sett.f_comp;
			if (temp_si->s_sett.f_comp_sd != -9) {
				normal_distribution<float> comp_dist;
				comp_dist = normal_distribution<float>(temp_si->s_sett.f_comp, temp_si->s_sett.f_comp_sd);
				dinfo.parent->comp_time[0][0] = dinfo.parent->comp_time[0][1] = comp_dist(eng);
			}
		}

		if (temp_si->s_sett.sexdep == 1) { //sex dependent, so need male values in a 2x2 matrix
			dinfo.parent->S0[1][0] = dinfo.parent->S0[1][1] = temp_si->s_sett.m_S0;
			
			if (temp_si->s_sett.densdep == 1) {
				dinfo.parent->alphaS[1][0] = dinfo.parent->alphaS[1][1] = temp_si->s_sett.m_alphaS;
				dinfo.parent->betaS[1][0] = dinfo.parent->betaS[1][1] = temp_si->s_sett.m_betaS;
			}
			if (temp_si->s_sett.m_sett_sd != -9) {
				normal_distribution<float> S0_dist(temp_si->s_sett.m_S0, temp_si->s_sett.m_sett_sd);
				rsamp = S0_dist(eng);
				if (rsamp < 0) { rsamp = 0; }
				if (rsamp > 1) { rsamp = 1; }
				dinfo.parent->S0[1][0] = dinfo.parent->S0[1][1] = rsamp;
				if (temp_si->s_sett.densdep == 1) {
					normal_distribution<float> alpha_dist(temp_si->s_sett.m_alphaS, temp_si->s_sett.m_sett_sd);
					normal_distribution<float> beta_dist(temp_si->s_sett.m_betaS, temp_si->s_sett.m_sett_sd);
					dinfo.parent->alphaS[0][0] = dinfo.parent->alphaS[0][1] = alpha_dist(eng);
					dinfo.parent->betaS[1][0] = dinfo.parent->betaS[1][1] = beta_dist(eng);
				}
			}
			if (mod->size_or_time == 0) {
				dinfo.parent->comp_size[1][0] = dinfo.parent->comp_size[1][1] = temp_si->s_sett.m_comp_size;
				if (temp_si->s_sett.m_comp_sd != -9) {
					normal_distribution<float> comp_dist;
					comp_dist = normal_distribution<float>(temp_si->s_sett.m_comp_size, temp_si->s_sett.m_comp_sd);
					dinfo.parent->comp_size[1][0] = dinfo.parent->comp_size[1][1] = comp_dist(eng);
				}
			}
			else {
				dinfo.parent->comp_time[1][0] = dinfo.parent->comp_time[1][1] = temp_si->s_sett.m_comp;
				if (temp_si->s_sett.m_comp_sd != -9) {
					normal_distribution<float> comp_dist;
					comp_dist = normal_distribution<float>(temp_si->s_sett.m_comp, temp_si->s_sett.m_comp_sd);
					dinfo.parent->comp_time[1][0] = dinfo.parent->comp_time[1][1] = comp_dist(eng);
				}
			}
		}
		

	}
	////interindividual variation is set separately, so need to check for it separately
	//if (first == true) {
	//	if (dinfo.sett_sd != -9) {
	//		normal_distribution<float> S0_dist(dinfo.S0, dinfo.sett_sd);
	//		normal_distribution<float> alpha_dist(dinfo.alphaS, dinfo.sett_sd);
	//		normal_distribution<float> beta_dist(dinfo.betaS, dinfo.sett_sd);
	//		dinfo.S0 = S0_dist(eng);
	//		dinfo.alphaS = alpha_dist(eng);
	//		dinfo.betaS = beta_dist(eng);
	//	}
	//	if (dinfo.pld_sd != -9) {
	//		normal_distribution<float> pld_dist(dinfo.pld, dinfo.pld_sd);
	//		dinfo.pld = pld_dist(eng);
	//	}
	//	if (dinfo.comp_sd != -9) {
	//		normal_distribution<float> comp_dist;
	//		if (mod->size_or_time == 0) {
	//			comp_dist = normal_distribution<float>(dinfo.comp_size, dinfo.comp_sd);
	//			dinfo.comp_size = comp_dist(eng);
	//		}
	//		else {
	//			comp_dist = normal_distribution<float>(dinfo.comp_time, dinfo.comp_sd);
	//			dinfo.comp_time = comp_dist(eng);
	//		}
	//	}
	//}
	if (dinfo.pld == 0) { cout << "pld is 0 for some reason" << endl; }
	/*if (pstage->stage == 0) {
		cout << " for indiv " << ID << ", before inheritance S0 is " << dinfo.S0 << ", with alphaS " << dinfo.alphaS << ", with betaS " << dinfo.betaS << endl;
	}*/

}


//void Individual::set_settle_info() {
//	if (pop_sinfo->stagedep == 1) { //if settlement is stage dependent
//		if (pop_sinfo->sexdep == 1 && sex == 0) { //and sex dependent, assign male info
//			dinfo.S0 = pstage->m_S0;
//			dinfo.alphaS = pstage->m_alphaS;
//			dinfo.betaS = pstage->m_betaS;
//			dinfo.sett_sd = pstage->m_sett_sd;
//			dinfo.pld = pstage->m_pld;
//			dinfo.pld_sd = pstage->m_pld_sd;
//			dinfo.comp_time = pstage->m_comp;
//			dinfo.comp_size = pstage->m_comp_size;
//			dinfo.comp_sd = pstage->m_comp_sd;
//			dinfo.findmate = pstage->m_findmate;
//		}
//		else { //else assign female info
//			dinfo.S0 = pstage->f_S0;
//			dinfo.alphaS = pstage->f_alphaS;
//			dinfo.betaS = pstage->f_betaS;
//			dinfo.sett_sd = pstage->f_sett_sd;
//			dinfo.pld = pstage->f_pld;
//			dinfo.pld_sd = pstage->f_pld_sd;
//			dinfo.comp_time = pstage->f_comp;
//			dinfo.comp_size = pstage->f_comp_size;
//			dinfo.comp_sd = pstage->f_comp_sd;
//			dinfo.findmate = pstage->f_findmate;
//		}
//		dinfo.S0_evo = pstage->S0_evo;
//		dinfo.S0_mut = pstage->S0_mut;
//		dinfo.S0_mut_size = pstage->S0_mut_size;
//		dinfo.alphaS_mut_size = pstage->alphaS_mut_size;
//		dinfo.betaS_mut_size = pstage->betaS_mut_size;
//
//
//	}
//	else { //if transfer is not stage dependent
//		if (pop_sinfo->sexdep == 1 && sex == 0) { //but it is sex dependent, assign male info
//			dinfo.S0 = pop_sinfo->m_S0;
//			dinfo.alphaS = pop_sinfo->m_alphaS;
//			dinfo.betaS = pop_sinfo->m_betaS;
//			dinfo.sett_sd = pop_sinfo->m_sett_sd;
//			dinfo.pld = pop_sinfo->m_pld;
//			dinfo.pld = pop_sinfo->m_pld_sd;
//			dinfo.comp_time = pop_sinfo->m_comp;
//			dinfo.comp_size = pop_sinfo->m_comp_size;
//			dinfo.comp_sd = pop_sinfo->m_comp_sd;
//			dinfo.findmate = pop_sinfo->m_findmate;
//		}
//		else { //else assign female info
//			dinfo.S0 = pop_sinfo->f_S0;
//			dinfo.alphaS = pop_sinfo->f_alphaS;
//			dinfo.betaS = pop_sinfo->f_betaS;
//			dinfo.sett_sd = pop_sinfo->f_sett_sd;
//			dinfo.pld = pop_sinfo->f_pld;
//			dinfo.pld_sd = pop_sinfo->f_pld_sd;
//			dinfo.comp_time = pop_sinfo->f_comp; //was already converted to hours in population set_settle info
//			dinfo.comp_size = pop_sinfo->f_comp_size;
//			dinfo.comp_sd = pop_sinfo->f_comp_sd;
//			dinfo.findmate = pop_sinfo->f_findmate;
//		}
//		dinfo.S0_evo = pop_sinfo->S0_evo;
//		dinfo.S0_mut = pop_sinfo->S0_mut;
//		dinfo.S0_mut_size = pop_sinfo->S0_mut_size;
//		dinfo.alphaS_mut_size = pop_sinfo->alphaS_mut_size;
//		dinfo.betaS_mut_size = pop_sinfo->betaS_mut_size;
//
//	}
//	//interindividual variation is set separately, so need to check for it separately
//
//	if (dinfo.sett_sd != -9) {
//		normal_distribution<float> S0_dist(dinfo.S0, dinfo.sett_sd);
//		normal_distribution<float> alpha_dist(dinfo.alphaS, dinfo.sett_sd);
//		normal_distribution<float> beta_dist(dinfo.betaS, dinfo.sett_sd);
//		dinfo.S0 = S0_dist(eng);
//		dinfo.alphaS = alpha_dist(eng);
//		dinfo.betaS = beta_dist(eng);
//	}
//	if (dinfo.pld_sd != -9) {
//		normal_distribution<float> pld_dist(dinfo.pld, dinfo.pld_sd);
//		dinfo.pld = pld_dist(eng);
//	}
//	if (dinfo.comp_sd != -9) {
//		normal_distribution<float> comp_dist;
//		if (mod->size_or_time == 0) {
//			comp_dist = normal_distribution<float>(dinfo.comp_size, dinfo.comp_sd);
//			dinfo.comp_size = comp_dist(eng);
//		}
//		else {
//			comp_dist = normal_distribution<float>(dinfo.comp_time, dinfo.comp_sd);
//			dinfo.comp_time = comp_dist(eng);
//		}
//	}
//}

bool Individual::check_settle(Cell* cell_sett, float pot_x, float pot_y, float pot_z) { //this checks the possible cell/patch to settle in before moving there
	//sample a random number:
	//cout << "checking cell for settlement, indiv status is " << status << endl;
	float ran_sett = unif_dist(eng);

	if (status == 9) { //forced settlement for a passive disperser. with some probability, the individual will settle, otherwise it continues (jsut to have some stochasticity)
		if (mod->patchmodel == 1 && cell_sett->patch == nullptr) { //cell might be suitable but isn't part of a patch so indivs cant settle there
			sett = false;
			return false;
		}
		if (dinfo.status == 8) {
			float s;
			if (pstage->s_sett.densdep) { //if settlement is density dependent, adjust settlement so it doesn't go above K?
				float n, k;
				if (mod->patchmodel == 1) { n = cell_sett->patch->p_subpop->size; k = cell_sett->patch->indepth_info.max_inds; }
				else { n = cell_sett->c_subpop->size; k = cell_sett->max_inds; }
				if (mod->stagestruct == 0) {

					s = dinfo.S0 / (1 + exp(-1 * (n / k - dinfo.betaS)*dinfo.alphaS));
				}
				else {
					//cout << "S0: " << dinfo.S0 << ", alpha S: " << dinfo.alphaS << ", betaS: " << dinfo.betaS << ", ";
					s = dinfo.S0 / (1 + exp(-1 * ((1 / k)*n - dinfo.betaS)*dinfo.alphaS)); //because it's stage structured, i am describing strength of density dependence (1/b) so i take 1/maxinds
				}

				if (ran_sett < s) { sett = true; }
				else {
					sett = false;
					status = 1;
					return false;
				}
			}
			else {
				sett = true;
			}
			 
		}
		else if (dinfo.status != 8 && ran_sett < 0.5) { sett = true; }//50% of the time, forced settlement will result in settlement
		else {
			sett = false;
			status = 1;
			return false;
		}
		if (sett == true) { 
			
			//cout << "forced settlement" << endl;
			//current_cell = cell_sett; //update current cell
			float new_x=-9, new_y=-9, new_z=-9;
			
			//float t, p_x, p_y, p_z;
			//bool north = false, south = false, east = false, west = false, up = false, check = true;
			//float diff_x, diff_y, diff_z, diff_r, diff_c, diff_l;
			//diff_r = cell_sett->r - current_cell->r;
			//diff_c = cell_sett->c - current_cell->c;
			//diff_l = current_cell->l - cell_sett->l;
			//if (diff_r > 0 && diff_c == 0 && diff_l == 0) { //north
			//	north = true;
			//}
			//else if (diff_r < 0 && diff_c == 0 && diff_l == 0) {//south
			//	south = true;
			//	//cout << "south, plane it intersects is y=" << cell_sett->cell_corners.lrc[1] << "< ";
			//}
			//else if (diff_r == 0 && diff_c > 0 && diff_l == 0) {//west
			//	west = true;
			//	
			//}
			//else if (diff_r == 0 && diff_c < 0 && diff_l == 0) { //east
			//	east = true;
			//	
			//}
			//else if (diff_r == 0 && diff_c == 0 && diff_l > 0) { //above
			//	up = true;
			//	
			//}
			//else { //not one of the cardinal directions or from above, so use east west and check coordinates
			//	float temp_phi = atan2(current_pos.y - cell_sett->midpoint[1], current_pos.x - cell_sett->midpoint[0]);
			//	temp_phi = fmod(temp_phi + M2_PI, M2_PI);
			//	if (temp_phi >= (4 * M_PI / 8) && temp_phi < (12 * M_PI / 8)) { //north west to southwest
			//	//the plane it intersects is x=min_x;
			//		west = true;
			//		//cout << "west, plane it intersects is x=" << cell_sett->cell_corners.llc[0] << "< ";
			//		t = (cell_sett->cell_corners.llc[0] - current_pos.x) / (pot_x - current_pos.x); //
			//		
			//	}
			//	else if (abs(temp_phi) > (12 * M_PI / 8) || temp_phi < (4 * M_PI / 8)) { //east
			//		//the plane it is intersection is x=max_x;
			//		east = true;
			//		//cout << "east, plane it intersects is x=" << cell_sett->cell_corners.lrc[0] << "< ";
			//		t = (cell_sett->cell_corners.lrc[0] - current_pos.x) / (pot_x - current_pos.x);
			//	}
			//	p_x = current_pos.x + t * (pot_x - current_pos.x);
			//	p_y = current_pos.y + t * (pot_y - current_pos.y);
			//	p_z = current_pos.z - t * (pot_z + -1 * current_pos.z);
			//	if (p_y > cell_sett->cell_corners.ulc[1]) { //if point of intersection is to the north
			//		north = true; west = false; east = false;
			//	}
			//	else if (p_y > cell_sett->cell_corners.llc[1]) {
			//		south = true; west = false; east = false;
			//	}
			//	if (p_z < cell_sett->min_depth) {
			//		up = true; east = false; west = false; north = false; south = false;
			//	}

			//}
			//if (north == true) {
			//	//cout << "north, plane it intersects is y=" << cell_sett->cell_corners.urc[1] << "< ";
			//	t = (cell_sett->cell_corners.ulc[1] - current_pos.y) / (pot_y - current_pos.y);
			//	if (lands->raster_3d[cell_sett->r-1][cell_sett->c][cell_sett->l]->habtype != 0) {
			//		cout << "should settle north but this is not a water cell!";
			//		check = false;
			//	}
			//}
			//else if (south == true) {
			//	//cout << "south, plane it intersects is y=" << cell_sett->cell_corners.lrc[1] << "< ";
			//	t = (cell_sett->cell_corners.llc[1] - current_pos.y) / (pot_y - current_pos.y);
			//	if (lands->raster_3d[cell_sett->r+1][cell_sett->c][cell_sett->l]->habtype != 0) {
			//		cout << "should settle south but this is not a water cell!";
			//		check = false;
			//	}
			//}
			//else if (west == true) {
			//	//cout << "west, plane it intersects is x=" << cell_sett->cell_corners.llc[0] << "< ";
			//	t = (cell_sett->cell_corners.llc[0] - current_pos.x) / (pot_x - current_pos.x); //
			//	if (lands->raster_3d[cell_sett->r][cell_sett->c - 1][cell_sett->l]->habtype != 0) {
			//		cout << "should settle west but this is not a water cell!";
			//		check = false;
			//	}
			//}
			//else if (east == true) {
			//	//cout << "east, plane it intersects is x=" << cell_sett->cell_corners.lrc[0] << "< ";
			//	t = (cell_sett->cell_corners.lrc[0] - current_pos.x) / (pot_x - current_pos.x);
			//	if (lands->raster_3d[cell_sett->r][cell_sett->c + 1][cell_sett->l]->habtype != 0) {
			//		cout << "should settle east but this is not a water cell!";
			//		check = false;
			//	}
			//}
			//else if (up == true) {
			//	//cout << "up, plane it intersects is z=" << cell_sett->min_depth << "< ";
			//	t = (current_pos.z - cell_sett->min_depth) / (current_pos.z - pot_z);
			//	if (lands->raster_3d[cell_sett->r][cell_sett->c ][cell_sett->l-1]->habtype != 0) {
			//		cout << "should settle up but this is not a water cell!";
			//		check = false;
			//	}
			//}
			//if (check == false) {
			//	cout << "after calcs, it would have settled on the wrong face, dead" << endl;
			//	dead = true;
			//	status = 4;
			//	sett = false;
			//	return false;
			//}
			//
			//p_x = current_pos.x + t * (pot_x - current_pos.x);
			//p_y = current_pos.y + t * (pot_y - current_pos.y);
			//p_z = current_pos.z - t * (pot_z + -1 * current_pos.z);
			////cout << " the point of intersection is " << p_x << ", " << p_y << ", " << p_z;

			//new_x = p_x; 
			//new_y = p_y; 
			//new_z = p_z;

			//cout << endl;
			
			
			//
			////calculate the edge of the cell for settlement
			//if (current_pos.y <= current_cell->cell_corners.ulc[1] && current_pos.y >= current_cell->cell_corners.llc[1]) { //if indiv is coming at the cell from the side (ie between y boundaries)
			//	if (current_pos.x <= current_cell->cell_corners.lrc[0] && current_pos.x >= current_cell->cell_corners.llc[0] && current_pos.z <= current_cell->min_depth) { //coming at it from above
			//		new_x = pot_x; new_y = pot_y; new_z = current_cell->min_depth; //make it settle on top of the cell
			//		//cout << "top ";
			//	}
			//	else if (current_pos.x < current_cell->cell_corners.llc[0]) { //coming at it from the left
			//		new_x = current_cell->cell_corners.llc[0]; //make new coordinate the minimum x
			//		new_y = pot_y; new_z = pot_z;
			//		//cout << "left ";
			//	}
			//	else if (current_pos.x > current_cell->cell_corners.lrc[0]) { //coming at it from the right
			//		new_x = current_cell->cell_corners.lrc[0]; //make new coordinate maximum x
			//		new_y = pot_y; new_z = pot_z;
			//		//cout << "right ";
			//	}
			//}
			//else if (current_pos.y > current_cell->cell_corners.ulc[1]) { //coming at the cell from north
			//	new_y = current_cell->cell_corners.ulc[1];
			//	new_x = pot_x; new_z = pot_z;
			//	//cout << "north ";
			//}
			//else if (current_pos.y < current_cell->cell_corners.llc[1]) { //coming at the cell from south
			//	new_y = current_cell->cell_corners.llc[1];
			//	new_x = pot_x; new_z = pot_z;
			//	//cout << "south ";
			//}
			//if (new_x==-9 || new_y==-9 || new_z==-9) {
			//	cout << "in check_settle, coordinate conditions not met, i've missed something." << endl;
			//	dead = true;
			//	status = 4;
			//	sett = false;
			//	return false;
			//}
			//cout << endl;

		
			//randomly assign a position in the cell that is sure to have an open water neighbour
			if (cell_sett->possibilities.size() == 0) { //ow neighbours haven't been assessed
				cell_sett->find_ow_neigh(lands);
			}

			//shuffle the possibilities so that when i pick one it is random
			random_shuffle(cell_sett->possibilities.begin(), cell_sett->possibilities.end());
			if (cell_sett->possibilities.size() == 0) { cout << "there were no open water neighbours..." << endl; }
			Cell::mat_index ind = cell_sett->possibilities[0];
			if (lands->raster_3d[ind.row][ind.col][ind.layer] == nullptr) {
				cout << "one of the possibilities in init_pos is a nullptr somehow!" << endl;
			}

			uniform_real_distribution<float> x_dist(cell_sett->cell_corners.llc[0] + 1, cell_sett->cell_corners.lrc[0] - 1); //dont allow it to be given the actual min and max
			uniform_real_distribution<float> y_dist(cell_sett->cell_corners.llc[1] + 1, cell_sett->cell_corners.ulc[1] - 1);
			uniform_real_distribution<float> z_dist(cell_sett->min_depth + 1, cell_sett->max_depth - 1);
			float x_coord;
			float y_coord;
			float z_coord;

			if (ind.row < cell_sett->r) { //north
				y_coord = cell_sett->cell_corners.ulc[1]; //maximum y coordinate
				x_coord = x_dist(eng); //random x coordinate
				z_coord = z_dist(eng); //random z coordinate
			}
			else if (ind.row > cell_sett->r) { //south
				y_coord = cell_sett->cell_corners.llc[1]; //minimum y coordinate
				x_coord = x_dist(eng);
				z_coord = z_dist(eng);
			}
			else { //same row
				if (ind.col < cell_sett->c) {//west
					y_coord = y_dist(eng);
					x_coord = cell_sett->cell_corners.ulc[0]; //minimum x coordinate
					z_coord = z_dist(eng); //random z coordinaet
				}
				else if (ind.col > cell_sett->c) { //east
					y_coord = y_dist(eng);
					x_coord = cell_sett->cell_corners.urc[0]; //maximum x coordinate
					z_coord = z_dist(eng); //random z coordinaet
				}
				else if (ind.layer < cell_sett->l) { //up
					y_coord = y_dist(eng);
					x_coord = x_dist(eng);
					z_coord = cell_sett->min_depth; //start on top of the cell
				}
			}


		//current_pos.x = x_coord; current_pos.y = y_coord; current_pos.z = z_coord;
			new_x = x_coord; new_y = y_coord; new_z = z_coord;

			dinfo.per_step_disp_dist += calc_distance(pot_x, pot_y, pot_z); //still use predicted so that i dont overestimate distance due to random position
			if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << "status 9 in check settle, per step is 0" << endl; }
			/*current_pos.x = new_x; current_pos.y = new_y; current_pos.z = new_z;
			m_x_movements.push_back(new_x); m_y_movements.push_back(new_y); m_z_movements.push_back(new_z);*/
			update_pos(new_x, new_y, new_z, cell_sett->r, cell_sett->c, cell_sett->l, true, false, false, true);
			//cout << "settling at " << new_x << ", " << new_y << ", " << new_z << endl;
			status = 9;
			if (dinfo.status == 8) { status = 3; } //if it was here because of buffer capture, give it the new status of successful settlement
			//current_patch = cell_sett->patch;
			current_cell = cell_sett;
			//current_patch = current_cell->patch;
			//if (current_patch == nullptr) { cout << "settled with status 9 but current patch is null..." << endl; }
			return true;
		}
		
	}

	//met habitat but not out of comp period yet
	if (dinfo.comp == false) { //if it has a competency period, check that it has dispersed long enough
		sett = false;
		status = 1;
		//cout << "outside of comp time" << endl;
		return false;
	}

	//met habitat but it's the natal patch or it has tried to settle there before (for active dispersers only, since passive wouldnt have a choice)
	if (dinfo.phase == 1) {
		if (mod->patchmodel == 0) {
			if (cell_sett->cell_ID == natal_cell) {
				//cout << "tried to settle in natal cell" << endl;
				sett = false;
				return false;
			}
			else if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), cell_sett->cell_ID)) {
				//cout << "already tried to settle here" << endl;
				sett = false;
				return false;
			}
		}
		else { //patch based
			if (cell_sett->patch == nullptr) { //cell might be suitable but isn't part of a patch so indivs cant settle there
				sett = false;
				return false;
			}
			//if (cell_sett->patch !=nullptr && cell_sett->patch->patch_ID == natal_patch) {
			//	//cout << "tried to settle in natal patch" << endl;
			//	sett = false;
			//	return false;
			//}
			if (cell_sett->patch != nullptr && count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), cell_sett->patch->patch_ID)) {
				//cout << "already tried to settle here" << endl;
				sett = false;
				return false;
			}
		}
	}

	//up to this point, passive dispersers with forced settlement, any indivs still within their comp time, and active disp
	//who have already tried to settle here previously, have been taken care of

	//check settlement conditions
	//if habitat suitability is the only requirement in a cell based model
	if (mod->patchmodel == 0 && pstage->s_sett.densdep == 0 && dinfo.findmate == 0 && cell_sett->K > 0) { 
		if (ran_sett < dinfo.S0) {
			status = 3;
			//cout << "settled" << endl;
			sett= true;
		}
		else { //if it doesn't settle
			//cout << "decided not to settle" << endl;			
		}
	}
	//if it's a patch-based model, habitat suitability is the only requirement
	else if (mod->patchmodel == 1 && pstage->s_sett.densdep == 0 && dinfo.findmate == 0 && cell_sett->patch->patch_info.K > 0) {
		if (ran_sett < dinfo.S0) {
			status = 3; sett= true;
		}
	}

	float s=-9; bool foundmate=false;

	if (pstage->s_sett.densdep == 1 && dinfo.findmate==0) { //if settlement is density dependent but no mating requirements
		float n, k;
		if (mod->patchmodel == 1) { n = cell_sett->patch->p_subpop->size; k = cell_sett->patch->indepth_info.max_inds; }
		else { n = cell_sett->c_subpop->size; k = cell_sett->max_inds;}
		if (mod->stagestruct == 0) {
			
			s = dinfo.S0 / (1 + exp(-1 * (n / k - dinfo.betaS)*dinfo.alphaS));
		}
		else {
			//cout << "S0: " << dinfo.S0 << ", alpha S: " << dinfo.alphaS << ", betaS: " << dinfo.betaS << ", ";
			s = dinfo.S0 / (1 + exp(-1 * ((1/k)*n - dinfo.betaS)*dinfo.alphaS)); //because it's stage structured, i am describing strength of density dependence (1/b) so i take 1/maxinds
		}

		if (ran_sett < s) { status = 3; sett= true; }
		
	}

	else if (pstage->s_sett.densdep == 0 && dinfo.findmate == 1) { //if mating requirements
	
		if (mod->patchmodel == 1) { //patch based model
			if (sex == 1) { //female
				if (cell_sett->patch->p_subpop->nmales > 0) { foundmate = true; }
				else { foundmate = false; }
			}
			else { //male
				if (cell_sett->patch->p_subpop->nfemales > 0) { foundmate = true; }
				else { foundmate = false; }
			}
		}
		else { //cell-based
			if (sex == 1) { //female
				if (cell_sett->c_subpop->nmales > 0) { foundmate = true; }
				else { foundmate = false; }
			}
			else { //male
				if (cell_sett->c_subpop->nfemales > 0) { foundmate = true; }
				else { foundmate = false; }
			}
		}
		//using S0 as a maximum settlement probability when density dependence is not in effect
		if (foundmate && ran_sett < dinfo.S0) { status = 3; sett= true; } 

	
	}

	else if (pstage->s_sett.densdep == 1 && dinfo.findmate == 1) { //if both are needed
		if (ran_sett < s && foundmate) { status = 3; sett= true; }
	}


	
	float new_x, new_y, new_z;

	if (sett == true) {
		//cout << "decided to settle!" << endl;
		current_cell = cell_sett; //update current cell
		//if (mod->patchmodel == 1) { current_patch = current_cell->patch; if (current_patch == nullptr) { cout << "after active check settle, tried assigning patch but it was null!" << endl; } }
		//calculate the edge of the cell for settlement
		//if (current_pos.y <= current_cell->cell_corners.ulc[1] && current_pos.y >= current_cell->cell_corners.llc[1]) { //if indiv is coming at the cell from the side (ie between y boundaries)
		//	if (current_pos.x <= current_cell->cell_corners.lrc[0] && current_pos.x >= current_cell->cell_corners.llc[0] && current_pos.z <= current_cell->min_depth) { //coming at it from above
		//		new_x = pot_x; new_y = pot_y; new_z = current_cell->min_depth; //make it settle on top of the cell
		//	}
		//	else if (current_pos.x < current_cell->cell_corners.llc[0]) { //coming at it from the left
		//		new_x = current_cell->cell_corners.llc[0]; //make new coordinate the minimum x
		//		new_y = pot_y; new_z = pot_z;
		//	}
		//	else { //coming at it from the right
		//		new_x = current_cell->cell_corners.lrc[0]; //make new coordinate maximum x
		//		new_y = pot_y; new_z = pot_z;
		//	}
		//}
		//else if (current_pos.y > current_cell->cell_corners.ulc[1]) { //coming at the cell from north
		//	new_y = current_cell->cell_corners.ulc[1];
		//	new_x = pot_x; new_z = pot_z;
		//}
		//else if (current_pos.y < current_cell->cell_corners.llc[1]) { //coming at the cell from south
		//	new_y = current_cell->cell_corners.llc[1];
		//	new_x = pot_x; new_z = pot_z;
		//}	

		//randomly assign a position in the cell that is sure to have an open water neighbour
		if (cell_sett->possibilities.size() == 0) { //ow neighbours haven't been assessed
			cell_sett->find_ow_neigh(lands);
		}

		//shuffle the possibilities so that when i pick one it is random
		random_shuffle(cell_sett->possibilities.begin(), cell_sett->possibilities.end());
		if (cell_sett->possibilities.size() == 0) { cout << "there were no open water neighbours..." << endl; }
		Cell::mat_index ind = cell_sett->possibilities[0];
		if (lands->raster_3d[ind.row][ind.col][ind.layer] == nullptr) {
			cout << "one of the possibilities in init_pos is a nullptr somehow!" << endl;
		}

		uniform_real_distribution<float> x_dist(cell_sett->cell_corners.llc[0] + 1, cell_sett->cell_corners.lrc[0] - 1); //dont allow it to be given the actual min and max
		uniform_real_distribution<float> y_dist(cell_sett->cell_corners.llc[1] + 1, cell_sett->cell_corners.ulc[1] - 1);
		uniform_real_distribution<float> z_dist(cell_sett->min_depth + 1, cell_sett->max_depth - 1);
		float x_coord;
		float y_coord;
		float z_coord;

		if (ind.row < cell_sett->r) { //north
			y_coord = cell_sett->cell_corners.ulc[1]; //maximum y coordinate
			x_coord = x_dist(eng); //random x coordinate
			z_coord = z_dist(eng); //random z coordinate
		}
		else if (ind.row > cell_sett->r) { //south
			y_coord = cell_sett->cell_corners.llc[1]; //minimum y coordinate
			x_coord = x_dist(eng);
			z_coord = z_dist(eng);
		}
		else { //same row
			if (ind.col < cell_sett->c) {//west
				y_coord = y_dist(eng);
				x_coord = cell_sett->cell_corners.ulc[0]; //minimum x coordinate
				z_coord = z_dist(eng); //random z coordinaet
			}
			else if (ind.col > cell_sett->c) { //east
				y_coord = y_dist(eng);
				x_coord = cell_sett->cell_corners.urc[0]; //maximum x coordinate
				z_coord = z_dist(eng); //random z coordinaet
			}
			else if (ind.layer < cell_sett->l) { //up
				y_coord = y_dist(eng);
				x_coord = x_dist(eng);
				z_coord = cell_sett->min_depth; //start on top of the cell
			}
		}


		//current_pos.x = x_coord; current_pos.y = y_coord; current_pos.z = z_coord;
		new_x = x_coord; new_y = y_coord; new_z = z_coord;
		//vector<int> indexes=calc_r_c_l(new_x, new_y, new_z);
		dinfo.per_step_disp_dist += calc_distance(new_x, new_y, new_z);
		if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " in check settle, per step is 0" << endl; }
		update_pos(new_x, new_y, new_z, cell_sett->r, cell_sett->c, cell_sett->l, true, false, false, true);
		//cout << "settling at " << new_x << ", " << new_y << ", " << new_z << endl;
		//current_pos.x = new_x; current_pos.y = new_y; current_pos.z = new_z;
		//m_x_movements.push_back(new_x); m_y_movements.push_back(new_y); m_z_movements.push_back(new_z);
		status = 3;
		return true;
	}
	else {
		//cout << "indiv " << ID << " decided not to settle";
		if (dinfo.phase == 1) {//active disperse
			//add it to the list of cells/patches that it has already tested
			if (mod->patchmodel == 1) {
				dinfo.tested_sett.push_back(cell_sett->patch->patch_ID);
				//cout << "added "<< cell_sett->patch->patch_ID << " to tested_sett" << endl;
				//if (dinfo.goal_bias == 1) { cout << ", it had goal bias"; }
				if (dinfo.goal_bias == 1 && dinfo.goal_focus == cell_sett) { //if it was locked onto a goal and this was it and it has tested it and decided not to settle, reset the goal bias
					//cout << ", and it was the cell that was tested for settlement" << endl;
					dinfo.goal_bias = 0; dinfo.goal_focus = nullptr;
				}
				else if(dinfo.goal_bias!=1 || dinfo.goal_focus !=cell_sett){
					//cout << ", it wasn't the cell tested for settlement, sett prob was " << s << " and indepth info for patch maxinds was "<< cell_sett->patch->indepth_info.max_inds << endl;
				}

			}
			else {
				dinfo.tested_sett.push_back(cell_sett->cell_ID);
			}
			if (dinfo.goal_bias == 1 && dinfo.goal_focus == cell_sett) { //if it was locked onto a goal and this was it and it has tested it and decided not to settle, reset the goal bias
				dinfo.goal_bias = 0; dinfo.goal_focus = nullptr;
			}
		}
		status = 1;
		return false;
	}
}

//movement model equations
float Individual::calc_cell_mag(int row, int col, int layer, float x, float y, float z) {
	//cout << "in calc_cell_mag, current pos is " << current_pos.x << ", " << current_pos.y << ", " << current_pos.z << endl;
	float cell_mag = 0;
	//to come up with a rule about cell_mag
	float res = lands->get_land_att(7);
	float dint = lands->get_land_att(8);
	float propr = dint / res;

	//if (layer - current_cell->l != 0) { //if the movement will require to change layers
	//	cell_mag = dint;
	//}
	//else {
	//	cell_mag = res;
	//}
	if (row - current_cell->r == 0 && col - current_cell->c == 0 && layer - current_cell->l != 0) { //if all that is changing is the layer
		cell_mag = dint; 
		return cell_mag;
	}




	//
	////to figure out which side of the cell i am moving towards, i need to use the cell's midpoint
	////float phi = atan2(y - current_pos.y, x - current_pos.x);
	//float phi = atan2(y - current_cell->midpoint[1], x - current_cell->midpoint[0]);

	////for this section, i am solving for a y=mx+b line using my current position, to then solve for the other coordinate
	//float m = (y - current_pos.y) / (x - current_pos.x); //calculate the slope of the line between the goal and current_pos
	//float b = current_pos.y - m * current_pos.x; //calculate the intercept of the y=mx+b
	//float x_val = 0, y_val = 0, z_val;
	//if (fabs(phi) > (6 * M_PI / 8)) { //the entire left column of cells will have the same x coordinate 
	//	x_val = current_cell->cell_corners.llc[0]-1; //move a bit further than cell boundary
	//	y_val = m * x_val + b; //simple y=mx+b formula
	//}
	//else if (fabs(phi) >= (2 * M_PI / 8) && fabs(phi) <= (6 * M_PI / 8)) { //row of cells above and below current position will all have the same y coordinate
	//	if (phi > 0) { // above
	//		y_val = current_cell->cell_corners.ulc[1] + 1; //move just beyond the upper y coordinate

	//	}
	//	else { //below
	//		y_val = current_cell->cell_corners.llc[1] - 1; //move just below the lower y coordinate
	//	}
	//	x_val = (y_val - b) / m; //y=mx+b solved for x

	//}
	//else if (fabs(phi) < (2 * M_PI / 8)) { //column of cells to the right of current pos
	//	x_val = current_cell->cell_corners.lrc[0] + 1; //move a bit further than cell boundary
	//	y_val = m * x_val + b; //simple y=mx+b formula

	//}
	//cell_mag = sqrt(pow(x_val - current_pos.x, 2) + pow(y_val - current_pos.y, 2));

	//if (layer - current_cell->l == 0) { //same layer, i want to multiply z by the lowest number so indiv isnt catapulted into a new layer, so if dint is smaller, use dint
	//	if (dint < cell_mag) { cell_mag = dint; }
	//}
	return cell_mag;
}

void Individual::update_pos(float new_x, float new_y, float new_z, int row, int col, int layer, bool mem, bool mem_empty,  bool c_cell, bool c_pos) {
	

	if (mem_empty == true) { 
		for (int x = 0; x < memory.size() + 1; x++) {
			memory.pop();
		}
		memory.push(current_pos); //i want to remember the movement from where i am to this next position, so i delete everything but my current position
	} 
	if (c_pos == true) { //update current pos
		current_pos.x = new_x; current_pos.y = new_y; current_pos.z = new_z; 
		m_x_movements.push_back(new_x); m_y_movements.push_back(new_y); m_z_movements.push_back(new_z);
	}
	if (mem == true) {
		locn temp; temp.x = new_x; temp.y = new_y; temp.z = new_z;
		if (memory.size() >= dinfo.memory) { memory.pop(); }
		memory.push(temp);
	}
	if (c_cell == true) {
		current_cell = lands->raster_3d[row][col][layer];
		if (sett == true && mod->patchmodel == 1) { 
			current_patch = current_cell->patch;
			if (current_cell->patch == nullptr) { cout << "in update pos, assigned patch but it was null!" << endl; }
			else if (current_cell->patch->p_subpop == nullptr) { cout << "in update pos, assigned patch but subpop pointer was null " << endl; }

		}
	}
}

vector<int> Individual::calc_r_c_l(float x, float y, float z){
	int row, col, layer;
	vector<int> indexes;

	col = floor((x - lands->get_land_att(9)) / lands->get_land_att(7)); //current position-minimum x value in the landscape / resolution will give you the column it's in
	
	//y is tricky because it goes from max_y at the top left
	if (y < lands->get_land_att(10)) { //if y is less than min_y (meaning, row is > nrows)
		row = floor(abs((y - lands->get_land_att(10) )/ lands->get_land_att(7)))+ lands->get_land_att(1); //add nrows because ymin is bottom left corner!
	}
	else { //otherwise
		row = floor((lands->get_land_att(12) - y) / lands->get_land_att(7)); //ymax-y will always give me a positive number for inside the landscape and a neg for outside, which is what i want
	}

	layer = floor((z - lands->get_land_att(13)) / lands->get_land_att(8)); //current z - minimum z / interval between layers
	
	//if (row == -1 || col == -1 || layer == -1) { indexes.push_back(row); indexes.push_back(col); indexes.push_back(layer); }
	//else{ indexes.push_back(row); indexes.push_back(col); indexes.push_back(layer); }
	indexes.push_back(row); indexes.push_back(col); indexes.push_back(layer);
	return indexes;
}

float Individual::calc_distance(float x, float y, float z) {
	float dist= sqrt(pow(x - current_pos.x, 2) + pow(y - current_pos.y, 2) + pow(z - current_pos.z, 2));
	return dist;
}


vector<float> Individual::evade(float phi) { //this function is used when unsuitable habitat is encountered and individual needs to move away from/around it
	//this function is called during dispersal
	//the point of this function is to be able to move away from unsuitable habitat and to track changes in depth
	//cout << "at the start of evade, the current cell is " << current_cell->r << ", " << current_cell->c << "," << current_cell->l << endl;
	//since this is no longer applicable to leaving the natal patch, going down isn't important anymore. same layer or up
	if (current_cell->cell_speed == -9999) { cout << "called evade and current cell speed is already -9999" << endl; }
	float to_go;//this is how far indiv is still allowed to travel
	if (dinfo.SL > 1) { to_go = (dinfo.rho*current_cell->cell_speed + (1-dinfo.rho)*dinfo.SL) - dinfo.per_step_disp_dist; } //if it's actively dispersing, it can go until its step length
	else { //if it's passive, it can go until it reaches the distance the current would have carried it
		to_go = current_cell->cell_speed - dinfo.per_step_disp_dist;
	}
	//cout << "in evade, can still travel " << to_go << "m distance" << endl;
	float new_x = current_pos.x, new_y = current_pos.y, new_z = current_pos.z;
	vector<float> coords;
	bool valid = false;
	//uniform_int_distribution<int> dir_unif;
	//if(current_cell->l==0) {dir_unif= uniform_int_distribution<int>(0, 1); } //if it's in the top layer, it can't go up! only left/right
	//else { dir_unif = uniform_int_distribution<int>(0, 2); } //otherwise it can go any direction
	//dir_unif = uniform_int_distribution<int>(0, 1);
	//float theta = acos(0); //because changing depth isnt soemthing i want to do in this function
	//uniform_real_distribution<double> turn_ang{ 0,2 };
	//while (valid == false) {
	//	new_x = 0; new_y = 0; new_z = 0; //every calculation needs to have a clean slate
	//	float try_phi = phi;
	//	try_phi += turn_ang(eng);

	//	//int index = dir_unif(eng);
	//	//if (index == 0 || index == 1) {
	//	//	if (index == 0) { //left
	//	//		cout << "moving left" << endl;
	//	//		phi += acos(0);
	//	//	}
	//	//	else if (index == 1) { //right
	//	//		cout << "moving right " << endl;
	//	//		phi -= acos(0);
	//	//	}
	//	//}
	//
	//	//else { //up, so move one layer up, then dont change depth, keep phi how it was
	//	//	cout << " moving up" << endl;
	//	//	float dist_to_layer = current_pos.z - ((current_cell->l - 1) * lands->get_land_att(8) + lands->get_land_att(8)); //this will make it move only as much as it needs to get to next layer
	//	//	new_z = dist_to_layer; //this will then be added to current pos
	//	//																													 //update position to have moved up one layer, update current cell, add to m_xyz and memory
	//	//	//update_pos(current_pos.x, current_pos.y, current_pos.z - dist_to_layer, current_cell->r, current_cell->c, current_cell->l - 1, true, false, true, true);
	//	//}
	//	new_x = current_pos.x + to_go * sin(theta)*cos(try_phi);
	//	new_y = current_pos.y + to_go * sin(theta)*sin(try_phi);
	//	new_z = current_pos.z - new_z; //(because new_z already has the dist_to layer, if indiv is moving up

	//	vector<int> indexes= calc_r_c_l(new_x, new_y, new_z);
	//	cout << "evade to " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << "at position " << new_x << ", " << new_y << ", " << new_z << endl;
	//	if (indexes[0] < 0 || indexes[0] >= lands->get_land_att(1) || indexes[1] < 0 || indexes[1] >= lands->get_land_att(2)) {
	//		valid = false;
	//	}
	//	else if (lands->raster_3d[indexes[0]][indexes[1]][current_cell->l] == nullptr) {
	//		valid = false;
	//	}
	//	else if (lands->raster_3d[indexes[0]][indexes[1]][current_cell->l]->K > 0) {
	//		valid = false;
	//	}
	//	else {
	//		valid = true;
	//	}

	//}

	vector<double> dirs{ -M_PI / 2, -M_PI/4, -3*M_PI/4, M_PI, M_PI/4, 3*M_PI/4, M_PI / 2 }; 

	//this will make the individual turn any direction except forward diagonal ( so i dont run into the problem of moving diagonally across an invalid cell and getting stuck)
	double try_phi, try_theta;
	//move clockwise?
	
	random_shuffle(dirs.begin(), dirs.end()); //so there is no bias which cell i test first
	int counter = 0;

	for (int c=0; c<dirs.size(); c++) {
		if (valid == true) { break; }
		try_phi = phi + dirs[counter]; //after it's been shuffled, take the first element
		
		//cout << "trying " << dirs[0] << endl;
		//if (c%2 == 0) {
		//	try_theta = 1.55; //go up a bit
		//}
		//else {
		//	try_theta = acos(0); //otherwise stay level
		//}
		
		/*uniform_real_distribution<double> up_dist{ 1.55, acos(0)};
		try_theta = up_dist(eng);*/
		try_theta = 1.57;
		new_x = current_pos.x + to_go * sin(try_theta)*cos(try_phi);
		new_y = current_pos.y + to_go * sin(try_theta)*sin(try_phi);
		new_z = current_pos.z - to_go * cos(try_theta); //(because new_z already has the dist_to layer, if indiv is moving up

		vector<int> indexes = calc_r_c_l(new_x, new_y, new_z);
		//cout << "evade to " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << "at position " << new_x << ", " << new_y << ", " << new_z << endl;
		if (indexes[0] < 0 || indexes[0] >= lands->get_land_att(1) || indexes[1] < 0 || indexes[1] >= lands->get_land_att(2)
			|| indexes[2] < 0 || indexes[2] >= lands->get_land_att(5)) {
			//cout << "outside " << endl;
			valid = false;
		}
		else if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]] == nullptr) {
			//cout << "was null" << endl;
			valid = false;
		}
		else if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]]->habtype != 0) {
			//cout << "was solid" << endl;
			valid = false;
		}
		else if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]]->cell_speed == -9999) {
			valid = false;
		}
		else {
			//cout << "valid" << endl;
			valid = true;
		}
		counter++;
		//if (counter >= dirs.size()) { //if this last round 

		//	cout << "in evade, none of the directions gave a valid cell, indiv's SL was " << dinfo.SL << endl; 
		//}
	}
	//while (valid == false) {
	//	try_phi = phi + dirs[counter]; //after it's been shuffled, take the first element
	//	//cout << "trying " << dirs[0] << endl;

	//	new_x = current_pos.x + to_go * sin(theta)*cos(try_phi);
	//	new_y = current_pos.y + to_go * sin(theta)*sin(try_phi);
	//	new_z = current_pos.z; //(because new_z already has the dist_to layer, if indiv is moving up

	//	vector<int> indexes = calc_r_c_l(new_x, new_y, new_z);
	//	//cout << "evade to " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << "at position " << new_x << ", " << new_y << ", " << new_z << endl;
	//	if (indexes[0] < 0 || indexes[0] >= lands->get_land_att(1) || indexes[1] < 0 || indexes[1] >= lands->get_land_att(2)) {
	//		//cout << "outside " << endl;
	//		valid = false;
	//	}
	//	else if (lands->raster_3d[indexes[0]][indexes[1]][current_cell->l] == nullptr) {
	//		//cout << "was null" << endl;
	//		valid = false;
	//	}
	//	else if (lands->raster_3d[indexes[0]][indexes[1]][current_cell->l]->habtype!= 0) {
	//		//cout << "was solid" << endl;
	//		valid = false;
	//	}
	//	else if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]]->cell_speed == -9999) {
	//		valid = false;
	//	}
	//	else {
	//		//cout << "valid" << endl;
	//		valid = true;
	//	}
	//	counter++;
	//}
	if (valid==false) { //if it cant move to the side, try moving up?
		//cout << "in evade, none of the directions gave a valid cell, trying upwards " << endl; 
		//cout << "current pos is " << current_pos.x << ", " << current_pos.y << ", " << current_pos.z << " and current cell has hab type " << current_cell->habtype << endl;
		//cout << "indiv has SL " << dinfo.SL << endl;
		
		if (current_cell->l != 0 && lands->raster_3d[current_cell->r][current_cell->c][current_cell->l - 1] != nullptr) {
			//cout << "hab type of cell above in layer " << current_cell->l-1 << " is " << lands->raster_3d[current_cell->r][current_cell->c][current_cell->l - 1]->habtype << endl;
			if (lands->raster_3d[current_cell->r][current_cell->c][current_cell->l - 1]->habtype == 0) { //if cell above is water
				valid = true;
				new_x = current_pos.x; new_y = current_pos.y; new_z = lands->raster_3d[current_cell->r][current_cell->c][current_cell->l - 1]->max_depth - 1;
			}
		}
		if(valid==false){ //if that didn't work
			//cout << "that didn't work, going back to previous position" << endl;
			new_x = m_x_movements[m_x_movements.size() - 1]; //go back to the last valid coordinate
			new_y = m_y_movements[m_y_movements.size() - 1];
			new_z = m_z_movements[m_z_movements.size() - 1];
		}
		//new_x = current_pos.x; new_y = current_pos.y; new_z = current_pos.z;
	}
	if (new_x == 0 || new_y == 0) { cout << "assigning zeros in evade" << endl; }

	vector<int> indexes = calc_r_c_l(new_x, new_y, new_z);
	//calculate the angle between current_pos and new cell
	if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]]->cell_speed == -9999) { cout << "did evade, even after conditions, it assigned a non water cell, to go was " << to_go  << endl; }
	
	coords.push_back(new_x); coords.push_back(new_y); coords.push_back(new_z);
	return coords;
}

bool Individual::check_cell(int row, int col, int layer, float x, float y, float z) {
	//if the step is outside the landscape xy boundaries
	bool test = true; 

	if (row< 0 || row >= lands->get_land_att(1) || col < 0 || col >= lands->get_land_att(2)) {
		dead = true; //individual has died
		status = 8; //indiv has been absorbed
		dinfo.per_step_disp_dist += calc_distance(x, y, z);
		update_pos(x, y, z, row, col, layer, false, false, false, true);
		//cout << "outside xy landscape" << row << ", " << col << ", " << layer << endl;
		//return false;
		test = false;
	}
	//if (layer < 0) { cout << "soemthing went wrong before check_cell, z hasn't been adjusted for surface" << endl; }
	else if (layer > lands->get_land_att(5)) { cout << "layer too deep for some reason, layer is " << layer << endl; }
	//if the cell is null
	else if (lands->raster_3d[row][col][layer] == nullptr) {
		dead = true; //individual has died
		status = 8; //indiv has been absorbed
		dinfo.per_step_disp_dist += calc_distance(x, y, z);
		update_pos(x, y, z, row, col, layer, false, false, false, true);
		//cout << "null cell" << row << ", " << col << ", " << layer << endl;
		//return false;
		test = false;
	}
	else if (dinfo.phase==0 && dinfo.buoy_max ==-9 && dinfo.buoy_min==-9){//only applies to passive individuals that have no buoyancy limits because these are treated as ABSORBING BOUNDARIES
		if (z > lands->get_land_att(14)) { //if indiv goes too deep
			dead = true; //individual has died
			status = 8; //indiv has been absorbed
			dinfo.per_step_disp_dist += calc_distance(x, y, z);
			update_pos(x, y, z, row, col, layer, false, false, false, true);
			//cout << "too deep" << row << ", " << col << ", " << layer << endl;
			//return false;
			test = false;
		}
		else if (z < lands->get_land_att(13) && lands->get_land_att(13) !=0) { //if the first layer in the model is not the surface
			dead = true; //passive indivs will still pass upwards of that limit
			status = 8;
			dinfo.per_step_disp_dist += calc_distance(x, y, z);
			update_pos(x, y, z, row, col, layer, false, false, false, true);
			//cout << "shallow" << row << ", " << col << ", " << layer << endl;
			//return false;
			test = false;
		}
	}
	/*else {
		return true;
	}*/
	return test;
	//note that this function doesn't deal with z boundaries (except for passive, where they die if they go too deep)
	//this function is made to see if the cell they are travelling into is valid, other movement rules are applied within traverse functions
}



void Individual::traverse2(float x, float y, float z) {
	int row, col, layer;
	vector<int> indexes;
	float max_dist;
	bool valid;
	if (dinfo.phase == 0) {
		max_dist = current_cell->cell_speed;
	}
	else {
		max_dist = calc_distance(x, y, z);
	}
	if (dead == true) {
		return;
	}
	//calculate row, col and layer of destination
	indexes = calc_r_c_l(x, y, z);
	row = indexes[0]; col = indexes[1]; layer = indexes[2];
	//cout << "indiv " << ID;
	//valid = check_cell(row, col, layer, x, y, z);
	//cout << " adjusting" << endl;
	//if from the very beginning the movement is within the same cell, 
	if (current_cell->r == row && current_cell->c == col && current_cell->l == layer) {
		//cout << "same cell" << endl;
		dinfo.per_step_disp_dist += calc_distance(x, y, z);
		update_pos(x, y, z, row, col, layer, true, false, false, true); //update current_pos, memory, not current_cell and dont empty memory
		//dont need to update cell because it hasnt left the cell
		//if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " same cell, per step was 0, max dist was " << max_dist << endl; }
		//cout << "perstep disp dist is now " << dinfo.per_step_disp_dist << endl;
		return;

	}


	float temp_x = current_pos.x; float temp_y = current_pos.y; float temp_z = current_pos.z;
	int temp_row = current_cell->r; int temp_col = current_cell->c; int temp_layer = current_cell->l;
	float diff_x, diff_y, diff_z;

	//calculate scalar for unit vector
	//calculate difference in destination to current position
	diff_x = x - current_pos.x; diff_y = y - current_pos.y; diff_z = -1 * (z - current_pos.z); //multiply by -1 so that i move in the correct direction! vector maths is assuming the z axis is the reverse of mine
	//calculate the magnitude of that vector
	float mag = sqrt(pow(diff_x, 2) + pow(diff_y, 2) + pow(diff_z, 2));
	float unit_x, unit_y, unit_z;
	//find its unit vector
	unit_x = diff_x / mag; unit_y = diff_y / mag; unit_z = diff_z / mag;
	//cout << "unit_z: " << unit_z << endl;
	//scale it by the magnitude of one cell's resolution
	float scaled_x, scaled_y, scaled_z;
	//scaled_x = unit_x * cell_mag; scaled_y = unit_y * cell_mag; scaled_z = unit_z * cell_mag;
	scaled_x = unit_x * lands->get_land_att(8); scaled_y = unit_y * lands->get_land_att(8); scaled_z = unit_z * lands->get_land_att(8);


	bool solid = false;
	bool too_far = false;
	int counter = 0;
	while (temp_row != row || temp_col != col || temp_layer != layer) {
		if (too_far == true) {  break; }
		//check whether the next step would take it past its max distance
		
		if ((dinfo.per_step_disp_dist + calc_distance(temp_x + scaled_x, temp_y + scaled_y, temp_z - scaled_z)) > max_dist) {
			
			//cout << "indiv " << ID << " would go too far at counter " << counter << endl;
			float scaled_mag = sqrt(pow(scaled_x, 2) + pow(scaled_y, 2) + pow(scaled_z, 2));
			if (max_dist <= scaled_mag) {
				temp_x = x; temp_y = y; temp_z = z; //to test the destination coordinates
				//cout << "indiv " << ID << " max_dist was " << max_dist << endl;
				too_far = true; //break; //I dont break yet because I want it to check the cell for validity and take the step etc
			}
			else {
				float to_go = max_dist - dinfo.per_step_disp_dist;
				temp_x += (unit_x * to_go); temp_y += (unit_y*to_go); temp_z += (unit_z*to_go);
				too_far = true;
				//cout << "indiv " << ID << " too far at counter " << counter << endl;
				//cout << "indiv has travelled " << dinfo.per_step_disp_dist << " in this timestep" << endl;
				//break;
			}

			//diff_x = (temp_x + scaled_x) - current_pos.x; diff_y = (temp_y + scaled_y) - current_pos.y; diff_z = (temp_z - scaled_z) + current_pos.z;
			////cout << "diff_x=" << diff_x << ", diff_y=" << diff_y << ", diff_z=" << diff_z << endl;
			//temp_x += diff_x; temp_y += diff_y; temp_z -= diff_z;
			//cout << "temp_x=" << temp_x << ", temp_y= " << temp_y << ", temp_z= " << temp_z << endl;
			
		}
		else { //if not, continue with the calculated unit vector
			//set up to check the next cell along
			temp_x += scaled_x; temp_y += scaled_y; temp_z -= scaled_z;
		}
		indexes = calc_r_c_l(temp_x, temp_y, temp_z);
		temp_row = indexes[0]; temp_col = indexes[1]; temp_layer = indexes[2];

		//buoyancy limits and the surface are treated as REFLECTIVE BOUNDARIES
		if (dinfo.buoy_min != -9 && temp_z < dinfo.buoy_min) { //if the indiv's minimum buoyancy depth is still below the surface, it shouldn't be allowed to reach the surface
			//cout << "hit upper buoy limit" << endl;
			temp_z = dinfo.buoy_min + 1;
			temp_layer = floor((temp_z - lands->get_land_att(13)) / lands->get_land_att(8)); //current z - minimum z / interval between layers
			//cout << "new temp_z is" << temp_z << ", and buoy layer is " << buoy_layer << endl;

			scaled_z = 0;
		}
		else if (dinfo.buoy_max != -9 && temp_z > dinfo.buoy_max) { //if indiv has gone deeper than buoyancy maximum
			//cout << "hit lower buoy limit" << endl;
			temp_z = dinfo.buoy_max - 1; //i want it to stay just above the buoyancy maximum
			temp_layer = floor((temp_z - lands->get_land_att(13)) / lands->get_land_att(8)); //current z - minimum z / interval between layers
			//cout << "new temp_z is" << temp_z << ", and buoy layer is " << buoy_layer << endl;

			scaled_z = 0;

		}
		else if (dinfo.buoy_min == -9 && temp_z < 0) { //if indiv has hit the surface and there is no buoyancy limit

			//cout << "hit surface" << endl;
			temp_z = 1;
			temp_layer = 0;
			scaled_z = 0;//make this 0 so that it doesnt increment anymore
		}
		//the case of an indiv going too low when buoy_max==-9 is covered in the check_cell function
		if (dinfo.diel_vert == 1 && dvm == true) { //check that it is old/big enough to undergo diel vertical migration

			int min_z, max_z, buoy_layer;
			//check the disp time. remember that it was released at night and there is a 12 hour day night cycle
			//so disp-time/12, if it's odd then it's day time? even then nighttime?
			//cout << "indiv has had " << dinfo.disp_time << " hours of dispersal" << endl;
			int portion = floor(dinfo.disp_time / 12);
			if (portion % 2 == 0) { //even=night

				min_z = dinfo.buoy_min;
				max_z = dinfo.buoy_min + dinfo.dv_range;
				if (max_z > current_cell->seafloor_depth) { //if this takes it too deep
					max_z = current_cell->seafloor_depth - 1;
				}
				//cout << "it's night so indiv between " << min_z <<" and " << max_z << endl;
			}
			else { //odd=day
				//cout << "buoy max is " << dinfo.buoy_max << " seafloor depth is " << current_cell->seafloor_depth;
				if (dinfo.buoy_max > current_cell->seafloor_depth) { //if the seafloor is shallower than the buoyancy max
					min_z = current_cell->seafloor_depth - dinfo.dv_range;
					max_z = current_cell->seafloor_depth - 1;
					if (min_z < 0) { min_z = 0; }
				}
				else {
					max_z = dinfo.buoy_max;
					min_z = dinfo.buoy_max - dinfo.dv_range;
				}
				//cout << "it's day so indiv between " << min_z <<" and " << max_z << endl;
			}
			if (temp_z < min_z) {
				temp_z = min_z + 1;
				temp_layer = floor((temp_z - lands->get_land_att(13)) / lands->get_land_att(8));;
			}
			if (temp_z > max_z) {
				temp_z = max_z - 1;
				temp_layer= floor((temp_z - lands->get_land_att(13)) / lands->get_land_att(8));
			}
		}
		else {
			
			if (temp_z > current_cell->seafloor_depth) {
				temp_z = current_cell->seafloor_depth - 1;
			}
		}

		//if (temp_z < 0) { cout << "when taking steps, but before conditions, z <0: " << temp_z << ", cell mag was " << mag << endl; }
		valid = check_cell(temp_row, temp_col, temp_layer, temp_x, temp_y, temp_z);
		//if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " after check cell 2, per step was 0" << endl; }

		if (valid == false && dead == true) {
			//cout << "indiv " << ID << " while loop " << counter << endl;
			return; 
		}
		//if (temp_z > lands->raster_3d[temp_row][temp_col][temp_layer]->seafloor_depth) { //THIS MIGHT NEED ADJUSTING IN FUTURE TO BE DEPENDENT ON COMPETENCY!
		//	temp_z = lands->raster_3d[temp_row][temp_col][temp_layer]->seafloor_depth - 1;
		//	//recalculate indexes because it might have jumped layers
		//	indexes = calc_r_c_l(temp_x, temp_y, temp_z);
		//	row = indexes[0]; col = indexes[1]; layer = indexes[2];
		//	bool valid = check_cell(temp_row, temp_col, temp_layer, temp_x, temp_y, temp_z);
		//	//if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " after check cell 1, per step was 0" << endl; }
		//	if (valid == false) {
		//		return;
		//	} //will only return false if it was out of xy bounds, lower z bounds or null
		//}
		if (valid == false) { cout << "cell invalid at " << temp_row << ", " << temp_col << ", " << temp_layer << endl; }
		if(lands->raster_3d[temp_row][temp_col][temp_layer]==nullptr){ cout << "cell invalid at " << temp_row << ", " << temp_col << ", " << temp_layer << endl; }
		if (lands->raster_3d[temp_row][temp_col][temp_layer]->habtype != 0) { //not water

			if (dinfo.phase == 0) {
				//if it hits the cell above where absolute seafloor depth is,
				if (temp_layer == 0) { solid = true; break; }
				if (temp_z < lands->raster_3d[temp_row][temp_col][temp_layer]->seafloor_depth && lands->raster_3d[temp_row][temp_col][temp_layer - 1]->habtype == 0) {
					temp_z = lands->raster_3d[temp_row][temp_col][temp_layer - 1]->max_depth - 1;
					temp_layer -= 1;
				}
				else {
					solid = true; break;
				}

			}
			else { //active disperser
				//if it's just the seafloor/unsuitable habitat, and above it is water, continue above it
				bool above = true;

				//under what circumstances do indivs want to check solid habitat?
				if (dinfo.comp == true) { //competent
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->K > 0) { //suitable
						if (mod->patchmodel == 1) {
							if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch != nullptr) { //part of a patch
								if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch->patch_ID != natal_patch) { //not natal
									above = false;
								}
								//if it's the natal patch, then it avoids
							}
							//if it's not part of a patch , then it avoids
						}
						else {
							if (lands->raster_3d[temp_row][temp_col][temp_layer]->cell_ID != natal_cell) { //not natal
								above = false;
							}
							//if it's the natal cell, then it avoids
						}
					}//if it's unsuitable, it avoids

				} //if not competent, then it avoids


				if (above == true) { //if it has passed the checks so far
					//check it's even possible
					if (temp_layer != 0 && lands->raster_3d[temp_row][temp_col][temp_layer - 1] != nullptr) {
						if (lands->raster_3d[temp_row][temp_col][temp_layer - 1]->habtype == 0) {
							//cout << "above is true, moving up ";
							temp_layer -= 1; //move up a layer
							temp_z = lands->raster_3d[temp_row][temp_col][temp_layer]->max_depth - 2; //2 meters above the maxdepth
							//cout << "temp_z is now " << temp_z << ", temp_layer is now " << temp_layer << endl;
							scaled_z = 0; //make this 0 so that it doesnt increment anymore
							z = temp_z; layer = temp_layer; //change the goal to not change depth anymore
						}
						else {
							above = false;
						}
					}
					else {
						above = false;
					}
				}
				if (above == false) { //wants to check or cant go above
					//cout << "above is false, check for settle ";
					solid = true; break;
				}

			}
		}
		//whatever the individual has done during transfer, if it did not break for whatever reason,
		//update how far it has travelled this timestep
		dinfo.per_step_disp_dist += calc_distance(temp_x, temp_y, temp_z);
		if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " met water, per step is 0" << endl; }

		//update current position in "realtime"
		current_pos.x = temp_x; current_pos.y = temp_y; current_pos.z = temp_z;
		//current_cell = lands->raster_3d[temp_row][temp_col][temp_layer];
		counter++;

		//cout << endl;
	}
	
	if (solid == true) {//if something solid stopped it
		//if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " before met solid, per step is 0" << endl; }

		if (dinfo.phase == 0) {
			status = 9; //forced settlement
			//suitable habitat, then it can just settle 
			if (lands->raster_3d[temp_row][temp_col][temp_layer]->K > 0) {
				check_settle(lands->raster_3d[temp_row][temp_col][temp_layer], temp_x, temp_y, temp_z);
				if (sett == true) {
					status = 9; //should still know that it was forced to settle
					//check_Settle has already updated position etc
					
					return;
				}
				else {
					if (dead == true) { return; }
					//cout << "evade 1" << endl;
					float evade_phi = atan2(temp_y - current_pos.y, temp_x - current_pos.x); //direction it was going to go
					vector<float> new_coords = evade(evade_phi);
					indexes = calc_r_c_l(new_coords[0], new_coords[1], new_coords[2]);
					dinfo.per_step_disp_dist += calc_distance(new_coords[0], new_coords[1], new_coords[2]);
					if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << "evaded, per step is 0" << endl; }
					//current_pos.x = new_coords[0]; current_pos.y = new_coords[1]; current_pos.z = new_coords[2];
					update_pos(new_coords[0], new_coords[1], new_coords[2], indexes[0], indexes[1], indexes[2], true, false, true, true); //reset memory so it doesnt affect dp
					sett = false;
					status = 1; //need to make sure it goes back to just being a disperser
					//cout << "indiv " << ID << " hit suitable habitat but didnt settle, evaded," << endl;
					return;
				}
			}
			//unsuitable (which means if they settle, they die
			else {
				float ran_sett = unif_dist(eng);
				//even with unsuitable habitat, i want stochasticity, so they only are forced to settle and die 80% of the time
				if (ran_sett < 0.5) {
					dinfo.per_step_disp_dist += calc_distance(temp_x, temp_y, temp_z);
					update_pos(temp_x, temp_y, temp_z, temp_row, temp_col, temp_layer, false, false, false, true);
					dead = true;
					sett = false,
					status = 4;
					return;
				}
				else { //otherwise, they evade 
					//cout << "evade" << endl;
					float evade_phi = atan2(temp_y - current_pos.y, temp_x - current_pos.x); //direction it was going to go
					vector<float> new_coords = evade(evade_phi);
					indexes = calc_r_c_l(new_coords[0], new_coords[1], new_coords[2]);
					dinfo.per_step_disp_dist += calc_distance(new_coords[0], new_coords[1], new_coords[2]);
					//if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " unsuit evade, per step is 0" << endl; }
					update_pos(new_coords[0], new_coords[1], new_coords[2], indexes[0], indexes[1], indexes[2], true, false, true, true); //reset memory so it doesnt affect dp
					status = 1; //need to make sure it goes back to just being a disperser
					sett = false;

					return;
				}
			}
		}
		else {
			if (lands->raster_3d[temp_row][temp_col][temp_layer]->K > 0 && lands->raster_3d[temp_row][temp_col][temp_layer]->patch != nullptr) { //if it's suitable habitat
			//cout << "suitable habitat" << endl;

				bool want_check = true; //does indiv even want to check settlement here?
				//cout << "checking want_check for patch " << lands->raster_3d[temp_row][temp_col][temp_layer]->patch->patch_ID << endl;
				if (dinfo.comp == false) { want_check = false; /*cout << "not comp yet" << endl;*/ }
				else if (mod->patchmodel == 0) {
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->cell_ID == natal_cell) { want_check = false; } //if its the natal cell
					else if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), lands->raster_3d[temp_row][temp_col][temp_layer]->cell_ID)) { //or if it's tried before
						want_check = false;
					}
				}
				else if (mod->patchmodel == 1) { //patch based
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch != nullptr && lands->raster_3d[temp_row][temp_col][temp_layer]->patch->patch_ID == natal_patch) {
						want_check = false; //cout << "natal" << endl;
					}
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch != nullptr && count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), lands->raster_3d[temp_row][temp_col][temp_layer]->patch->patch_ID)) {
						want_check = false; //cout << "tested before" << endl;
					}
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch == nullptr) { //suitable but not belonging to a patch
						want_check = false; //cout << "doesnt belong to a patch" << endl;
					}

				}

				if (want_check == true) {
					//cout << "checking settlement" << endl;
					check_settle(lands->raster_3d[temp_row][temp_col][temp_layer], temp_x, temp_y, temp_z);
				}

				if (sett == true) { //if it decided to settle
					status = 3;
					if (current_cell->patch == nullptr) { cout << "in active traverse, sett is true but patch was null" << endl; }
					return;
				}
				else {
					//if it didn't settle, i want it to evade, change direction
					float evade_phi = atan2(temp_y - current_pos.y, temp_x - current_pos.x); //direction it was going to go
					vector<float> new_coords = evade(evade_phi);
					indexes = calc_r_c_l(new_coords[0], new_coords[1], new_coords[2]);
					dinfo.per_step_disp_dist += calc_distance(new_coords[0], new_coords[1], new_coords[2]);
					if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << "evaded from suit, per step is 0" << endl; }
					update_pos(new_coords[0], new_coords[1], new_coords[2], indexes[0], indexes[1], indexes[2], true, false, true, true); //reset memory so it doesnt affect dp
					status = 1;
					if (current_cell->cell_speed == -9999) {
						cout << "indiv " << ID << " hit unsuitable habitat, evaded, current cell speed is  " << current_cell->cell_speed << endl;
					}

					sett = false;
				}

			}
			else { //unsuitable habitat

				float evade_phi = atan2(temp_y - current_pos.y, temp_x - current_pos.x); //direction it was going to go
				vector<float> new_coords = evade(evade_phi);
				indexes = calc_r_c_l(new_coords[0], new_coords[1], new_coords[2]);
				dinfo.per_step_disp_dist += calc_distance(new_coords[0], new_coords[1], new_coords[2]);
				if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " evaded from unsuit, per step is 0" << endl; }
				//cout << "in while loop: indexes after evade were " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << endl;
				if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]] == nullptr) { cout << "in while loop: indexes after evade were " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << "and cell was null! " << endl; }

				update_pos(new_coords[0], new_coords[1], new_coords[2], indexes[0], indexes[1], indexes[2], true, false, true, true); //reset memory so it doesnt affect dp
			}
		}


	}
	else {
		//reached its destination and it was water
		if (temp_row == row && temp_col == col) {
			temp_x = x; temp_y = y; //always use temp_z because that is most likely to have changed
		}
		dinfo.per_step_disp_dist += calc_distance(temp_x, temp_y, temp_z);
		if (dinfo.per_step_disp_dist == 0) { 
			cout << "indiv " << ID << " met distination, per step is 0" << endl; 
			if (solid == true) { cout << "solid was true" << endl; }
			if (too_far == true) { cout << "too far was true at counter " << counter << " with max dist " << max_dist << endl; }
		}
		
		update_pos(temp_x, temp_y, temp_z, temp_row, temp_col, temp_layer, true, false, true, true);
		return;
	}
	
	return;
}

void Individual::traverse3(float x, float y, float z) {
	int row, col, layer;
	vector<int> indexes;
	float max_dist;
	bool valid;
	if (dinfo.phase == 0) {
		max_dist = current_cell->cell_speed;
	}
	else {
		max_dist = calc_distance(x, y, z);
	}
	if (dead == true) {
		return;
	}
	//calculate row, col and layer of destination
	indexes = calc_r_c_l(x, y, z);
	row = indexes[0]; col = indexes[1]; layer = indexes[2];
	//cout << "indiv " << ID;
	//valid = check_cell(row, col, layer, x, y, z);
	//cout << " adjusting" << endl;
	//if from the very beginning the movement is within the same cell, 
	if (current_cell->r == row && current_cell->c == col && current_cell->l == layer) {
		//cout << "same cell" << endl;
		dinfo.per_step_disp_dist += calc_distance(x, y, z);
		update_pos(x, y, z, row, col, layer, true, false, false, true); //update current_pos, memory, not current_cell and dont empty memory
		//dont need to update cell because it hasnt left the cell
		if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " same cell, per step was 0, max dist was " << max_dist << endl; }
		//cout << "perstep disp dist is now " << dinfo.per_step_disp_dist << endl;
		return;

	}


	float temp_x = current_pos.x; float temp_y = current_pos.y; float temp_z = current_pos.z;
	int temp_row = current_cell->r; int temp_col = current_cell->c; int temp_layer = current_cell->l;
	float diff_x, diff_y, diff_z;

	//calculate scalar for unit vector
	//calculate difference in destination to current position
	diff_x = x - current_pos.x; diff_y = y - current_pos.y; diff_z = -1 * (z - current_pos.z); //multiply by -1 so that i move in the correct direction! vector maths is assuming the z axis is the reverse of mine
	//calculate the magnitude of that vector
	float mag = sqrt(pow(diff_x, 2) + pow(diff_y, 2) + pow(diff_z, 2));
	float unit_x, unit_y, unit_z;
	//find its unit vector
	unit_x = diff_x / mag; unit_y = diff_y / mag; unit_z = diff_z / mag;
	//cout << "unit_z: " << unit_z << endl;
	//scale it by the magnitude of one cell's resolution
	float scaled_x, scaled_y, scaled_z;
	//scaled_x = unit_x * cell_mag; scaled_y = unit_y * cell_mag; scaled_z = unit_z * cell_mag;
	scaled_x = unit_x * lands->get_land_att(8); scaled_y = unit_y * lands->get_land_att(8); scaled_z = unit_z * lands->get_land_att(8);


	bool solid = false;
	bool too_far = false;
	int counter = 0;
	while (temp_row != row || temp_col != col || temp_layer != layer) {
		if (too_far == true) { break; }
		//check whether the next step would take it past its max distance

		if ((dinfo.per_step_disp_dist + calc_distance(temp_x + scaled_x, temp_y + scaled_y, temp_z - scaled_z)) > max_dist) {

			//cout << "indiv " << ID << " would go too far at counter " << counter << endl;
			float scaled_mag = sqrt(pow(scaled_x, 2) + pow(scaled_y, 2) + pow(scaled_z, 2));
			if (max_dist <= scaled_mag) {
				temp_x = x; temp_y = y; temp_z = z; //to test the destination coordinates
				//cout << "indiv " << ID << " max_dist was " << max_dist << endl;
				too_far = true; //break; //I dont break yet because I want it to check the cell for validity and take the step etc
			}
			else {
				float to_go = max_dist - dinfo.per_step_disp_dist;
				temp_x += (unit_x * to_go); temp_y += (unit_y*to_go); temp_z += (unit_z*to_go);
				too_far = true;
				//cout << "indiv " << ID << " too far at counter " << counter << endl;
				//cout << "indiv has travelled " << dinfo.per_step_disp_dist << " in this timestep" << endl;
				//break;
			}

			//diff_x = (temp_x + scaled_x) - current_pos.x; diff_y = (temp_y + scaled_y) - current_pos.y; diff_z = (temp_z - scaled_z) + current_pos.z;
			////cout << "diff_x=" << diff_x << ", diff_y=" << diff_y << ", diff_z=" << diff_z << endl;
			//temp_x += diff_x; temp_y += diff_y; temp_z -= diff_z;
			//cout << "temp_x=" << temp_x << ", temp_y= " << temp_y << ", temp_z= " << temp_z << endl;

		}
		else { //if not, continue with the calculated unit vector
			//set up to check the next cell along
			temp_x += scaled_x; temp_y += scaled_y; temp_z -= scaled_z;
		}
		indexes = calc_r_c_l(temp_x, temp_y, temp_z);
		temp_row = indexes[0]; temp_col = indexes[1]; temp_layer = indexes[2];

		//buoyancy limits and the surface are treated as REFLECTIVE BOUNDARIES
		if (dinfo.buoy_min != -9 && temp_z < dinfo.buoy_min) { //if the indiv's minimum buoyancy depth is still below the surface, it shouldn't be allowed to reach the surface
			//cout << "hit upper buoy limit" << endl;
			temp_z = dinfo.buoy_min + 1;
			temp_layer = floor((temp_z - lands->get_land_att(13)) / lands->get_land_att(8)); //current z - minimum z / interval between layers
			//cout << "new temp_z is" << temp_z << ", and buoy layer is " << buoy_layer << endl;

			scaled_z = 0;
		}
		else if (dinfo.buoy_max != -9 && temp_z > dinfo.buoy_max) { //if indiv has gone deeper than buoyancy maximum
			//cout << "hit lower buoy limit" << endl;
			temp_z = dinfo.buoy_max - 1; //i want it to stay just above the buoyancy maximum
			temp_layer = floor((temp_z - lands->get_land_att(13)) / lands->get_land_att(8)); //current z - minimum z / interval between layers
			//cout << "new temp_z is" << temp_z << ", and buoy layer is " << buoy_layer << endl;

			scaled_z = 0;

		}
		else if (dinfo.buoy_min == -9 && temp_z < 0) { //if indiv has hit the surface and there is no buoyancy limit

			//cout << "hit surface" << endl;
			temp_z = 1;
			temp_layer = 0;
			scaled_z = 0;//make this 0 so that it doesnt increment anymore
		}
		//the case of an indiv going too low when buoy_max==-9 is covered in the check_cell function
		//check seafloor depth
		float sf = 0; int sf_layer = 0;
		if (lands->raster_3d[temp_row][temp_col][0] != nullptr) {
			
			sf_layer = floor((lands->raster_3d[temp_row][temp_col][0]->seafloor_depth - lands->get_land_att(13)) / lands->get_land_att(8));
			if (lands->raster_3d[temp_row][temp_col][0]->seafloor_depth >= (sf_layer * lands->get_land_att(8) + (lands->get_land_att(8) / 2))) { //if seafloor is deeper than halfway through this cell
				//then the cell is water
				sf = lands->raster_3d[temp_row][temp_col][0]->seafloor_depth;
			}
			else { //the cell is habitat
				sf = lands->raster_3d[temp_row][temp_col][sf_layer]->min_depth - 1;
			}
		}
		else { //otherwise use info from current cell
			sf_layer = floor((current_cell->seafloor_depth - lands->get_land_att(13)) / lands->get_land_att(8));
			if (lands->raster_3d[current_cell->r][current_cell->c][sf_layer]->habtype != 0) { //it's solid, that means seafloor depth was less than the halfway depth of the cell, so whole cell is solid
				sf = lands->raster_3d[current_cell->r][current_cell->c][sf_layer]->min_depth - 1;
			}
			else { //water, so can use seafloor depth as a measure
				sf = current_cell->seafloor_depth - 1;
			}
		}
		
		if (dinfo.diel_vert == 1 && dvm == true) { //check that it is old/big enough to undergo diel vertical migration

			int min_z, max_z, buoy_layer;
			//check the disp time. remember that it was released at night and there is a 12 hour day night cycle
			//so disp-time/12, if it's odd then it's day time? even then nighttime?
			//cout << "indiv has had " << dinfo.disp_time << " hours of dispersal" << endl;
			int portion = floor(dinfo.disp_time / 12);
			if (portion % 2 == 0) { //even=night so at the surface

				min_z = dinfo.buoy_min;
				max_z = dinfo.buoy_min + dinfo.dv_range;

				//check seafloor
				if (max_z > sf) { //if this takes it too deep
					max_z = sf;
				}
				//cout << "it's night so indiv between " << min_z <<" and " << max_z << endl;
			}
			else { //odd=day
				max_z = dinfo.buoy_max;
				min_z = dinfo.buoy_max - dinfo.dv_range;
				//cout << "buoy max is " << dinfo.buoy_max << " seafloor depth is " << current_cell->seafloor_depth;
				if (dinfo.buoy_max >= sf) {
					min_z = sf - dinfo.dv_range;
					max_z = sf;
				}
				//cout << "it's day so indiv between " << min_z <<" and " << max_z << endl;
			}
			if (temp_z < min_z) {
				temp_z = min_z + 1;
				temp_layer = floor((temp_z - lands->get_land_att(13)) / lands->get_land_att(8));
			}
			if (temp_z > max_z) {
				temp_z = max_z - 1;
				temp_layer = floor((temp_z - lands->get_land_att(13)) / lands->get_land_att(8));
			}
		}
		else {
			if (temp_z > sf) {
				temp_z = sf - 1;
				temp_layer = floor((sf - lands->get_land_att(13)) / lands->get_land_att(8));
			}
		}

		//if (temp_z < 0) { cout << "when taking steps, but before conditions, z <0: " << temp_z << ", cell mag was " << mag << endl; }
		valid = check_cell(temp_row, temp_col, temp_layer, temp_x, temp_y, temp_z);
		//if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " after check cell 2, per step was 0" << endl; }

		if (valid == false && dead == true) {
			//cout << "indiv " << ID << " while loop " << counter << endl;
			return;
		}
		//if (temp_z > lands->raster_3d[temp_row][temp_col][temp_layer]->seafloor_depth) { //THIS MIGHT NEED ADJUSTING IN FUTURE TO BE DEPENDENT ON COMPETENCY!
		//	temp_z = lands->raster_3d[temp_row][temp_col][temp_layer]->seafloor_depth - 1;
		//	//recalculate indexes because it might have jumped layers
		//	indexes = calc_r_c_l(temp_x, temp_y, temp_z);
		//	row = indexes[0]; col = indexes[1]; layer = indexes[2];
		//	bool valid = check_cell(temp_row, temp_col, temp_layer, temp_x, temp_y, temp_z);
		//	//if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " after check cell 1, per step was 0" << endl; }
		//	if (valid == false) {
		//		return;
		//	} //will only return false if it was out of xy bounds, lower z bounds or null
		//}

		if (lands->raster_3d[temp_row][temp_col][temp_layer]->habtype != 0) { //not water

			if (dinfo.phase == 0) {
				//if it hits the cell above where absolute seafloor depth is, 
				if (temp_z < lands->raster_3d[temp_row][temp_col][temp_layer]->seafloor_depth && lands->raster_3d[temp_row][temp_col][temp_layer-1]->habtype==0) {
					temp_z = lands->raster_3d[temp_row][temp_col][temp_layer - 1]->max_depth - 1;
					temp_layer -= 1;
				}
				else {
					solid = true; break;
				}
				
			}
			else { //active disperser
				//if it's just the seafloor/unsuitable habitat, and above it is water, continue above it
				bool above = true;

				//under what circumstances do indivs want to check solid habitat?
				if (dinfo.comp == true) { //competent
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->K > 0) { //suitable
						if (mod->patchmodel == 1) {
							if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch != nullptr) { //part of a patch
								if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch->patch_ID != natal_patch) { //not natal
									above = false;
								}
								//if it's the natal patch, then it avoids
							}
							//if it's not part of a patch , then it avoids
						}
						else {
							if (lands->raster_3d[temp_row][temp_col][temp_layer]->cell_ID != natal_cell) { //not natal
								above = false;
							}
							//if it's the natal cell, then it avoids
						}
					}//if it's unsuitable, it avoids

				} //if not competent, then it avoids


				if (above == true) { //if it has passed the checks so far
					//check it's even possible
					if (temp_layer != 0 && lands->raster_3d[temp_row][temp_col][temp_layer - 1] != nullptr) {
						if (lands->raster_3d[temp_row][temp_col][temp_layer - 1]->habtype == 0) {
							//cout << "above is true, moving up ";
							temp_layer -= 1; //move up a layer
							temp_z = lands->raster_3d[temp_row][temp_col][temp_layer]->max_depth - 2; //2 meters above the maxdepth
							//cout << "temp_z is now " << temp_z << ", temp_layer is now " << temp_layer << endl;
							scaled_z = 0; //make this 0 so that it doesnt increment anymore
							z = temp_z; layer = temp_layer; //change the goal to not change depth anymore
						}
						else {
							above = false;
						}
					}
					else {
						above = false;
					}
				}
				if (above == false) { //wants to check or cant go above
					//cout << "above is false, check for settle ";
					solid = true; break;
				}

			}
		}
		//whatever the individual has done during transfer, if it did not break for whatever reason,
		//update how far it has travelled this timestep
		dinfo.per_step_disp_dist += calc_distance(temp_x, temp_y, temp_z);
		if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " met water, per step is 0" << endl; }

		//update current position in "realtime"
		current_pos.x = temp_x; current_pos.y = temp_y; current_pos.z = temp_z;
		//current_cell = lands->raster_3d[temp_row][temp_col][temp_layer];
		counter++;

		//cout << endl;
	}

	if (solid == true) {//if something solid stopped it
		//if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " before met solid, per step is 0" << endl; }

		if (dinfo.phase == 0) {

			status = 9; //forced settlement
			//suitable habitat, then it can just settle 
			if (lands->raster_3d[temp_row][temp_col][temp_layer]->K > 0) {
				check_settle(lands->raster_3d[temp_row][temp_col][temp_layer], temp_x, temp_y, temp_z);
				if (sett == true) {
					status = 9; //should still know that it was forced to settle
					//check_Settle has already updated position etc

					return;
				}
				else {
					if (dead == true) { return; }
					//cout << "evade 1" << endl;
					float evade_phi = atan2(temp_y - current_pos.y, temp_x - current_pos.x); //direction it was going to go
					vector<float> new_coords = evade(evade_phi);
					indexes = calc_r_c_l(new_coords[0], new_coords[1], new_coords[2]);
					dinfo.per_step_disp_dist += calc_distance(new_coords[0], new_coords[1], new_coords[2]);
					if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << "evaded, per step is 0" << endl; }
					//current_pos.x = new_coords[0]; current_pos.y = new_coords[1]; current_pos.z = new_coords[2];
					update_pos(new_coords[0], new_coords[1], new_coords[2], indexes[0], indexes[1], indexes[2], true, false, true, true); //reset memory so it doesnt affect dp
					sett = false;
					status = 1; //need to make sure it goes back to just being a disperser
					//cout << "indiv " << ID << " hit suitable habitat but didnt settle, evaded," << endl;
					return;
				}
			}
			//unsuitable (which means if they settle, they die
			else {
				float ran_sett = unif_dist(eng);
				//even with unsuitable habitat, i want stochasticity, so they only are forced to settle and die 80% of the time
				if (ran_sett < 0.5) {
					dinfo.per_step_disp_dist += calc_distance(temp_x, temp_y, temp_z);
					update_pos(temp_x, temp_y, temp_z, temp_row, temp_col, temp_layer, false, false, false, true);
					dead = true;
					sett = false,
					status = 4;
					return;
				}
				else { //otherwise, they evade 
					//cout << "evade" << endl;
					float evade_phi = atan2(temp_y - current_pos.y, temp_x - current_pos.x); //direction it was going to go
					vector<float> new_coords = evade(evade_phi);
					indexes = calc_r_c_l(new_coords[0], new_coords[1], new_coords[2]);
					dinfo.per_step_disp_dist += calc_distance(new_coords[0], new_coords[1], new_coords[2]);
					if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " unsuit evade, per step is 0" << endl; }
					update_pos(new_coords[0], new_coords[1], new_coords[2], indexes[0], indexes[1], indexes[2], true, false, true, true); //reset memory so it doesnt affect dp
					status = 1; //need to make sure it goes back to just being a disperser
					sett = false;

					return;
				}
			}
		}
		else {
			if (lands->raster_3d[temp_row][temp_col][temp_layer]->K > 0 && lands->raster_3d[temp_row][temp_col][temp_layer]->patch != nullptr) { //if it's suitable habitat
			//cout << "suitable habitat" << endl;

				bool want_check = true; //does indiv even want to check settlement here?
				//cout << "checking want_check for patch " << lands->raster_3d[temp_row][temp_col][temp_layer]->patch->patch_ID << endl;
				if (dinfo.comp == false) { want_check = false; /*cout << "not comp yet" << endl;*/ }
				else if (mod->patchmodel == 0) {
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->cell_ID == natal_cell) { want_check = false; } //if its the natal cell
					else if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), lands->raster_3d[temp_row][temp_col][temp_layer]->cell_ID)) { //or if it's tried before
						want_check = false;
					}
				}
				else if (mod->patchmodel == 1) { //patch based
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch != nullptr && lands->raster_3d[temp_row][temp_col][temp_layer]->patch->patch_ID == natal_patch) {
						want_check = false; //cout << "natal" << endl;
					}
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch != nullptr && count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), lands->raster_3d[temp_row][temp_col][temp_layer]->patch->patch_ID)) {
						want_check = false; //cout << "tested before" << endl;
					}
					if (lands->raster_3d[temp_row][temp_col][temp_layer]->patch == nullptr) { //suitable but not belonging to a patch
						want_check = false; //cout << "doesnt belong to a patch" << endl;
					}

				}

				if (want_check == true) {
					//cout << "checking settlement" << endl;
					check_settle(lands->raster_3d[temp_row][temp_col][temp_layer], temp_x, temp_y, temp_z);
				}

				if (sett == true) { //if it decided to settle
					status = 3;
					if (current_cell->patch == nullptr) { cout << "in active traverse, sett is true but patch was null" << endl; }
					return;
				}
				else {
					//if it didn't settle, i want it to evade, change direction
					float evade_phi = atan2(temp_y - current_pos.y, temp_x - current_pos.x); //direction it was going to go
					vector<float> new_coords = evade(evade_phi);
					indexes = calc_r_c_l(new_coords[0], new_coords[1], new_coords[2]);
					dinfo.per_step_disp_dist += calc_distance(new_coords[0], new_coords[1], new_coords[2]);
					if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << "evaded from suit, per step is 0" << endl; }
					update_pos(new_coords[0], new_coords[1], new_coords[2], indexes[0], indexes[1], indexes[2], true, false, true, true); //reset memory so it doesnt affect dp
					status = 1;
					if (current_cell->cell_speed == -9999) {
						cout << "indiv " << ID << " hit unsuitable habitat, evaded, current cell speed is  " << current_cell->cell_speed << endl;
					}

					sett = false;
				}

			}
			else { //unsuitable habitat

				float evade_phi = atan2(temp_y - current_pos.y, temp_x - current_pos.x); //direction it was going to go
				vector<float> new_coords = evade(evade_phi);
				indexes = calc_r_c_l(new_coords[0], new_coords[1], new_coords[2]);
				dinfo.per_step_disp_dist += calc_distance(new_coords[0], new_coords[1], new_coords[2]);
				if (dinfo.per_step_disp_dist == 0) { cout << "indiv " << ID << " evaded from unsuit, per step is 0" << endl; }
				//cout << "in while loop: indexes after evade were " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << endl;
				if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]] == nullptr) { cout << "in while loop: indexes after evade were " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << "and cell was null! " << endl; }

				update_pos(new_coords[0], new_coords[1], new_coords[2], indexes[0], indexes[1], indexes[2], true, false, true, true); //reset memory so it doesnt affect dp
			}
		}


	}
	else {
		//reached its destination and it was water
		if (temp_row == row && temp_col == col) {
			temp_x = x; temp_y = y; //always use temp_z because that is most likely to have changed
		}
		dinfo.per_step_disp_dist += calc_distance(temp_x, temp_y, temp_z);
		if (dinfo.per_step_disp_dist == 0) {
			cout << "indiv " << ID << " met distination, per step is 0" << endl;
			if (solid == true) { cout << "solid was true" << endl; }
			if (too_far == true) { cout << "too far was true at counter " << counter << " with max dist " << max_dist << endl; }
		}

		update_pos(temp_x, temp_y, temp_z, temp_row, temp_col, temp_layer, true, false, true, true);
		return;
	}

	return;
}

void Individual::leave_natal() {
	//individual will take direction it took to move away from natal patch, which is stored in memory
	if (current_cell->cell_speed == -9999) {
		cout << "at leave natal, current cell has no speed" << endl;
	}
	locn start, end, predicted;
	float theta, phi;
	vector<int> indexes;
	int temp_row, temp_col, temp_layer;
	start = memory.front(); end = memory.back();
	predicted.z = end.z; //it doesnt matter whether initial move was up, down, sideways. i dont want them to change depth yet
	theta = acos(0);

	//for (int b = 0; b < pop_sinfo->buffer_xycells; b++) { //so that you move out of the range of the buffer cell for the natal patch
		//end = memory.back();
		if (end.z > lands->get_land_att(14)) {
			cout << " in leave natal, start: " << start.x << " , " << start.y << ", " << start.z;
			cout << ", end: " << end.x << ", " << end.y << ", " << end.z << endl;

		}

		phi = atan2((end.y - start.y), (end.x - start.x));
		if (phi == 0) { //if individual went straight up at juv_release, so the only change was in the z coordinate
			phi = current_cell->cell_angle; //take the angle of the current
		}
		bool valid = false;
		while (valid == false) {
			if (dinfo.phase == 0) { //if it's passive, the step length will depend on current speed
				predicted.x = end.x + (current_cell->cell_speed * sin(theta) * cos(phi));
				predicted.y = end.y + (current_cell->cell_speed * sin(theta) * sin(phi));
			}
			else { //otherwise it depends on swimming ability
				predicted.x = end.x + (dinfo.SL * sin(theta) * cos(phi));
				predicted.y = end.y + (dinfo.SL * sin(theta) * sin(phi));
			}

			indexes = calc_r_c_l(predicted.x, predicted.y, predicted.z);
			temp_row = indexes[0]; temp_col = indexes[1]; temp_layer = indexes[2];
			//cout << "predicted: " << predicted.x << ", " << predicted.y << ", " << predicted.z << endl;

			if (predicted.x < 0 || predicted.x >= lands->get_land_att(11) || predicted.y < 0 || predicted.y >= lands->get_land_att(12)) {
				//cout << "left landscape" << endl;
				dead = true;
				dinfo.per_step_disp_dist += calc_distance(predicted.x, predicted.y, predicted.z);
				status = 8;
				return;
			}
			if (temp_row < 0 || temp_row >= lands->get_land_att(1) || temp_col < 0 || temp_col >= lands->get_land_att(2)) {
				dead = true;
				dinfo.per_step_disp_dist += calc_distance(predicted.x, predicted.y, predicted.z);
				status = 8;
				return;
			}
			//don't need to check z because it isnt changing
			if (lands->raster_3d[temp_row][temp_col][temp_layer] == nullptr ) { //if the chosen cell is a nullptr
				cout << "in leave_natal, cell " << temp_row << ", " << temp_col << ", " << temp_layer << " was null" << endl;
				valid = false; //cout << "turning";
				normal_distribution<float> turn(0, 2);
				float angle = fmod(acos(0) / 2 + turn(eng), M2_PI);
				phi += angle; 
				phi = fmod(phi, M2_PI);
				//cout << phi << " degrees" << endl;
			}
			else if (lands->raster_3d[temp_row][temp_col][temp_layer]->habtype > 0) { //if it's solid
				valid = false; //cout << "turning";
				normal_distribution<float> turn(0, 2);
				float angle = fmod(acos(0) / 2 + turn(eng), M2_PI);
				phi += angle;
				phi = fmod(phi, M2_PI);
			}
			else {
				valid = true;
			}
		}
		dinfo.per_step_disp_dist += calc_distance(predicted.x, predicted.y, predicted.z);
		update_pos(predicted.x, predicted.y, predicted.z, temp_row, temp_col, temp_layer, true, false, true, true);
		dinfo.disp_time++;
	//}
		if (current_cell->cell_speed == -9999) {
			cout << "indiv " << ID << " after leave natal, current cell has no speed" << endl;
		}
}

vector<float> Individual::calc_move() {
	//this is no longer responsible for moving away from unsuitable habitat, it just calculates movements
	vector<float> coordinates;
	float phi, theta;
	cauchy_distribution<double> cauchy, curr_cauchy;
	float new_x, new_y, new_z;
	float temp_x, temp_y, temp_z;
	float dist = 0;
	float propr = lands->get_land_att(8) / lands->get_land_att(7);

	Cell* temp_goal;
	//cout << "this indiv's SL is " << dinfo.SL << endl;
	bool avoid = false;

	//passive, only dependent on the current, can't use a buffer
	if (dinfo.phase == 0) { 
		dist = current_cell->cell_speed;
		theta = acos(current_cell->w / dist);
		//cout << "cell speed is " << dist;
		//try {
		//	if (current_cell->w == -9999) {
		//		throw 9999;
		//	}
		//}
		//catch (int x) {
		//	cout << "settled in an invalid cell" << endl;
		//}

		cauchy = cauchy_distribution<double>(current_cell->cell_angle, -log(.9)); //i am using rho=0.9 for passive individuals
		phi = fmod(cauchy(eng), M2_PI); //random direction in wrapped cauchy distribution, divide by 2pi to put it on the unit circl= "wrap" the cauchy distribution
		//cout << "phi is " << phi << ", theta is " << theta << endl;
		temp_x = dist * sin(theta) * cos(phi);
		temp_y = dist * sin(theta) * sin(phi);
		temp_z = dist * cos(theta); 
		//cout << "temp_x: " << temp_x << ", temp_y: " << temp_y << ", temp_z: " << temp_z << endl;
		new_x = current_pos.x + temp_x;
		new_y = current_pos.y + temp_y;
		new_z = current_pos.z - temp_z; //because my z axis is flipped, a negative temp_z denotes increasing depth, so i need to subtract to switch the sign
		if (new_z > lands->get_land_att(14)) {
			cout << "phi is " << phi << ", theta is " << theta << endl;
			cout << "cell speed is " << dist;
			cout << "current cell habtype is " << current_cell->habtype << endl;
			cout << "starting x " << current_pos.x << " starting y" << current_pos.y << " starting z " << current_pos.z << endl;
			cout << "starting row " << current_cell->r << ", starting col " << current_cell->c << ", starting layer " << current_cell->l << endl;
			cout << "goal x " << new_x << " goal y " << new_y << " goal z " << new_z << endl;

		}


		coordinates.push_back(new_x); coordinates.push_back(new_y); coordinates.push_back(new_z);
		if (calc_distance(new_x, new_y, new_z) == 0) {
			cout << "distance to travel was 0" << endl;
		}
		//cout << "using the current, distance is: " << calc_distance(new_x, new_y, new_z)  << endl;

		return coordinates;
	}
	
	

	//active, dependent on current but can respond to environmental cues, can use a buffer
	else {
		//where is the current going?
		curr_cauchy = cauchy_distribution<double>(current_cell->cell_angle, .9); //add some stochasticity
		float curr_phi = curr_cauchy(eng);
		float curr_theta = acos(current_cell->w / current_cell->cell_speed);
		float curr_x = current_cell->cell_speed * sin(curr_theta) * cos(curr_phi);
		float curr_y = current_cell->cell_speed * sin(curr_theta) * sin(curr_phi);
		float curr_z = current_cell->cell_speed * cos(curr_theta);

		//where is the individual going?
		float dp_theta, dp_phi;
		
		if (dinfo.goal_bias == 1) { //if it locked onto a goal before i want it to continue to that goal, not switch goals
			dist = calc_distance(dinfo.goal_focus->midpoint[0], dinfo.goal_focus->midpoint[1], dinfo.goal_focus->midpoint[2]);
			//cout << "was locked onto a goal, staying with that goal" << endl;
			dp_phi = atan2((dinfo.goal_focus->midpoint[1] - current_pos.y), (dinfo.goal_focus->midpoint[0] - current_pos.x));
			dp_theta = acos(-(dinfo.goal_focus->midpoint[2] - current_pos.z) / dist); //acos(z/r)
			if (dist >= dinfo.SL) { //if the pos goal is closer than the SL
				dist = dinfo.SL; //it still shouldnt travel past its step length
			}
		}
		else if (dinfo.comp == true && pstage->s_sett.buffer==1) {//step 1: explore surroundings, they can only seek suitable habitat if they are competent
			if (current_cell->bufferpos == -9999) { //if it hasnt been assessed yet,
				//cout << "assessing cell buffer" << endl;
				current_cell->assess_pos_buffer(pstage->s_sett.buffer_xycells, pstage->s_sett.buffer_zlayers, lands);
			}
			if (current_cell->bufferpos == 1) { //not just open water 
				//cout << "found pos buffer" << endl;
				for (int p = 0; p < current_cell->buffer_focus.size(); p++) {

					if (mod->patchmodel == 0) {//if indiv has checked that cell for settlement before /if that cell is natal
						if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), current_cell->buffer_focus[p]->cell_ID)) {
							avoid = true;
						}
					}
					else { //if that indiv has tested that patch before / if that patch is natal
						if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), current_cell->buffer_focus[p]->patch->patch_ID)) {
							avoid = true;
						}
					}
				}
				if (avoid == false) {
					//cout << "positive goal bias" << endl;
					random_shuffle(current_cell->buffer_focus.begin(), current_cell->buffer_focus.end());
					temp_goal = current_cell->buffer_focus[0];
					dist = calc_distance(temp_goal->midpoint[0], temp_goal->midpoint[1], temp_goal->midpoint[2]);
					//cout << "chose pos bias temp goal" << endl;
					dinfo.goal_bias == 1; dinfo.goal_focus = temp_goal;

					dp_phi = atan2((temp_goal->midpoint[1] - current_pos.y), (temp_goal->midpoint[0] - current_pos.x));
					//if (temp_goal->midpoint[2] - current_pos.z > 0) { //need to go down
					//	theta = 1.6;
					//}
					//else { //need to go up
					//	theta = 1.55;
					//}
					dp_theta = acos(-(temp_goal->midpoint[2] - current_pos.z) / dist); //acos(z/r)
					if (dist >= dinfo.SL) { //if the pos goal is closer than the SL
						dist = dinfo.SL; //it still shouldnt travel past its step length
					}
				}
			}
		}
		if (dinfo.comp==false || pstage->s_sett.buffer == 0 || current_cell->bufferpos == 0 || avoid == true) {//open water, use memory and current
				//cout << "using memory and current" << endl;
				//float dp_theta, dp_phi;
				//float curr_theta, curr_phi;
			dist = dinfo.SL;
			locn prev = memory.front();
			//cout << "memory is " << memory.size() << " steps long" << endl;

			if (memory.size() >= 2) { //use both memory and current
				//current data
				//cout << "using prev y " << prev.y << " and prev x " << prev.x << endl;
				dp_phi = atan2((current_pos.y - prev.y), (current_pos.x - prev.x)); //angle indiv was travelling at
				//cout << "dp_phi is " << dp_phi << endl;
				dp_theta = acos((current_pos.z - prev.z) / calc_distance(prev.x, prev.y, prev.z)); //acos(z/r)

			}
			else { //sample a new direction using current and rho
				curr_cauchy = cauchy_distribution<double>(current_cell->cell_angle, dinfo.rho);
				dp_phi = curr_cauchy(eng);
				dp_theta = curr_theta;

				
			}
			

		}
		if (dinfo.down_bias == 1) {
			//cout << "down bias" << endl;
			dp_theta = 1.6; //always move downwards
		}
		else {
			if (dp_theta > 1.6) { //need to go down
				dp_theta = 1.6;
			}
			else if (dp_theta < 1.55) {
				dp_theta = 1.55;
			}
		}

		float dp_x = dist * sin(dp_theta) * cos(dp_phi);
		float dp_y = dist * sin(dp_theta) * sin(dp_phi);
		float dp_z = dist * cos(dp_theta);

		/*float av_x = ((1 - dinfo.rho)*dp_x + dinfo.rho*curr_x) / 2;
		float av_y = ((1 - dinfo.rho)*dp_y + dinfo.rho*curr_y) / 2;
		float av_z = ((1 - dinfo.rho)*dp_z + dinfo.rho*curr_z) / 2;*/
		float av_x = ((1 - dinfo.rho)*dp_x + dinfo.rho*curr_x);
		float av_y = ((1 - dinfo.rho)*dp_y + dinfo.rho*curr_y);
		float av_z = ((1 - dinfo.rho)*dp_z + dinfo.rho*curr_z);

		//calculate new coordinates
		new_x = current_pos.x + av_x;
		new_y = current_pos.y + av_y;
		new_z = current_pos.z - av_z; //because temp_z is going to reflect negative or positive, subtracting here matches whether it goes up or down
		//check seafloor?
		/*vector<int>indexes = calc_r_c_l(new_x, new_y, new_z);
		if (lands->raster_3d[indexes[0]][indexes[1]][0] != nullptr) {
			if (new_z >= lands->raster_3d[indexes[0]][indexes[1]][0]->seafloor_depth) { new_z = lands->raster_3d[indexes[0]][indexes[1]][0]->seafloor_depth; }
		}
		else {
			if (new_z >= current_cell->seafloor_depth) { new_z = current_cell->seafloor_depth; }
		}*/
		if (new_z >= current_cell->seafloor_depth) { new_z = current_cell->seafloor_depth; }
		coordinates.push_back(new_x); coordinates.push_back(new_y); coordinates.push_back(new_z);
		if (calc_distance(new_x, new_y, new_z) == 0) {
			cout << "distance to travel was 0" << endl;
		}
		//cout << "using the current, distance would have been: " << calc_distance(current_pos.x + curr_x, current_pos.y + curr_y, current_pos.z + curr_z) << " but indiv is moving " << calc_distance(new_x, new_y, new_z) << endl;
		//cout << "with SL= " << dinfo.SL << ", rho " << dinfo.rho << ", phi " << dp_phi << ", theta " << dp_theta << ", coordinates are " << new_x << ", " << new_y << ", " << new_z << endl;
		return coordinates;
		
		//if (pstage->s_sett.buffer == 1) { // if indiv uses buffer 
		//	if (dinfo.goal_bias == 1) { //if it locked onto a goal before i want it to continue to that goal, not switch goals
		//		dist = calc_distance(dinfo.goal_focus->midpoint[0], dinfo.goal_focus->midpoint[1], dinfo.goal_focus->midpoint[2]);
		//		//cout << "was locked onto a goal, staying with that goal" << endl;
		//		dp_phi = atan2((dinfo.goal_focus->midpoint[1] - current_pos.y), (dinfo.goal_focus->midpoint[0] - current_pos.x));
		//		dp_theta = acos(-(dinfo.goal_focus->midpoint[2] - current_pos.z) / dist); //acos(z/r)
		//	}
		//	else if (current_cell->bufferpos == -9999) { //if it hasnt been assessed yet,
		//		//cout << "assessing cell buffer" << endl;
		//		current_cell->assess_pos_buffer(pstage->s_sett.buffer_xycells, pstage->s_sett.buffer_zlayers, lands);
		//	}
		//	//if a positive bias was found 

		//	if (current_cell->bufferpos == 1) { //not just open water 
		//		//cout << "found pos buffer" << endl;
		//		for (int p = 0; p < current_cell->buffer_focus.size(); p++) {

		//			if (mod->patchmodel == 0) {//if indiv has checked that cell for settlement before /if that cell is natal
		//				if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), current_cell->buffer_focus[p]->cell_ID)) {
		//					avoid = true;
		//				}
		//			}
		//			else { //if that indiv has tested that patch before / if that patch is natal
		//				if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), current_cell->buffer_focus[p]->patch->patch_ID)) {
		//					avoid = true;
		//				}
		//			}
		//		}
		//		if (avoid == false) {
		//			//cout << "positive goal bias" << endl;
		//			random_shuffle(current_cell->buffer_focus.begin(), current_cell->buffer_focus.end());
		//			temp_goal = current_cell->buffer_focus[0];
		//			dist = calc_distance(temp_goal->midpoint[0], temp_goal->midpoint[1], temp_goal->midpoint[2]);
		//			//cout << "chose pos bias temp goal" << endl;
		//			dinfo.goal_bias == 1; dinfo.goal_focus = temp_goal;

		//			dp_phi = atan2((temp_goal->midpoint[1] - current_pos.y), (temp_goal->midpoint[0] - current_pos.x));
		//			//if (temp_goal->midpoint[2] - current_pos.z > 0) { //need to go down
		//			//	theta = 1.6;
		//			//}
		//			//else { //need to go up
		//			//	theta = 1.55;
		//			//}
		//			dp_theta = acos(-(temp_goal->midpoint[2] - current_pos.z) / dist); //acos(z/r)
		//			if (dist >= dinfo.SL) { //if the pos goal is closer than the SL
		//				dist = dinfo.SL; //it still shouldnt travel past its step length
		//			}
		//		}
		//	}
		//	if (current_cell->bufferpos == 0 || avoid == true) { //open water, use memory and current
		//		//cout << "using memory and current" << endl;
		//		//float dp_theta, dp_phi;
		//		//float curr_theta, curr_phi;
		//		dist = dinfo.SL;
		//		locn prev = memory.front();
		//		//cout << "memory is " << memory.size() << " steps long" << endl;

		//		if (memory.size() >= 2) { //use both memory and current
		//			//current data
		//			//cout << "using prev y " << prev.y << " and prev x " << prev.x << endl;
		//			dp_phi = atan2((current_pos.y - prev.y), (current_pos.x - prev.x)); //angle indiv was travelling at
		//			//cout << "dp_phi is " << dp_phi << endl;
		//			dp_theta = acos((current_pos.z - prev.z) / calc_distance(prev.x, prev.y, prev.z)); //acos(z/r)


		//		}
		//		else { //sample a new direction using current and rho
		//			curr_cauchy = cauchy_distribution<double>(current_cell->cell_angle, dinfo.rho);
		//			dp_phi = curr_cauchy(eng);
		//			dp_theta = curr_theta;

		//			if (dinfo.down_bias == 1) {
		//				//cout << "down bias" << endl;
		//				theta = 1.6; //always move downwards
		//			}
		//			else {
		//				if (dp_theta > 1.6) { //need to go down
		//					dp_theta = 1.6;
		//				}
		//				else if (dp_theta < 1.55) {
		//					dp_theta = 1.55;
		//				}
		//			}
		//		}
		//		
		//	}
		//	float dp_x = dist * sin(dp_theta) * cos(dp_phi);
		//	float dp_y = dist * sin(dp_theta) * sin(dp_phi);
		//	float dp_z = dist * cos(dp_theta);

		//	/*float av_x = ((1 - dinfo.rho)*dp_x + dinfo.rho*curr_x) / 2;
		//	float av_y = ((1 - dinfo.rho)*dp_y + dinfo.rho*curr_y) / 2;
		//	float av_z = ((1 - dinfo.rho)*dp_z + dinfo.rho*curr_z) / 2;*/
		//	float av_x = ((1 - dinfo.rho)*dp_x + dinfo.rho*curr_x) ;
		//	float av_y = ((1 - dinfo.rho)*dp_y + dinfo.rho*curr_y) ;
		//	float av_z = ((1 - dinfo.rho)*dp_z + dinfo.rho*curr_z) ;

		//	//calculate new coordinates
		//	new_x = current_pos.x + av_x;
		//	new_y = current_pos.y + av_y;
		//	new_z = current_pos.z - av_z; //because temp_z is going to reflect negative or positive, subtracting here matches whether it goes up or down
		//	coordinates.push_back(new_x); coordinates.push_back(new_y); coordinates.push_back(new_z);
		//	if (calc_distance(new_x, new_y, new_z) == 0) {
		//		cout << "distance to travel was 0"  << endl;
		//	}
		//	//cout << "using the current, distance would have been: " << calc_distance(current_pos.x + curr_x, current_pos.y + curr_y, current_pos.z + curr_z) << " but indiv is moving " << calc_distance(new_x, new_y, new_z) << endl;
		//	return coordinates;

		//}
		////cout << "phi is " << phi << ", theta is " << theta;
		////cout << "dist is " << dist << endl;
		//if (dist < dinfo.SL) { //if the pos goal is closer than the SL
		//	temp_x = dist * sin(theta) * cos(phi); //make the distance to goal the step length
		//	temp_y = dist * sin(theta) * sin(phi);
		//	if (dinfo.goal_bias == 1) { //if it's heading to a suitable cell
		//		temp_z = dist * cos(theta); //then i dont want to change the goal to be proportionate
		//	}
		//	else {
		//		temp_z = /*propr * */dist * cos(theta);
		//	}
		//}
		//else if (dist >= dinfo.SL) {//it'll take more than one step to reach the pos goal or we dont have a goal and jsut moving by SL	
		//	temp_x = dinfo.SL * sin(theta) * cos(phi); //keep SL 
		//	temp_y = dinfo.SL * sin(theta) * sin(phi);
		//	//for the z velocity, dont take the whole difference, take a bit of it so that it still takes more than one step to get there
		//	if (current_cell->bufferpos == 1 && avoid == false) { //if it's heading to a suitable cell
		//		temp_z = dinfo.SL * cos(theta); //then i dont want to change the goal to be proportionate
		//	}
		//	else {
		//		temp_z = /*propr * */dinfo.SL * cos(theta);
		//	}
		//}

		//float x_actual = temp_x + current_cell->u;
		//float y_actual = temp_y + current_cell->v;
		//float z_actual = temp_z + current_cell->w; //dont know if i need this

		////calculate new coordinates
		//new_x = current_pos.x + x_actual;
		//new_y = current_pos.y + y_actual;
		//new_z = current_pos.z - z_actual; //because temp_z is going to reflect negative or positive, subtracting here matches whether it goes up or down

		//coordinates.push_back(new_x); coordinates.push_back(new_y); coordinates.push_back(new_z);
		//return coordinates;
		////cout << "done with calc_move()" << endl;
	}
}

bool Individual::buffer_capture() {
	//this function will only be called if the individual is competent, uses buffers and has buffer_cap=T
	
	if (current_cell->bufferpos == -9999) { //if it hasnt been assessed yet,
		//cout << "assessing cell buffer" << endl;
		current_cell->assess_pos_buffer(pstage->s_sett.buffer_xycells, pstage->s_sett.buffer_zlayers, lands);
	}

	if (current_cell->bufferpos == 1) { //not just open water 
		//cout << "called buffer_capture, found buffer, ";
		//cout << "found pos buffer" << endl;
		
		//cout << "positive goal bias" << endl;
		random_shuffle(current_cell->buffer_focus.begin(), current_cell->buffer_focus.end());
		dinfo.status = 8; //considering buffer capture (if it doesnt settle, this will be reverted back to dinfo.status=0)
		//if (dinfo.phase == 0) { status = 9; cout << "indiv is passive, "; } //if individual is passive, treat it like forced settlement
		bool cap = check_settle(current_cell->buffer_focus[0], current_cell->buffer_focus[0]->midpoint[0], current_cell->buffer_focus[0]->midpoint[1], current_cell->buffer_focus[0]->midpoint[2]);
		if (cap == false) {
			dinfo.status = 0; //reset the status
			//cout << "indiv wasnt captured" << endl;
			return false;
		}
		else {
			//cout << "indiv was captured" << endl;
			return true;
		}
	}
	else {
		
		return false;
	}
	
}

vector<float> Individual::calc_move2() {
	if (current_cell->cell_speed == -9999) { cout << "indiv " << ID << " started calc_mov with cell speed -9999 at pld " << dinfo.pld << endl; }
	//this is no longer responsible for moving away from unsuitable habitat, it just calculates movements
	vector<float> coordinates;
	float phi, theta;
	cauchy_distribution<double> cauchy, curr_cauchy;
	float new_x, new_y, new_z;
	float temp_x, temp_y, temp_z;
	float dist = 0;
	float propr = lands->get_land_att(8) / lands->get_land_att(7); // dint/resolution

	Cell* temp_goal;
	//cout << "this indiv's SL is " << dinfo.SL << endl;
	bool avoid = false;

	//passive, only dependent on the current, can't use a buffer
	if (dinfo.phase == 0) {
		dist = current_cell->cell_speed;
		theta = acos(current_cell->w / dist);
		//cout << "cell speed is " << dist;
		//try {
		//	if (current_cell->w == -9999) {
		//		throw 9999;
		//	}
		//}
		//catch (int x) {
		//	cout << "settled in an invalid cell" << endl;
		//}

		cauchy = cauchy_distribution<double>(current_cell->cell_angle, -log(.9)); //i am using rho=0.9 for passive individuals
		phi = fmod(cauchy(eng), M2_PI); //random direction in wrapped cauchy distribution, divide by 2pi to put it on the unit circl= "wrap" the cauchy distribution
		//cout << "phi is " << phi << ", theta is " << theta << endl;
		temp_x = dist * sin(theta) * cos(phi);
		temp_y = dist * sin(theta) * sin(phi);
		temp_z = dist * cos(theta);
		//cout << "temp_x: " << temp_x << ", temp_y: " << temp_y << ", temp_z: " << temp_z << endl;
		new_x = current_pos.x + temp_x;
		new_y = current_pos.y + temp_y;
		new_z = current_pos.z - temp_z; //because my z axis is flipped, a negative temp_z denotes increasing depth, so i need to subtract to switch the sign
		if (new_z > lands->get_land_att(14)) { //z_max
			cout << "phi is " << phi << ", theta is " << theta << endl;
			cout << "cell speed is " << dist;
			cout << "current cell habtype is " << current_cell->habtype << endl;
			cout << "starting x " << current_pos.x << " starting y" << current_pos.y << " starting z " << current_pos.z << endl;
			cout << "starting row " << current_cell->r << ", starting col " << current_cell->c << ", starting layer " << current_cell->l << endl;
			cout << "goal x " << new_x << " goal y " << new_y << " goal z " << new_z << endl;
		}
		


		//coordinates.push_back(new_x); coordinates.push_back(new_y); coordinates.push_back(new_z);
		//if (calc_distance(new_x, new_y, new_z) == 0) {
		//	cout << "distance to travel was 0" << endl;
		//}
		////cout << "using the current, distance is: " << calc_distance(new_x, new_y, new_z)  << endl;

		//return coordinates;
	}



	//active, dependent on current but can respond to environmental cues, can use a buffer
	else {
		//where is the current going?
		curr_cauchy = cauchy_distribution<double>(current_cell->cell_angle, .9); //add some stochasticity
		float curr_phi = curr_cauchy(eng);
		float curr_theta = acos(current_cell->w / current_cell->cell_speed);
		float curr_x = current_cell->cell_speed * sin(curr_theta) * cos(curr_phi);
		float curr_y = current_cell->cell_speed * sin(curr_theta) * sin(curr_phi);
		float curr_z = current_cell->cell_speed * cos(curr_theta);

		//where is the individual going?
		float dp_theta, dp_phi;

		if (dinfo.goal_bias == 1) { //if it locked onto a goal before i want it to continue to that goal, not switch goals
			dist = calc_distance(dinfo.goal_focus->midpoint[0], dinfo.goal_focus->midpoint[1], dinfo.goal_focus->midpoint[2]);
			//cout << "was locked onto a goal, staying with that goal" << endl;
			dp_phi = atan2((dinfo.goal_focus->midpoint[1] - current_pos.y), (dinfo.goal_focus->midpoint[0] - current_pos.x));
			dp_theta = acos(-(dinfo.goal_focus->midpoint[2] - current_pos.z) / dist); //acos(z/r)
			if (dist >= dinfo.SL) { //if the pos goal is closer than the SL
				dist = dinfo.SL; //it still shouldnt travel past its step length
			}
		}
		else if (dinfo.comp == true && pstage->s_sett.buffer == 1) {//step 1: explore surroundings, they can only seek suitable habitat if they are competent
			if (current_cell->bufferpos == -9999) { //if it hasnt been assessed yet,
				//cout << "assessing cell buffer" << endl;
				current_cell->assess_pos_buffer(pstage->s_sett.buffer_xycells, pstage->s_sett.buffer_zlayers, lands);
			}
			if (current_cell->bufferpos == 1) { //not just open water 
				//cout << "found pos buffer" << endl;
				for (int p = 0; p < current_cell->buffer_focus.size(); p++) {

					if (mod->patchmodel == 0) {//if indiv has checked that cell for settlement before /if that cell is natal
						if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), current_cell->buffer_focus[p]->cell_ID)) {
							avoid = true;
						}
					}
					else { //if that indiv has tested that patch before / if that patch is natal
						if (count(dinfo.tested_sett.begin(), dinfo.tested_sett.end(), current_cell->buffer_focus[p]->patch->patch_ID)) {
							avoid = true;
						}
					}
				}
				if (avoid == false) {
					//cout << "positive goal bias" << endl;
					random_shuffle(current_cell->buffer_focus.begin(), current_cell->buffer_focus.end());
					temp_goal = current_cell->buffer_focus[0];
					dist = calc_distance(temp_goal->midpoint[0], temp_goal->midpoint[1], temp_goal->midpoint[2]);
					//cout << "chose pos bias temp goal" << endl;
					dinfo.goal_bias == 1; dinfo.goal_focus = temp_goal;

					dp_phi = atan2((temp_goal->midpoint[1] - current_pos.y), (temp_goal->midpoint[0] - current_pos.x));
					//if (temp_goal->midpoint[2] - current_pos.z > 0) { //need to go down
					//	theta = 1.6;
					//}
					//else { //need to go up
					//	theta = 1.55;
					//}
					dp_theta = acos(-(temp_goal->midpoint[2] - current_pos.z) / dist); //acos(z/r)
					if (dist >= dinfo.SL) { //if the pos goal is closer than the SL
						dist = dinfo.SL; //it still shouldnt travel past its step length
					}
				}
			}
		}
		if (dinfo.comp == false || pstage->s_sett.buffer == 0 || current_cell->bufferpos == 0 || avoid == true) {//open water, use memory and current
				//cout << "using memory and current" << endl;
				//float dp_theta, dp_phi;
				//float curr_theta, curr_phi;
			dist = dinfo.SL;
			locn prev = memory.front();
			//cout << "memory is " << memory.size() << " steps long" << endl;

			if (memory.size() >= 2) { //use both memory and current
				//current data
				//cout << "using prev y " << prev.y << " and prev x " << prev.x << endl;
				dp_phi = atan2((current_pos.y - prev.y), (current_pos.x - prev.x)); //angle indiv was travelling at
				//cout << "dp_phi is " << dp_phi << endl;
				dp_theta = acos((current_pos.z - prev.z) / calc_distance(prev.x, prev.y, prev.z)); //acos(z/r)

			}
			else { //sample a new direction using current and rho
				curr_cauchy = cauchy_distribution<double>(current_cell->cell_angle, dinfo.rho);
				dp_phi = curr_cauchy(eng);
				dp_theta = curr_theta;


			}


		}
		if (dinfo.down_bias == 1) {
			//cout << "down bias" << endl;
			dp_theta = 1.6; //always move downwards
		}
		else {
			if (dp_theta > 1.6) { //need to go down
				dp_theta = 1.6;
			}
			else if (dp_theta < 1.55) {
				dp_theta = 1.55;
			}
		}

		float dp_x = dist * sin(dp_theta) * cos(dp_phi);
		float dp_y = dist * sin(dp_theta) * sin(dp_phi);
		float dp_z = dist * cos(dp_theta);

		/*float av_x = ((1 - dinfo.rho)*dp_x + dinfo.rho*curr_x) / 2;
		float av_y = ((1 - dinfo.rho)*dp_y + dinfo.rho*curr_y) / 2;
		float av_z = ((1 - dinfo.rho)*dp_z + dinfo.rho*curr_z) / 2;*/
		//float av_x = ((1 - dinfo.rho)*dp_x + dinfo.rho*curr_x);
		//float av_y = ((1 - dinfo.rho)*dp_y + dinfo.rho*curr_y);
		//float av_z = ((1 - dinfo.rho)*dp_z + dinfo.rho*curr_z);
		float av_x = dp_x + curr_x;
		float av_y = dp_y + curr_y;
		float av_z =  dp_z + curr_z;

		//calculate new coordinates
		new_x = current_pos.x + av_x;
		new_y = current_pos.y + av_y;
		new_z = current_pos.z - av_z; //because temp_z is going to reflect negative or positive, subtracting here matches whether it goes up or down

		//check seafloor?
		/*if (new_z > current_cell->seafloor_depth) {
			new_z = current_cell->seafloor_depth - 1;
		}*/
		//vector<int> indexes = calc_r_c_l(new_x, new_y, new_z);
		/*if (indexes[0] > lands->get_land_att(1) || indexes[0] < 0) {
			cout << "indiv " << ID << " with SL " << dinfo.SL << " and comp time " << dinfo.comp_time << " calculated an invalid move at " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << endl;
			cout << "now: " << current_pos.x << ", " << current_pos.y << ", " << current_pos.z << endl;
			cout << "curr: " << curr_x << ", " << curr_y << ", " << curr_z << endl;
			cout << "dp: " << dp_x << ", " << dp_y << ", " << dp_z << endl;
			cout << "av: " << av_x << ", " << av_y << ", " << av_z << endl;
			cout << "new: " << new_x << ", " << new_y << ", " << new_z << endl;
		}*/

		/*vector<int>indexes = calc_r_c_l(new_x, new_y, new_z);
		if (lands->raster_3d[indexes[0]][indexes[1]][0] != nullptr) {
			if (new_z >= lands->raster_3d[indexes[0]][indexes[1]][0]->seafloor_depth) { new_z = lands->raster_3d[indexes[0]][indexes[1]][0]->seafloor_depth; }
		}
		else {
			if (new_z >= current_cell->seafloor_depth) { new_z = current_cell->seafloor_depth; }
		}*/
		//if (new_z >= current_cell->seafloor_depth) { new_z = current_cell->seafloor_depth; }
		//coordinates.push_back(new_x); coordinates.push_back(new_y); coordinates.push_back(new_z);
		//if (calc_distance(new_x, new_y, new_z) == 0) {
		//	cout << "distance to travel was 0" << endl;
		//}
		////cout << "using the current, distance would have been: " << calc_distance(current_pos.x + curr_x, current_pos.y + curr_y, current_pos.z + curr_z) << " but indiv is moving " << calc_distance(new_x, new_y, new_z) << endl;
		////cout << "with SL= " << dinfo.SL << ", rho " << dinfo.rho << ", phi " << dp_phi << ", theta " << dp_theta << ", coordinates are " << new_x << ", " << new_y << ", " << new_z << endl;
		//return coordinates;

	}

	//adjust destination for buoyancy and diel vertical migration
	//check buoyancy
	if (dinfo.buoy_min != -9) { //if there are buoyancy limits

		if (new_z < dinfo.buoy_min) {
			//cout << "hit upper buoy limit" << endl;
			new_z = dinfo.buoy_min + 1;
		}
		else if (new_z > dinfo.buoy_max) {
			//cout << "hit lower buoy limit" << endl;
			new_z = dinfo.buoy_max - 1; //i want it to stay just above the buoyancy maximum
		}
	}
	else if (dinfo.buoy_min == -9 && new_z < 0) { //if indiv has hit the surface and there is no buoyancy limit

		//cout << "hit surface" << endl;
		new_z = 1;
		//i don't update position or anything because i want to know if that new cell is water 
	}

	//the case of an indiv going too low when buoy_max==-9 is covered in the check_cell function
	if (dinfo.diel_vert == 1 && dvm == true) { //check that it is old/big enough to undergo diel vertical migration

		int min_z, max_z, buoy_layer;
		//int sf_layer = floor((current_cell->seafloor_depth - lands->get_land_att(13)) / lands->get_land_att(8));
		//float sf = 0;
		//if (lands->raster_3d[current_cell->r][current_cell->c][sf_layer]->habtype != 0) { //it's solid, that means seafloor depth was less than the halfway depth of the cell, so whole cell is solid
		//	sf = lands->raster_3d[current_cell->r][current_cell->c][sf_layer]->min_depth - 1;
		//}
		//else { //water, so can use seafloor depth as a measure
		//	sf = current_cell->seafloor_depth - 1;
		//}
		//check the disp time. remember that it was released at night and there is a 12 hour day night cycle
		//so disp-time/12, if it's odd then it's day time? even then nighttime?
		//cout << "indiv has had " << dinfo.disp_time << " hours of dispersal" << endl;
		int portion = floor(dinfo.disp_time / 12);
		if (portion % 2 == 0) { //even=night

			min_z = dinfo.buoy_min;
			max_z = dinfo.buoy_min + dinfo.dv_range;
			//if (max_z >= sf) { //if this takes it too deep
			//	max_z = sf;
			//}
			if (max_z >= current_cell->seafloor_depth) { //if this takes it too deep
				max_z = current_cell->seafloor_depth - 1;
			}


			//cout << "it's night so indiv between " << min_z <<" and " << max_z << endl;
		}
		else { //odd=day
			max_z = dinfo.buoy_max;
			min_z = dinfo.buoy_max - dinfo.dv_range;
			//cout << "buoy max is " << dinfo.buoy_max << " seafloor depth is " << current_cell->seafloor_depth;
			if (dinfo.buoy_max >= current_cell->seafloor_depth) { //if the seafloor is shallower than the buoyancy max
				min_z = current_cell->seafloor_depth- dinfo.dv_range;
				max_z = current_cell->seafloor_depth - 1;
				if (min_z < 0) { min_z = 0; }
			}
			//cout << "it's day so indiv between " << min_z <<" and " << max_z << endl;
		}
		uniform_int_distribution<int>dv_dist(min_z, max_z);
		new_z = dv_dist(eng);
		//if (portion % 2 != 0) { cout << " z is " << z << endl; }
		//if it's daytime, cant go above dinfo.buoy_min or below dinfo.buoy_min +dv_range
		//if it's nighttime, cant go above dinfo.buoy_max-dv_range or below dinfo.buoy_max		
	}

	vector<int> indexes = calc_r_c_l(new_x, new_y, new_z);
	/*if (indexes[0] > lands->get_land_att(1) || indexes[0] < 0 ) {
		cout << "indiv " << ID << " with SL " << dinfo.SL << " and comp time " << dinfo.comp_time << " calculated an invalid move at " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << endl;
	}*/
	coordinates.push_back(new_x); coordinates.push_back(new_y); coordinates.push_back(new_z);
	if (calc_distance(new_x, new_y, new_z) == 0) {
		cout << "distance to travel was 0" << endl;
	}
	//cout << "using the current, distance is: " << calc_distance(new_x, new_y, new_z)  << endl;

	return coordinates;
}
//mat333 Individual::calc_dp() {
//	mat333 dp;
//	float phi;
//
//	locn prev = memory.front(); //even if memory==1, this holds the previous step
//	phi = atan2(current_pos.y - prev.y, current_pos.x - prev.x); //calculate the angle the indiv has been travelling in
//
//	if (phi > 0) { fill_mat333(&dp, dinfo.dp, phi); }
//	else{ 	//initialise structs to -9
//		for (int i = 0; i < 3; i++) {
//			for (int j = 0; j < 3; j++) {
//				for (int k = 0; k < 3; k++) {
//					dp.neigh[i][j][k] = 1;
//				}
//			}
//		}
//	}
//	return dp;
//}

bool Individual::transfer_to_settle() { //this function only returns true if it is still in transit! all other cases, it returns false (ie death or settlement)
	//have all indivs take at least one step
	//cout << "taking steps for indiv " << ID << endl;
	//cout << "starting cell for indiv " << ID << " is " << current_cell->r << ", " << current_cell->c << ", " << current_cell->l << endl;
	//dinfo.per_step_disp_dist = 0; //reset perstep distance so there is no remnant of anything (for whatever reason)
	float before = dinfo.disp_distance;
	//cout << "for indiv " << ID << " before disp, distance was " << dinfo.disp_distance << endl;
	bool captured = false;
	if (pstage->s_sett.buffer_capture == 1 && dinfo.comp == true ) {
		captured = buffer_capture();
	}
	if (captured == true) {
		//cout << "individual has been captured by settlement buffer" << endl;
		return false; //because it has settled
	}
	vector<float> new_coords = calc_move2();
	//cout << "starting traverse" << endl;
	//cout << "indiv " << ID << " has SL " << dinfo.SL << endl;
	/*if (dinfo.SL == 1) { passive_traverse(new_coords[0], new_coords[1], new_coords[2]); }
	else if (dinfo.SL > 1) { traverse(new_coords[0], new_coords[1], new_coords[2]); }
	else {
		cout << "indiv " << ID << "has weird SL: " << dinfo.SL << endl;
	}*/
	traverse2(new_coords[0], new_coords[1], new_coords[2]);
	//traverse3(new_coords[0], new_coords[1], new_coords[2]);
	//cout << "for indiv " << ID << " per_step was " << dinfo.per_step_disp_dist << endl;
	if (dead == true) { return false; } //i am returning true so that it will be subtracted from the matrix population and not considered in transit anymore
	
	//if (dinfo.per_step_disp_dist <= 0) { cout << "after traverse, indiv " << ID << "has per step of " << dinfo.per_step_disp_dist << endl; }
	//cout << "done with traverse" << endl;
	//if (dinfo.per_step_disp_dist == 0 && dead == false) { cout << "indiv " << ID << " after traverse, per step is 0, status is " << status << endl; }
	dinfo.disp_distance += dinfo.per_step_disp_dist;
	//if (dinfo.disp_distance == 0) { cout << "indiv" << ID << " after traverse, disp distance=0" << endl; }
	dinfo.per_step_disp_dist = 0;
	//cout << "for indiv " << ID << "disp dist is now " << dinfo.disp_distance << endl;
	float after = dinfo.disp_distance;
	//if(after<=before){ cout << "indiv" << ID << " after traverse, disp distance is less than before " << endl; }
	if (dead == false) { dispersal_mort(); } //if it hasn't been absorbed by a boundary, apply dispersal mortality
	
	if (dead == true) { return false; }
	if (sett == true) { return false; }
	else { return true; } //this means it is still in transit!
	
}


//survival and development/aging
bool Individual::survival() { //this assumes a stage-structured population, because in non-overlapping gens models, all indivs of previous gen die
	/*default_random_engine unif_generator;
	uniform_real_distribution<float> surv_dist(0, 1);*/
	float ran_surv = unif_dist(eng); //sample randomly from a uniform distribution

	if (age > pop_dyn->max_age) { //if its too old, it dies
		dead = true; status = 7; //cout << "too old ";
		return false;
	}

	float surv;

	if (sex == 0) { //male
		if (pop_dyn->survdensdep == 1) { //if survival is density dependent
			surv = pstage->m_surv*exp(-pop_dyn->survdenscoef * pop_dyn->bc * get_localN());
		}
		else { //if not
			surv = pstage->m_surv;
		}
	}

	else { //female
		if (pop_dyn->survdensdep == 1) { //if survival is density dependent
			surv = pstage->f_surv * exp(-pop_dyn->survdenscoef * pop_dyn->bc * get_localN());
		}
		else { //if not
			surv = pstage->f_surv;
		}
	}
	//if (ran_surv < surv) { /*cout << "survives ";*/ return true; } //individual survives
	//else { dead = true; status = 6; /*cout << "died ";*/ return false; } //if random number is greater than survival rate, individual die

	if (ran_surv >= surv) {
		dead = true; status = 6; return false;
	}
	else {
		return true;
	}
	
	
}
bool Individual::fishing_mort() {
	//if (pop_dyn->fishing_mort != -9) {
	if (pstage->stage == 0) {
		return false;
	}
	float ran_fish = unif_dist(eng);
	//cout << "ran fish is " << ran_fish << " and mort prob is " << pop_dyn->fishing_mort;
	if (ran_fish <= pop_dyn->fishing_mort) {
			
		dead = true; status = 11; //cout << ", died, status is " << status;
		return true;
	}
	else {
		//cout << ", escaped, status is " << status  ;
		return false;
	}
	//}
	
	
}


//this is for YEARLY development to the next stage in a stage structured population
void Individual::stage_development() {//this assumes a stage-structured population, because in non-overlapping gens models, all indivs of previous gen die
	if (status < 0) { cout << "in stage_development, status is weird" << endl; }
	//if disp is stage-dependent, then it can't have iiv or evo, so everything can be overwritten without worries, as all values will be the same anyway
	//if disp is stage-independent and there might be iiv or evo, need to be careful what i overwrite when I update stage information
	//also need to be careful of any calculations that need to happen, ie growth 
	if (pstage->stage == mod->stages-1) { //if individual is already in the last stage, dont bother with stage dev
		return;
	}

	//step1: check if they are old enough to transition
	float ran_dev = unif_dist(eng); //sample randomly from a uniform distribution
	stage_infob* temp_si = (*subpop_stages)[pstage->stage +1]; //get min age of next stage
	float dev = 0;
	if (sex == 0 && age >= temp_si->m_min_age) { //male; if they would be old enough to transition
		if (pop_dyn->devdensdep == 1) { //if density dependence is true
			dev = pstage->m_trans * exp(-pop_dyn->devdenscoef * pop_dyn->bc * get_localN());
		}
		else {
			dev = pstage->m_trans;
		}
	}
	else if (sex == 1 && age >= temp_si->f_min_age) { //female; 
		if (pop_dyn->devdensdep == 1) {//if density dependence is true
			dev = pstage->f_trans * exp(-pop_dyn->devdenscoef * pop_dyn->bc * get_localN());
		}
		else {
			//cout << "indiv's age is " << age << ", min age to transition is " << temp_si->f_min_age;
			dev = pstage->f_trans;
		}
	}

	if (ran_dev < dev) { //if it develops
		
		pstage = temp_si; //update stage information
		if (emig == false && pstage->s_emig.stagedep==1 ) {//if stage dependent, then no iiv/evo so can just take population values of that stage
			
			set_emig_info(false);
			set_trans_info(false);
			set_settle_info(false);
			
			//if stage independent so values have already been calculated at inheritance, so don't need to update dispersal info
		}
		else if (emig == false && pstage->s_emig.stagedep == 0 && pstage->stage == pstage->s_emig.emigstage) { //if it's stage indep and this indiv is now in the stage that emigrates
			////emgiration
			//dinfo.ep_evo = pstage->s_emig.ep_evo;
			//dinfo.ep_mut = pstage->s_emig.ep_mut;
			//dinfo.ep_mut_size = pstage->s_emig.ep_mut_size;
			//dinfo.alpha_mut_size = pstage->s_emig.alpha_mut_size;
			//dinfo.beta_mut_size = pstage->s_emig.beta_mut_size;
			////transfer
			//dinfo.active_evo = pstage->s_trans.active_evo;
			//dinfo.active_mut = pstage->s_trans.active_mut;
			//dinfo.active_mut_size = pstage->s_trans.active_mut_size;
			//dinfo.buoy_min = pstage->s_trans.buoy_min;
			//dinfo.buoy_max = pstage->s_trans.buoy_max;
			//dinfo.diel_vert = pstage->s_trans.diel_vert;
			//dinfo.dv_range = pstage->s_trans.dv_range;
			//dinfo.min_dv_size = pstage->s_trans.min_dv_size;
			//dinfo.min_dv_time = pstage->s_trans.min_dv_time;
			//dinfo.down_bias = pstage->s_trans.down_bias;
			//dinfo.memory = pstage->s_trans.memory;
			//dinfo.grow_info.method = pstage->s_trans.pop_grow.method;
			//dinfo.grow_info.growth_evo = pstage->s_trans.pop_grow.growth_evo;
			//dinfo.grow_info.growth_mut = pstage->s_trans.pop_grow.growth_mut;
			//if (dinfo.min_dv_size == -9) { cout << "min dv_size not set right at stage dev" << endl; }
			////settlement
			//dinfo.S0_evo = pstage->s_sett.S0_evo;
			//dinfo.S0_mut = pstage->s_sett.S0_mut;
			//dinfo.S0_mut_size = pstage->s_sett.S0_mut_size;
			//dinfo.alphaS_mut_size = pstage->s_sett.alphaS_mut_size;
			//dinfo.betaS_mut_size = pstage->s_sett.betaS_mut_size;
			//dinfo.comp_evo = pstage->s_sett.comp_evo;
			//dinfo.comp_mut = pstage->s_sett.comp_mut;
			//dinfo.comp_mut_size = pstage->s_sett.comp_mut_size;

			set_emig_info(false);
			set_trans_info(false);
			set_settle_info(false);
			disp_gene_expression(); //give it its inherited genes from parents
		}

		
	}

}

//this is for growth and development DURING DISPERSAL, where changes in behaviour happen depending on size
void Individual::larval_growth() {
	//cout << "calling larval growth, ";
	float r;
	

	if (dinfo.grow_info.method==0) { //linear
		if (mod->temp_dep == 1) {
				r = dinfo.grow_info.b + dinfo.grow_info.m*current_cell->cell_temp; //r=0.02+0.03*T comes from Lett et al 2008
		}
		else {
			r = dinfo.grow_info.m;
		}
		size= (r/24)*dinfo.disp_time + size_at_birth; //the input data for rate is given in cm/d so i need to make it cm/hr
	}
	else if (dinfo.grow_info.method==1) { //gompertz
		//int t = dinfo.pld - dinfo.disp_time;
		//size = size_at_birth + dinfo.grow_info.Linf * (exp(-exp(-dinfo.grow_info.G_K * (-t)))); //here, t is the time until PLD is reached
		int t = floor(dinfo.disp_time/24); //the equation works in days, so need to put disp time back into days
		//cout << "t in larval growth is " << t << " and disp time is " << dinfo.disp_time << ", Ti is " << dinfo.grow_info.Ti ;
		size = size_at_birth + dinfo.grow_info.Linf * (exp(-exp(-dinfo.grow_info.G_K * (t-dinfo.grow_info.Ti)))); //here, t is the time until PLD is reached
		//cout << " indiv " << ID << "is " << size << " mm big " << endl;
	}
	else if (dinfo.grow_info.method == 2) { //von bertalanffy

	}
	
	//after size has been calculated, are there any behavioural changes?
	//diel vertical migration

	if (dinfo.diel_vert==1 && dinfo.min_dv_size !=-9 && size >= dinfo.min_dv_size) {
		//if (dvm == false) { cout << "indiv " << ID << " is now undergoing diel vertical migration at size " << size << endl; }
		dvm = true;
			
	}
	calc_SL_rho();
	
	//the only time calc_SL_rho isnt called is if individuals are below the minimum active size because then SL remains =1 and rho remains =.9
	if (dinfo.mode !=0 && size >= dinfo.min_active_size) { //as long as it's not passive
		if (dinfo.dv_active == 0) {
			dvm = false;
		}
		dinfo.phase = 1; //transition to active phase
		//cout << "called calc_SL_rho, SL is now: " << dinfo.SL << " and rho is " << dinfo.rho << endl;
	
	}
	if (dinfo.comp_size != -9 && size >= dinfo.comp_size) { //if they are the minimum size they need to reach before they can settle
			//if (dinfo.comp != true) { cout << "individual " << ID << " is now competent for settlement at size " << size << endl; } //if they have just now become competent
		/*calc_SL_rho();*/ //if this hasnt been called above, that means that min_active_size==-9 so rho and SL never change from 0.9 and 1 respectively 
		dinfo.comp = true;

	}
	
	
	
	//by the end of this function, the indiv should have grown and any behavioural changes should be made
}