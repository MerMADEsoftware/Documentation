#include "classes.h"
//default constructor
Subpopulation::Subpopulation() {

}

Subpopulation::~Subpopulation()
{
}

//constructor for cell-based subpop
Subpopulation::Subpopulation(Cell* c, Model* m, Landscape* la, pop_info* pd, matrix* pm, /*emig_info* e, trans_info* t, sett_info*s,*/ vector<stage_infob*> st) {
	cell_num = c; //assign the pointer to that cell to the subpop
	mod = m;
	lands = la;
	pop_dyn = pd;
	pop_matrix = pm;
	/*pop_einfo = e;
	pop_tinfo = t;
	pop_sinfo = s;*/
	pop_stageinfo = st;
	size = 0;
}

//constructor for patch-based subpop
Subpopulation::Subpopulation(Patch* p, Model* m, Landscape* la, pop_info* pd, matrix* pm, /*emig_info* e, trans_info* t, sett_info*s,*/ vector<stage_infob*> st) {
	patch_num = p; //assign the pointer to that cell to the subpop
	mod = m;
	lands = la;
	pop_dyn = pd;
	pop_matrix = pm;
	/*pop_einfo = e;
	pop_tinfo = t;
	pop_sinfo = s;*/
	pop_stageinfo = st;
	size = 0;
}

//general use functions
int Subpopulation::assign_sex() {
	int sex;
	if (mod->repro == 0) { sex = 1; } //if it's an asexual model, sex of individual is female
	else {
		/*default_random_engine s_generator;
		uniform_real_distribution<double> uni_dist(0.0, 1.0);*/
		double sgen = unif_dist(eng);
		if (sgen < mod->propmales) { //assign a certain proportion the male sex
			sex = 0;
		}
		else { sex = 1; }
	}
	return sex;
}

void Subpopulation::leave(Individual* pind) { //takes care of updating subpop numbers and removing individuals from vectors
	//size--; //reduce overall size
	if (mod->repro > 0) {//if its a sexual model 
		if (mod->stagestruct == 1) { //stagestructured model
			if (pind->sex == 0 ) {//and the individual is a reproductively mature male
				nmales--;
				//rep_males.erase(remove_if(rep_males.begin(), rep_males.end(), [](Individual* i, Individual * pind) {return i==pind; }), rep_males.end());
				if (pind->pstage->m_fec > 0) {
					for (int m = 0; m < rep_males.size(); m++) {
						if (rep_males[m] == pind) { rep_males.erase(rep_males.begin() + (m - 1)); break; } //once you find it, stop searching
					}
				}
			}
			else if (pind->sex == 1) {
				nfemales--;
				//rep_females.erase(remove_if(rep_females.begin(), rep_females.end(), [](Individual* i, Individual * pind) {return i == pind; }), rep_females.end());
				if (pind->pstage->f_fec > 0) {
					for (int m = 0; m < rep_females.size(); m++) {
						if (rep_females[m] == pind) { rep_females.erase(rep_females.begin() + (m - 1)); break; } //once you find it, stop searching
					}
				}

			}
		}
		else { //nonoverlapping generations
			if (pind->sex == 0) {
				nmales--;
				for (int m = 0; m < rep_males.size(); m++) {
					if (rep_males[m] == pind) { rep_males.erase(rep_males.begin() + (m - 1)); break; } //once you find it, stop searching
				}
			}
			else {
				nfemales--;
				for (int m = 0; m < rep_females.size(); m++) {
					if (rep_females[m] == pind) { rep_females.erase(rep_females.begin() + (m - 1)); break; } //once you find it, stop searching
				}
			}
		}
		
	}
	//don't need to do anything else for an asexual model because they are only present in the indivs_vector, which will be cleared
	if (mod->patchmodel == 1) { pind->current_patch = nullptr; }
	pop_matrix->matrix_indivs.push_back(pind); //add them to the vector of indivs that undergo transfer phase
	pop_matrix->total++;
}


void Subpopulation::enter(Individual* pind) { //does the opposite of leave()
	//cout << "entering subpop " << patch_num->patch_ID << endl;
	size++;
	//cout << "natal patch was " << pind->natal_patch;
	if (mod->stagestruct == 1) { //if the model is stagestructured
		if (mod->repro > 0) { //if it's a sexual model
			if (pind->sex == 0 && pind->pstage->m_fec > 0) { nmales++; rep_males.push_back(pind); } //and the individual is a reproductively mature male
			else if (pind->sex == 1 && pind->pstage->f_fec > 0) { nfemales++; rep_females.push_back(pind); } //or female
		}
		//nothing else needed for an asexual model because there are no possible mating requirements
	}
	else {//non overlapping generations
		if (mod->repro > 0) {
			if (pind->sex == 0) { nmales++; rep_males.push_back(pind); }
			else { nfemales++; rep_females.push_back(pind); }
		}
		//nothing else needed for an asexual model because there are no possible mating requirements
	}
	surv_indivs.push_back(pind); //add it to the vector of surviving individuals
	if (mod->patchmodel == 1) {
		pind->current_patch = patch_num; //update location info
	}
	else {
		pind->current_cell = cell_num;
	}
	pop_matrix->total--; //subtract one from the matrix population
	//cout << " and is now " << pind->current_patch->patch_ID << endl;
}

Population::stage_infob* Subpopulation::get_stage_info(int s) {
	return(pop_stageinfo[s]);
}
//initialisation functions
void Subpopulation::init_age(Individual* pind, int s) {
	//default_random_engine a_generator; //to assign age
	
	if (mod->initage == 0) { //if age initialisation is at lowest possible
		if (pind->sex == 0) { pind->age = pind->pstage->m_min_age; } //if it's male
		else { pind->age = pind->pstage->f_min_age; } //if it's female
	} 
	else {					//if it is randomised
		int agen; float min, max;
		if (mod->stagestruct == 1) { //if it's a stage structured model
			if (pind->sex == 0) { //if it's male
				min = pind->pstage->m_min_age; max = pind->pstage->m_max_age;
			}
			else { //if it's female
				min = pind->pstage->f_min_age; max = pind->pstage->f_max_age;
			}
		}
		else { //if its non overlapping generations
			min = 0; max = pop_dyn->max_age; //population max age
		}
		uniform_real_distribution<float> age_dist(min, max);
		agen = age_dist(eng);

		pind->age = floor(agen); //make age a whole number
	}
}

void Subpopulation::init_size(Individual* pind) {
	normal_distribution<float> size_dist;
	float n_size;
	if (pop_dyn->iiv_juv_size == 1) { //if there is inter-individual variability in juvenile size
		size_dist = normal_distribution<float>(pop_dyn->juv_size, pop_dyn->juv_size_sd);
		n_size = size_dist(eng);
	}
	else {
		n_size = pop_dyn->juv_size;
	}
	pind->size_at_birth = n_size; //size at birth doesnt get updated, so we can record it
	pind->size = n_size; //this will get updated as it grows throughout dispersal
}

void Subpopulation::init_position(Individual* pind) {
	uniform_int_distribution<int> distribution(0, 1); //i know this is duplicating the global uniform dist, but since i need it to be modified, i have to put this here
	int index; 

	if (mod->patchmodel == 0) { //if cell-based
		pind->current_cell = cell_num; //
	}
	else { //patch-based
		pind->natal_patch = patch_num->patch_ID;
		if (mod->seedtype == 0 || mod->seedtype==2) { //if not initialised from species distribution
			distribution = uniform_int_distribution<int>(0, patch_num->indepth_cells.size() - 1); //choose position from cells in the right depth range
			index = distribution(eng);
			pind->current_cell = patch_num->indepth_cells[index];//assign individual that cell
			if (pind->current_cell == nullptr) { cout << "assigning cell didn't work!" << endl; }
		}
		else if(mod->seedtype==1) { //if initialised from species dist
			distribution= uniform_int_distribution<int>(0, patch_num->in_spdist_cells.size() - 1); //choose position from cells in species distribution
			index = distribution(eng);
			pind->current_cell = patch_num->in_spdist_cells[index];//assign individual that cell

		}
		pind->current_patch = patch_num; //update cell's patch pointer
	}
	
	//struct mat_index { int row, col, layer; };
	//
	//std::deque<mat_index> possibilities;

	//if ((pind->current_cell->r - 1) >=0 && lands->raster_3d[pind->current_cell->r - 1][pind->current_cell->c][pind->current_cell->l]!=nullptr && lands->raster_3d[pind->current_cell->r -1][pind->current_cell->c ][pind->current_cell->l]->habtype == 0) { //north
	//	mat_index ind;
	//	ind.row = pind->current_cell->r - 1; ind.col = pind->current_cell->c; ind.layer=pind->current_cell->l; //record that location
	//	possibilities.push_back(ind); //and add it to the list of possibilities
	//}

	//if ((pind->current_cell->r +1) < lands->get_land_att(1) && lands->raster_3d[pind->current_cell->r + 1][pind->current_cell->c][pind->current_cell->l]!=nullptr && lands->raster_3d[pind->current_cell->r + 1][pind->current_cell->c][pind->current_cell->l]->habtype == 0) { //south
	//	mat_index ind;
	//	ind.row = pind->current_cell->r + 1; ind.col = pind->current_cell->c; ind.layer=pind->current_cell->l; //record that location
	//	possibilities.push_back(ind); //and add it to the list of possibilities
	//}

	//if ((pind->current_cell->c - 1) >= 0  && lands->raster_3d[pind->current_cell->r][pind->current_cell->c - 1][pind->current_cell->l]!=nullptr && lands->raster_3d[pind->current_cell->r][pind->current_cell->c-1][pind->current_cell->l]->habtype == 0) { //west
	//	mat_index ind;
	//	ind.row = pind->current_cell->r ; ind.col = pind->current_cell->c-1; ind.layer=pind->current_cell->l; //record that location
	//	possibilities.push_back(ind); //and add it to the list of possibilities
	//}

	//if ((pind->current_cell->c + 1) < lands->get_land_att(2) && lands->raster_3d[pind->current_cell->r][pind->current_cell->c + 1][pind->current_cell->l]!=nullptr && lands->raster_3d[pind->current_cell->r][pind->current_cell->c+1][pind->current_cell->l]->habtype == 0) { //east
	//	mat_index ind;
	//	ind.row = pind->current_cell->r ; ind.col = pind->current_cell->c+1; ind.layer=pind->current_cell->l; //record that location
	//	possibilities.push_back(ind); //and add it to the list of possibilities
	//}

	//if ((pind->current_cell->l - 1) >= 0 && lands->raster_3d[pind->current_cell->r][pind->current_cell->c][pind->current_cell->l - 1]!=nullptr && lands->raster_3d[pind->current_cell->r][pind->current_cell->c][pind->current_cell->l-1]->habtype == 0) { //up
	//	mat_index ind;
	//	ind.row = pind->current_cell->r; ind.col = pind->current_cell->c; ind.layer=pind->current_cell->l-1; //record that location
	//	possibilities.push_back(ind); //and add it to the list of possibilities
	//}
	//i'm not going to allow them to go down
	if(pind->current_cell->possibilities.size()==0){
		pind->current_cell->find_ow_neigh(lands);
	}
	//shuffle the possibilities so that when i pick one it is random
	random_shuffle(pind->current_cell->possibilities.begin(), pind->current_cell->possibilities.end());
	if (pind->current_cell->possibilities.size() == 0) { cout << "there were no open water neighbours...hab type is " << cell_num->habtype << endl; }
	Cell::mat_index ind = pind->current_cell->possibilities[0];
	if (lands->raster_3d[ind.row][ind.col][ind.layer] == nullptr) {
		cout << "one of the possibilities in init_pos is a nullptr somehow!" << endl;
	}
	
	uniform_real_distribution<float> x_dist(pind->current_cell->cell_corners.llc[0]+1, pind->current_cell->cell_corners.lrc[0]-1); //dont allow it to be given the actual min and max
	uniform_real_distribution<float> y_dist(pind->current_cell->cell_corners.llc[1] +1, pind->current_cell->cell_corners.ulc[1] -1);
	uniform_real_distribution<float> z_dist(pind->current_cell->min_depth +1, pind->current_cell->max_depth -1);
	float x_coord;
	float y_coord;
	float z_coord;

	if (ind.row < pind->current_cell->r) { //north
		y_coord = pind->current_cell->cell_corners.ulc[1]; //maximum y coordinate
		x_coord = x_dist(eng); //random x coordinate
		z_coord = z_dist(eng); //random z coordinate
	}
	else if( ind.row> pind->current_cell->r){ //south
		y_coord = pind->current_cell->cell_corners.llc[1]; //minimum y coordinate
		x_coord = x_dist(eng);
		z_coord = z_dist(eng);
	}
	else { //same row
		if (ind.col < pind->current_cell->c) {//west
			y_coord = y_dist(eng);
			x_coord = pind->current_cell->cell_corners.ulc[0]; //minimum x coordinate
			z_coord = z_dist(eng); //random z coordinaet
		}
		else if (ind.col > pind->current_cell->c) { //east
			y_coord = y_dist(eng);
			x_coord = pind->current_cell->cell_corners.urc[0]; //maximum x coordinate
			z_coord = z_dist(eng); //random z coordinaet
		}
		else if (ind.layer < pind->current_cell->l) { //up
			y_coord = y_dist(eng);
			x_coord = x_dist(eng);
			z_coord = pind->current_cell->min_depth; //start on top of the cell
		}
	}


	pind->current_pos.x = x_coord; pind->current_pos.y = y_coord; pind->current_pos.z = z_coord; 
	//pind->memory.push(pind->current_pos); //create the first memory entry

	pind->m_x_movements.push_back(x_coord);
	pind->m_y_movements.push_back(y_coord);
	pind->m_z_movements.push_back(z_coord); //add to coordinate record
}

int Subpopulation::init_number() {//initialise only enough individuals in the patch/cells that are in the right depth range!

	int init;

	if (mod->patchmodel == 1) {

		//int mindepth = pop_dyn->Depthmin; int maxdepth = pop_dyn->Depthmax; //get the population's depth range
		//cout << "patch has " << patch_num->included_cells.size() << " included cells" << endl;
		//for (int p = 0; p < patch_num->included_cells.size(); p++) {
		//	if (mindepth <= patch_num->included_cells[p]->min_depth && patch_num->included_cells[p]->max_depth <= maxdepth) {
		//		patch_num->indepth_cells.push_back(patch_num->included_cells[p]);
		//	}
		//}
		//patch_num->calc_K("indepth_info"); //recalculate maximum individuals for how many cells are suitable
		//cout << "patch's K is " << patch_num->patch_info.K;
		if (mod->initdens == 0) { //at K
			init = patch_num->indepth_info.max_inds;
			//cout << "at carrying capacity: " << init << " indivs " << endl;
		}
		else if (mod->initdens == 1) { //at K/2
			init = floor((patch_num->indepth_info.max_inds) / 2);
		}
		else if (mod->initdens == 2) { //at specified density
			init = floor(mod->indsha* patch_num->indepth_info.size_ha);
		}
		
	}
	else { //cell-based
		if (mod->initdens == 0) { //at K
			init = cell_num->max_inds; 
		}
		else if (mod->initdens == 1) { //at K/2
			init = floor((cell_num->max_inds)) / 2;
		}
		else if (mod->initdens == 2) { //at specified density
			init = floor(mod->indsha * cell_num->ha); 
		}
	}

	return init;
}

void Subpopulation::init_indivs() {
	//poisson_distribution<int>poisson(init_number());
	//int init=poisson(eng); //calculate how many idividuals to initialise
	int init = init_number(); //calculate how many idividuals to initialise
	//cout << "should initialise " << init << " individuals in subpop " << cell_num->cell_ID << "which has " << cell_num->ha << " hectares and " << cell_num->max_inds << endl;
	if (mod->stagestruct == 1) { //if it's stagestructured
		for (int s = 1; s < mod->stages; s++) { //for each stage (starts at 1 because you dont initialise juveniles)
			float prop = mod->stage_props[s]; //what proportion of the initialised indivs are in this stage
			int n = floor(prop*init); //how many indivs is that?
			for (int i = 0; i < n; i++) { //for each of those individuals
				//Individual * pind = new Individual(mod, pop_dyn, lands, pop_einfo, pop_tinfo, pop_sinfo, get_stage_info(s), assign_sex(), pop_dyn->total_births); //create a new stage structured individual
				Individual * pind;
				if (mod->repro > 0) {
					pind = new Individual(mod, pop_dyn, lands, get_stage_info(s), assign_sex(), pop_dyn->total_births, true, &pop_stageinfo); //create a new stage structured individual
				}
				else { //if it's an asexual model, i don't really need the genetic structure of the parent 
					pind = new Individual(mod, pop_dyn, lands, get_stage_info(s),1, pop_dyn->total_births, true, &pop_stageinfo); //create a new stage structured individual
				}
				//cout << "initialised, ";
				pop_dyn->total_births++;
				if (mod->patchmodel == 1) { pind->natal_patch = patch_num->patch_ID;  }
				else { pind->natal_cell = cell_num->cell_ID; }

				if (mod->repro > 0) { //if it's a sexual model
					if (pind->sex == 0 ){
						nmales++;//increment total number of males
						if (pind->pstage->m_fec > 0) { //if this male is reproductively mature, add to rep_males
							rep_males.push_back(pind);
						}
						
					} 
					else if (pind->sex==1 ) { //check female fecundity
						nfemales++; //increment the total number of females in the subpop if the female is reproductively mature
						if (pind->pstage->f_fec > 0) {
							rep_females.push_back(pind); //add to rep_females
						}
					}
				}
				//none of this is required for an asexual model because there are no possible mating requirements 
				//cout << "ready to ";
				init_age(pind, s);
				//cout << "got past age, ";
				init_size(pind);
				//cout << " and size, ";
				init_position(pind);
				//cout << "and position, ";
				indivs.push_back(pind);
				size++; //increase subpopulation size by 1
			}
		}
	}
	else { //non-overlapping generations
		for (int i = 0; i < init; i++) { //for each of those individuals
			//Individual * pind = new Individual(mod, pop_dyn, lands, pop_einfo, pop_tinfo, pop_sinfo, assign_sex(), pop_dyn->total_births); //create a new individual
			Individual * pind;
			if (mod->repro > 0) {
				pind = new Individual(mod, pop_dyn, lands, get_stage_info(0), assign_sex(), pop_dyn->total_births, true, &pop_stageinfo); //create a new individual
			}
			else {
				pind = new Individual(mod, pop_dyn, lands, get_stage_info(0), 1, pop_dyn->total_births, true, &pop_stageinfo); //create a new individual
			}
			pop_dyn->total_births++;
			if (mod->patchmodel == 1) { pind->natal_patch = patch_num->patch_ID; }
			else { pind->natal_cell = cell_num->cell_ID; }
			//current_cell is set in init_position()

			if (mod->repro > 0) { //if it's a sexual model
				if (pind->sex == 0) {
					nmales++;//increment total number of males
					rep_males.push_back(pind);
				}
				else {
					nfemales++;
					rep_females.push_back(pind);
				}
			}
			//none of this is needed for an asexual model with non overlapping generations because no mating requirements are possible

			init_position(pind);
			init_size(pind);
			indivs.push_back(pind);
			size++; //increase subpopulation size by 1
		}
	}
	//cout << "initialised " << size << " indivs in subpop " << cell_num->cell_ID << endl;
}



//reproduction

void Subpopulation::reproduction() {
	if (mod->repro == 0) { //asexual
		for (int i = 0; i < indivs.size(); i++) {
			int noff;
			noff = indivs[i]->reproduction(); //how many offspring will it produce?
			if (noff != 0 && noff != -9) { //if it produced any
				for (int n = 0; n < noff; n++) {
					Individual * pind; //create a new individual for each offspring
					if (mod->stagestruct == 1) { //if it's stagestructured
						pind = new Individual(mod, pop_dyn, lands, get_stage_info(0), 1, pop_dyn->total_births, false, &pop_stageinfo); //sex=1 because it's an asexual model
						//pind = new Individual(mod, pop_dyn, lands, pop_einfo, pop_tinfo, pop_sinfo, get_stage_info(0), 1, pop_dyn->total_births); //sex=1 because it's an asexual model
						//init_size(pind);
						pop_dyn->total_births++;
					}
					else { //non overlapping generations
						pind = new Individual(mod, pop_dyn, lands, get_stage_info(0), 1, pop_dyn->total_births, false, &pop_stageinfo); //sex=1 because it's an asexual model
						//pind = new Individual(mod, pop_dyn, lands, pop_einfo, pop_tinfo, pop_sinfo, 1, pop_dyn->total_births);
						pop_dyn->total_births++;
					}
					//inheritance, this may override some of the previous step if evolution is allowed
					pind->inheritance(indivs[i], nullptr); //get all of parent's information (also the pop_dyn and pop_e/t/sinfo pointers), as well as position and patch/cell pointers
					if (mod->stagestruct == 0) { /*size++;*/ juv_indivs.push_back(pind); } //if it's non overlapping gens, offspring go to juvs_vector 
					else { /*size++;*/ indivs.push_back(pind); } //otherwise join whole subpopulation

				}
			}
		}
	}

	if (mod->repro > 0) { if (rep_males.size() == 0 || rep_females.size() == 0) { return; } } //if it's a sexual model and there are either no males or no females, reproduction can't happen

	if (mod->repro == 1) {//broadcast spawning/promiscuity: offspring are sired by multiple males
		/*default_random_engine m_generator;*/
		uniform_int_distribution<int> sex_dist(0.0, rep_males.size() - 1);
		for (int f = 0; f < rep_females.size(); f++) {
			int noff;

			noff = rep_females[f]->reproduction(); //how many offspring will it produce?

			if (noff != 0 && noff != -9) { //if it produced any
				for (int n = 0; n < noff; n++) {
					Individual * pind; //create a new individual for each offspring
					if (mod->stagestruct == 1) { //if it's stagestructured
						pind = new Individual(mod, pop_dyn, lands, get_stage_info(0), assign_sex(), pop_dyn->total_births, false, &pop_stageinfo);
						//pind = new Individual(mod, pop_dyn, lands, pop_einfo, pop_tinfo, pop_sinfo, get_stage_info(0), assign_sex(), pop_dyn->total_births);
						pop_dyn->total_births++;
					}
					else { //non overlapping generations
						pind = new Individual(mod, pop_dyn, lands, get_stage_info(0), assign_sex(), pop_dyn->total_births, false, &pop_stageinfo);
						//pind = new Individual(mod, pop_dyn, lands, pop_einfo, pop_tinfo, pop_sinfo, assign_sex(), pop_dyn->total_births);
						pop_dyn->total_births++;
					}

					int mate_index = sex_dist(eng);
					pind->inheritance(rep_females[f], rep_males[mate_index]);//get all of parent's information (also the pop_dyn and pop_einfo pointers), as well as position and patch/cell pointers
					if (mod->stagestruct == 0) { /*size++; */juv_indivs.push_back(pind); } //if it's non overlapping gens, offspring go to juvs_vector 
					else { /*size++;*/ indivs.push_back(pind); } //otherwise join whole subpopulation
				}
			}
		}
	}
	
	else if (mod->repro == 2) {//harem or monogamy, offspring are sired by only one male
		//temporary vectors to record which individuals have mated
		vector<Individual*> ran_males;
		copy(rep_males.begin(), rep_males.end(), ran_males.begin());
		vector<Individual*> ran_females;
		copy(rep_females.begin(), rep_females.end(), ran_females.begin());
		random_device rd;
		mt19937 g(rd());

		shuffle(ran_males.begin(), ran_males.end(), g);
		shuffle(ran_females.begin(), ran_females.end(), g);

		int f = 0; //keep track of which female in the vector
		for (int m = 0; m < ran_males.size(); m++) {
			int mates = 0; //how many mates has this male had
			while (mates != pop_dyn->Hsize &&  f < ran_females.size()) { //each male only gets Hsize females to mate with, and if there are no females left, stop
				int noff = ran_females[f]->reproduction(); //check if she reproduces/how many offspring she produces
				if (noff != -9) { //if the female reproduces
					for (int n = 0; n < noff; n++) {
						Individual * pind; //create a new individual for each offspring
						if (mod->stagestruct == 1) { //if it's stagestructured
							pind = new Individual(mod, pop_dyn, lands, get_stage_info(0), assign_sex(), pop_dyn->total_births, false, &pop_stageinfo); 
							//pind = new Individual(mod, pop_dyn, lands, pop_einfo, pop_tinfo, pop_sinfo, get_stage_info(0), assign_sex(), pop_dyn->total_births); 
							pop_dyn->total_births++;
						}
						else { //non overlapping generations
							pind = new Individual(mod, pop_dyn, lands, get_stage_info(0), assign_sex(), pop_dyn->total_births, false,  &pop_stageinfo);
							pop_dyn->total_births++;
						}
						
						pind->inheritance(ran_females[f], ran_males[m]);//get all of parent's information (also the pop_dyn and pop_einfo pointers), as well as position and patch/cell pointers
						if (mod->stagestruct == 0) { /*size++;*/ juv_indivs.push_back(pind); } //if it's non overlapping gens, offspring go to juvs_vector 
						else { /*size++;*/ indivs.push_back(pind); } //otherwise join whole subpopulation
					}
					mates++; //increment the number of mates this male has
				}
				f++; //increment the female index
			}
		}
	}
	
	if (mod->stagestruct == 0) { //if it's non-overlapping generations
		nmales -= rep_males.size();
		rep_males.clear(); //clear the reproductively mature males
		nfemales -= rep_females.size();
		rep_females.clear(); //and females
	} //this is so arriving individuals into a patch don't think they've found a mate when the adult generation just hasnt died yet
	//for stage structured pops, these vectors need to remain in tact for settlement
}

//emigration

void Subpopulation::juv_release(Individual* pind) { //this function already assumes that the parent of this individual is settled on a cell side accessible by water!
	if (pind->current_cell == nullptr) {
		cout << "before doing juv release, current cell was nullptr" << endl;
	}
	//cout << "releasing at " << pind->current_pos.x << ", " << pind->current_pos.y << ", " << pind->current_pos.z << endl;
	int drow=0, dcol=0, dlayer=0;
	int prev_r=pind->current_cell->r, prev_c=pind->current_cell->c, prev_l=pind->current_cell->l;
	if (pind->current_pos.z == pind->current_cell->min_depth) { //if it's settled on top of the cell
		dlayer = -1; //nearest neighbour for release is above it
		//cout << "up";
	}
	else if (pind->current_pos.x == pind->current_cell->cell_corners.llc[0]) { //if it's settled on the minimum x face of the cell
		dcol = -1; //move one column left
		//cout << "left ";
	}
	else if (pind->current_pos.x == pind->current_cell->cell_corners.lrc[0]) { //if it's settled on the maximum x face of the cell
		dcol = 1; //move one column right
		//cout << "right ";
	}
	else if (pind->current_pos.y == pind->current_cell->cell_corners.llc[1]) { //if it's settled on the minimum y face of the cell
		drow = 1; //move one row down
		//cout << "south ";
	}
	else if (pind->current_pos.y == pind->current_cell->cell_corners.ulc[1]) { //if it's settled on the maximum y face of the cell
		drow = -1; //move one row up
		//cout << "north";
	}
	else {
		cout << "in juv_release, coordinate conditions not met, i've missed something" << endl;
	}
	//cout << endl;

	//update current cell
	float new_x=pind->current_pos.x, new_y = pind->current_pos.y, new_z = pind->current_pos.z;
	//cout << "in juv release, the habtype of new cell is :" << lands->raster_3d[prev_r + drow][prev_c + dcol][prev_l + dlayer]->habtype << endl;
	//move 2m into next cell
	if (drow > 0) { new_y -= 2; }
	else if (drow < 0) { new_y += 2; }
	if (dcol > 0) { new_x += 2; }
	else if (dcol < 0) { new_x -= 2; }
	if (dlayer < 0) { new_z -= 2; }
	//for depth changes, it depends on buoyancy limits and diel vertical migration
	//if (pop_tinfo->diel_vert == 1) { //if there is diel vertical migration, assume it is release at night, so upper dv range
	//if (pop_stageinfo[0]->s_trans.diel_vert == 1) { //if there is diel vertical migration, assume it is release at night, so upper dv range
	//	int min_dv = pop_stageinfo[0]->s_trans.buoy_min;
	//	//int min_dv = pop_tinfo->buoy_min;
	//	int max_dv = min_dv + pop_stageinfo[0]->s_trans.dv_range;
	//	//int max_dv = min_dv + pop_tinfo->dv_range;
	//	uniform_int_distribution<int> dv_dist(min_dv, max_dv);
	//	if (new_z < min_dv || new_z > max_dv) { //if current depth is outside those boundaries, 
	//		new_z = dv_dist(eng); //give it a random location in the dv range
	//	}
	//	//otherwise don't change the depth
	//	//cout << "due to dv migration, new_z is " << new_z << endl;
	//}
	//for depth changes, depends on juv_release input parameter
	if (pop_dyn->juv_rel != -9) {
		new_z -= pop_dyn->juv_rel;
	}
	if(pop_stageinfo[0]->s_trans.buoy_min != -9) { //if buoyancy limits are given but there is NO DIEL VERTICAL MIGRATION, else if(pop_tinfo->buoy_min != -9) {
		if (new_z < pop_stageinfo[0]->s_trans.buoy_min) { //if position is above minimum, new_z < pop_tinfo->buoy_min
			new_z = pop_stageinfo[0]->s_trans.buoy_min + 1; //put it within the buoyancy range, pop_tinfo->buoy_min + 1;
		}
		if (new_z > pop_stageinfo[0]->s_trans.buoy_max) { //if position is below the maximum , pop_tinfo->buoy_max
			new_z = pop_stageinfo[0]->s_trans.buoy_max - 1; //put it within the buoyancy range
		}
		//cout << "due to buoy limits, new_z is" << new_z << endl;
	}

	
	//new_z = (pind->current_cell->max_depth- pind->current_cell->min_depth)/2; //TRIAL: starting them in the midpoint depth so that downward vector movements work better
	vector<int>indexes = pind->calc_r_c_l(new_x, new_y, new_z);

	
	if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]] == nullptr) { 
		cout << "in juv release trying to assign cell " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << endl;

		cout << "assigning cell didn't work in juv_release" << endl; 
	}
	else if (lands->raster_3d[indexes[0]][indexes[1]][indexes[2]]->habtype != 0) {
		cout << "juv release should release into water but hits habitat at " << indexes[0] << ", " << indexes[1] << ", " << indexes[2] << ", at depth " << new_z << endl;
		cout << "coming from cell " << pind->current_cell->r << ", " << pind->current_cell->c << ", " << pind->current_cell->l << 
			"at coordinates " << pind->current_pos.x << ", " << pind->current_pos.y << ", " << pind->current_pos.z <<", drow=" << drow << ", dcol="<< dcol << ", dlay=" << dlayer<< endl;
	}
	pind->current_pos.x = new_x; pind->current_pos.y = new_y; pind->current_pos.z = new_z;
	pind->m_x_movements.push_back(new_x); pind->m_y_movements.push_back(new_y); pind->m_z_movements.push_back(new_z);
	pind->current_cell = lands->raster_3d[indexes[0]][indexes[1]][indexes[2]]; //need to update its current cell
	//cout << "after emig, z is " << pind->current_pos.z  << "and habtype is " << pind->current_cell->habtype << "with w value " << pind->current_cell->w << " and v value " << pind->current_cell->v << endl;
	//pind->memory.pop();
	//pind->memory.push(pind->current_pos); //add it to memory
}

void Subpopulation::emigration() {
	if (mod->stagestruct == 0) { //if non overlapping generations, only go through new offspring for emigration probability
		for (int j = 0; j < juv_indivs.size(); j++) {
			bool e = juv_indivs[j]->check_emig();
			if (e==true) {//if it emigrates
				juv_indivs[j]->emig = true;
				juv_indivs[j]->status = 1;
				if (mod->size_or_time == 0) { juv_indivs[j]->calc_SL_rho(); /*cout << "indiv " << juv_indivs[j]->ID << " has SL " << juv_indivs[j]->dinfo.SL << endl;*/ } //if behaviour is size-dependent, they need to calculate SL and rho
				juv_release(juv_indivs[j]);
				leave(juv_indivs[j]); 

			}
			else { //if it stays 
				juv_indivs[j]->status = 2;
				//implementing density dependence in self-recruitment
				float s = 1 / (1 + exp(-1 * (size / patch_num->patch_info.K - 1)*-6)); //numerator is supposed to be S0, but using S0=1, alpha=-6, beta=1
				float ran_s = unif_dist(eng);
				
				if (ran_s < s) { //let it settle
					//cout << "self recruiter has settled" << endl;
					size++;
					surv_indivs.push_back(juv_indivs[j]); //add to surviving indivs vector
					if (juv_indivs[j]->current_cell->habtype == 0) {
						cout << "local recruit but has current cell water type" << endl;
					}
				}
				else { //otherwise, it's dead
					//size--;
					juv_indivs[j]->dead = true;
					delete juv_indivs[j]; //it will be cleared anyway so it can be deleted
					//cout << "due to density dependence, self recruiter has died" << endl;
				}	
			}
		} //adults will die anyway because there are no overlapping generations so no need to go through indivs vector
		juv_indivs.clear();
	} 
	else { //if it is stage structured, go through everyone. juvs_vector is empty
		//cout << "there are " << indivs.size() << " indivs in subpop " << patch_num->patch_ID << " which matches size " << size << endl;
		for (int i = 0; i < indivs.size(); i++) {
			bool e;
			if (pop_stageinfo[0]->s_emig.stagedep==0 && indivs[i]->pstage->stage != pop_stageinfo[0]->s_emig.emigstage) { //if only one stage emigrates, and this indiv isnt in that stage
				surv_indivs.push_back(indivs[i]);
				continue; //skip it

			}
			else if (pop_stageinfo[0]->s_emig.stagedep == 1 ) { //more than one stage can disperse
				if (pop_stageinfo[0]->s_emig.sexdep == 1 && indivs[i]->sex==0) { //sex dependent and it's male
					if (indivs[i]->pstage->s_emig.m_emig_prob == -9) { //if this stage/sex doesn't emigrate
						surv_indivs.push_back(indivs[i]);
						continue;
					}
				}
				else { //sex independent or female
					if (indivs[i]->pstage->s_emig.f_emig_prob == -9) { //if this stage/sex doesn't emigrate
						surv_indivs.push_back(indivs[i]);
						continue;
					}
				}
			}
			//if it hasn't continued by now, then it emigrates
			e = indivs[i]->check_emig();
			
			if (e==true) {//if it emigrates
				indivs[i]->emig = true;
				indivs[i]->status = 1;
				if (mod->size_or_time == 0) { indivs[i]->calc_SL_rho();} //if behaviour is size-dependent, they need to calculate SL and rho
				juv_release(indivs[i]);
				leave(indivs[i]);
				
			}
			else { //if it stays 
				indivs[i]->status = 2;
				bool dies = false;
				int max_inds = 0;
				if (mod->patchmodel == 0) { //cell based
					max_inds = cell_num->max_inds;
					if (size > cell_num->max_inds) {
						dies = true;
					}
				}
				else {
					max_inds = patch_num->patch_info.max_inds;
					if (size > patch_num->patch_info.max_inds) {
						//size--;
						dies = true;
					}
				}
				if (dies == true) {
					indivs[i]->dead = true;
					delete(indivs[i]);
				}
				else {
					//cout << "size is " << size << ", K is " << patch_num->patch_info.K << ", max inds is " << patch_num->patch_info.max_inds << ", ";
					//otherwise do density dependence
					float s = 1 / (1 + exp(-1 * (size / max_inds - 1)*-6)); //numerator is supposed to be S0, but using S0=1, alpha=-6, beta=1
					float ran_s = unif_dist(eng);
					//cout << "size is " << size << ", K is " << patch_num->patch_info.K << ", max inds is " << patch_num->patch_info.max_inds<< ", ";
					//cout << "s is " << s << " and ran_s is " << ran_s;
					if (ran_s < s) {
						//cout << "self recruiter has settled, " ;
						size++;
						surv_indivs.push_back(indivs[i]); //add to surviving indivs vector
						//cout << "size is " << size << endl;
					}
					else {
						//cout << "s is " << s << " and ran_s is " << ran_s;
						//size--;
						indivs[i]->dead = true;
						delete(indivs[i]);
						//cout << "due to density dependence, self recruiter has died and size is now " << size  << endl;
					}
				}
				
				
			}
			
		}
	}
	//cout << "size of surv_indivs is " << surv_indivs.size() << " in subpop " << patch_num->patch_ID << endl;
	//cout << "size of indivs is " << indivs.size() << " matching size " << size << endl;
	//cout << "there are " << pop_matrix->matrix_indivs.size() << " individuals in the matrix, matches total " << pop_matrix->total << endl;
}


//survival
void Subpopulation::survival() {
	//cout << "at annual survival, there are " << surv_indivs.size() << " indivs in subpop, matching size " << size << endl;
	if (mod->stagestruct == 0) { //non overlapping generations
		for (int i = 0; i < indivs.size(); i++) { 
			indivs[i]->dead = true;//previous adult generation all dies
			size--;
		}
	} //there is no explicit survival phase for newly settled individuals

	else { //stage structured pops
		for (int i = 0; i < surv_indivs.size(); i++) { //all individuals are now in surv_indivs

			if (surv_indivs[i]->survival()==false) { //if individual dies
				surv_indivs[i]->dead = true;
				size--; //update size of subpopulation
			}
			else {
				if (pop_dyn->fishing_mort != -9) { //if it didn't die before and fishing mortality is present
					if (surv_indivs[i]->fishing_mort() == true) { //if it dies by fishing
						surv_indivs[i]->dead = true; 
						size--;
					}
					if (surv_indivs[i]->status < 0) { cout << "after fishing mort in subpop surv, status is weird " << endl; }
				}
			}
			
			if (surv_indivs[i]->status < 0) { cout << "at subpop survival, status is weird " << endl; }
			
		}

	}
}

void Subpopulation::development(bool aging) {
	//cout << "starting development" << endl;
	//now that theyve all settled, i can update the rep_males/females vectors with only surviving individuals
	nmales = 0;
	rep_males.clear(); //clear the reproductively mature males
	nfemales = 0;
	rep_females.clear(); //and females

	for (int i = 0; i < surv_indivs.size(); i++) {
		if (surv_indivs[i]->status < 0) { cout << "in development, status is weird" << endl; }
		if (surv_indivs[i]->dead == true) { continue; } //if indiv is dead, skip it
		if (aging) { surv_indivs[i]->age++; } //if aging=true that means this development is happening at the end of the year, so we want to age +1
		if (surv_indivs[i]->pstage->stage == mod->stages) {
			if (mod->repro > 0) {
				if (surv_indivs[i]->sex == 0 ) {
					nmales++;
					if (surv_indivs[i]->pstage->m_fec > 0) {
						rep_males.push_back(surv_indivs[i]);
					}
				}
				else {
					nfemales++;
					if (surv_indivs[i]->pstage->f_fec > 0) {
						rep_females.push_back(surv_indivs[i]);
					}
				}
			}
			
			continue;
		}
		else {
			surv_indivs[i]->stage_development(); //see if they develop to next stage
			if (mod->repro > 0) {
				if (surv_indivs[i]->sex == 1) {
					nfemales++;
					if (surv_indivs[i]->pstage->f_fec > 0) { rep_females.push_back(surv_indivs[i]); } //are they reproductively mature?
				}
				else {
					nmales++;
					if (surv_indivs[i]->pstage->m_fec > 0) { rep_males.push_back(surv_indivs[i]); }
				}
			}
		}
	}
}


void Subpopulation::clean_up() {
	if (mod->stagestruct == 0) { //non overlapping generations
		for (int i = 0; i < indivs.size(); i++) { //adult generation all dies, so clear indivs vector (dispersers that entered were added to surv_indivs)
			delete indivs[i]; //only for non overlapping generations can indivs be deleted like this, because the pointers are not used anywhere else
		}
		//rep_males and rep_females only contain juveniles that have just arrived (because i cleared rep males/females at end of repro for nonstructured pops), so leave them in tact
		size = surv_indivs.size(); //update size of population
		indivs.clear();
		indivs = surv_indivs; //there is no explicit survival phase for non-overlapping generations
		surv_indivs.clear();
	}
	else { //stage structured pop
		indivs.clear(); //this will delete empty pointers, but for stage structured this doesnt deallocate memory!

		nmales = 0;
		rep_males.clear(); //clear the reproductively mature males
		nfemales = 0;
		rep_females.clear(); //and females

		for (int i = 0; i < surv_indivs.size(); i++) { //all individuals are now in surv_indivs
			if (surv_indivs[i]->status < 0) { cout << "in cleanup, status is weird " << endl; }
			if (surv_indivs[i]->dead==false) { //if individual hasn't died
				indivs.push_back(surv_indivs[i]); //add to individuals vector
				
				if (mod->repro > 0) {
					if (surv_indivs[i]->sex == 0) {
						nmales++;
						if (surv_indivs[i]->pstage->m_fec > 0) {
							rep_males.push_back(surv_indivs[i]);
						}
					}
					else {
						nfemales++;
						if (surv_indivs[i]->pstage->f_fec > 0) {
							rep_females.push_back(surv_indivs[i]);
						}
					}
				}

			}
			else {
				delete surv_indivs[i]; //otherwise deallocate the memory for that individual
			}
		}
		surv_indivs.clear(); //clear the surv_indivs vector ready for the next reproductive season

	}
	if (lands->patch_extinct == 1) {
		if (lands->p_ext_patch == patch_num->patch_ID) {
			cout << "size of surv_indivs is " << surv_indivs.size() << " and indivs is " << indivs.size() << endl;
		}
	}
}

//saving outputs
vector<float> Subpopulation::summarise() {
	vector<float> values;
	if (mod->patchmodel == 0) { //if cell-based
		values.push_back(cell_num->cell_ID); values.push_back(cell_num->r); values.push_back(cell_num->c); values.push_back(cell_num->l); values.push_back(-9);
	}
	else { //patch-based
		values.push_back(-9); values.push_back(patch_num->patch_ID); 
	}
	values.push_back(size);
	if (mod->repro == 1) { //if sexual reproduction
		values.push_back(nfemales); values.push_back(nmales);
	}
	else { //asexual, no need to fill these
		values.push_back(-9); values.push_back(-9);
	}

	return(values);

}

vector<float> Subpopulation::mean_evolution(string param, bool first) { //param can be emigration, transfer, growth, or settlement
	//cout << "called mean evolution function, ";
	vector<float> values;
	vector<Individual*>* vec = &indivs;
	//if (first == true) {
	//	//cout << "in year0, there are " << indivs.size() << " surviving indivs in subpop " << patch_num->patch_ID;
	//	vec = &indivs;
	//	//cout << "this matches vec size " << vec->size() << endl;
	//}
	//else {
	//	vec = &surv_indivs;
	//}

	if (param == "emigration") {
		float ep_sum = 0;
		float alpha_sum = 0;
		float beta_sum = 0;
		for (int ind = 0; ind < vec->size(); ind++) {
			if (pop_stageinfo[0]->s_emig.densdep == 1) { //density dependent
				if (mod->repro == 0) {//asexual
					ep_sum += ((*vec)[ind]->dinfo.parent->D0[0][0]);
					alpha_sum += ((*vec)[ind]->dinfo.parent->alpha[0][0]);
					beta_sum += ((*vec)[ind]->dinfo.parent->beta[0][0]);
				}
				else {
					if (pop_stageinfo[0]->s_emig.sexdep == 1 && (*vec)[ind]->sex == 0) { //sex dependent
						ep_sum += ((*vec)[ind]->dinfo.parent->D0[1][0] + (*vec)[ind]->dinfo.parent->D0[1][1]) / 2;
						alpha_sum += ((*vec)[ind]->dinfo.parent->alpha[1][0] + (*vec)[ind]->dinfo.parent->alpha[1][1]) / 2;
						beta_sum += ((*vec)[ind]->dinfo.parent->beta[1][0] + (*vec)[ind]->dinfo.parent->beta[1][1]) / 2;
					}
					else { //female
						ep_sum += ((*vec)[ind]->dinfo.parent->D0[0][0] + (*vec)[ind]->dinfo.parent->D0[0][1]) / 2;
						alpha_sum += ((*vec)[ind]->dinfo.parent->alpha[0][0] + (*vec)[ind]->dinfo.parent->alpha[0][1]) / 2;
						beta_sum += ((*vec)[ind]->dinfo.parent->beta[0][0] + (*vec)[ind]->dinfo.parent->beta[0][1]) / 2;
					}
				}
				
				//ep_sum += (*vec)[ind]->dinfo.D0;
				//alpha_sum += (*vec)[ind]->dinfo.alpha;
				//beta_sum += (*vec)[ind]->dinfo.beta;
			}
			else {
				if (mod->repro == 0) {//asexual
					ep_sum += ((*vec)[ind]->dinfo.parent->emig_prob[0][0]);
				}
				else {
					if (pop_stageinfo[0]->s_emig.sexdep == 1 && (*vec)[ind]->sex == 0) { //sex dependent
						ep_sum += ((*vec)[ind]->dinfo.parent->emig_prob[1][0] + (*vec)[ind]->dinfo.parent->emig_prob[1][1]) / 2;
					}
					else {
						ep_sum += ((*vec)[ind]->dinfo.parent->emig_prob[0][0] + (*vec)[ind]->dinfo.parent->emig_prob[0][1]) / 2;
					}
				}
				
			}
		}
		values.push_back(ep_sum / vec->size());
		if (pop_stageinfo[0]->s_emig.densdep == 1) {
			values.push_back(alpha_sum / vec->size());
			values.push_back(beta_sum / vec->size());
		}
		return values;
	}
	else if (param == "transfer") {
		float active_sum = 0, dv_sum = 0;
		for (int ind = 0; ind < vec->size(); ind++) {
			if (mod->repro == 0) {//asexual
				if (mod->size_or_time == 1) { //time depedent
					active_sum += ((*vec)[ind]->dinfo.parent->min_active_time[0][0]);
					dv_sum += ((*vec)[ind]->dinfo.parent->min_dv_time[0][0]);
				}
				else {
					active_sum += ((*vec)[ind]->dinfo.parent->min_active_size[0][0]);
					dv_sum += ((*vec)[ind]->dinfo.parent->min_dv_size[0][0]);
				}
			}
			else {
				if (pop_stageinfo[0]->s_trans.sexdep == 1 && (*vec)[ind]->sex == 0) { //sex dependent
					if (mod->size_or_time == 1) { //time depedent
						active_sum += ((*vec)[ind]->dinfo.parent->min_active_time[1][0] + (*vec)[ind]->dinfo.parent->min_active_time[1][1]) / 2;
						dv_sum += ((*vec)[ind]->dinfo.parent->min_dv_time[1][0] + (*vec)[ind]->dinfo.parent->min_dv_time[1][1]) / 2;
					}
					else {
						active_sum += ((*vec)[ind]->dinfo.parent->min_active_size[1][0] + (*vec)[ind]->dinfo.parent->min_active_size[1][1]) / 2;
						dv_sum += ((*vec)[ind]->dinfo.parent->min_dv_size[1][0] + (*vec)[ind]->dinfo.parent->min_dv_size[1][1]) / 2;
					}
				}
				else {
					if (mod->size_or_time == 1) { //time depedent
						active_sum += ((*vec)[ind]->dinfo.parent->min_active_time[0][0] + (*vec)[ind]->dinfo.parent->min_active_time[1][1]) / 2;
						dv_sum += ((*vec)[ind]->dinfo.parent->min_dv_time[0][0] + (*vec)[ind]->dinfo.parent->min_dv_time[1][1]) / 2;
					}
					else {
						active_sum += ((*vec)[ind]->dinfo.parent->min_active_size[0][0] + (*vec)[ind]->dinfo.parent->min_active_size[0][1]) / 2;
						dv_sum += ((*vec)[ind]->dinfo.parent->min_dv_size[0][0] + (*vec)[ind]->dinfo.parent->min_dv_size[0][1]) / 2;
					}
				}
			}
			
			
		}
		values.push_back(active_sum / vec->size());
		values.push_back(dv_sum / vec->size());
		return values;
	}
	else if (param == "growth") {
		if (pop_stageinfo[0]->s_trans.pop_grow.method == 0) {
			float m_sum = 0;
			float b_sum = 0;
			
			for (int ind = 0; ind < vec->size(); ind++) {
				if (mod->repro == 0) {
					m_sum += ((*vec)[ind]->dinfo.parent->m[0][0]);
				}
				else {
					if (pop_stageinfo[0]->s_trans.pop_grow.sex_dep == 1 && (*vec)[ind]->sex == 0) { //sex dependent
						m_sum += ((*vec)[ind]->dinfo.parent->m[1][0] + (*vec)[ind]->dinfo.parent->m[1][1]) / 2;
					}
					else {
						m_sum += ((*vec)[ind]->dinfo.parent->m[0][0] + (*vec)[ind]->dinfo.parent->m[0][1]) / 2;
					}
				}
				
			}
			values.push_back(m_sum / vec->size());
			values.push_back(b_sum / vec->size());
		}
		else if (pop_stageinfo[0]->s_trans.pop_grow.method == 1) {
			float linf_sum = 0;
			float k_sum = 0;
			float ti_sum = 0;
			for (int ind = 0; ind < vec->size(); ind++) {
				if (mod->repro == 0) {
					linf_sum += ((*vec)[ind]->dinfo.parent->Linf[0][0]);
					k_sum += ((*vec)[ind]->dinfo.parent->G_K[0][0]);
					ti_sum += ((*vec)[ind]->dinfo.parent->Ti[0][0]);
				}
				else {
					if (pop_stageinfo[0]->s_trans.pop_grow.sex_dep == 1 && (*vec)[ind]->sex == 0) { //sex dependent
						linf_sum += ((*vec)[ind]->dinfo.parent->Linf[1][0] + (*vec)[ind]->dinfo.parent->Linf[1][1]) / 2;
						k_sum += ((*vec)[ind]->dinfo.parent->G_K[1][0] + (*vec)[ind]->dinfo.parent->G_K[1][1]) / 2;
						ti_sum += ((*vec)[ind]->dinfo.parent->Ti[1][0] + (*vec)[ind]->dinfo.parent->Ti[1][1]) / 2;
					}
					else {
						linf_sum += ((*vec)[ind]->dinfo.parent->Linf[0][0] + (*vec)[ind]->dinfo.parent->Linf[0][1]) / 2;
						k_sum += ((*vec)[ind]->dinfo.parent->G_K[0][0] + (*vec)[ind]->dinfo.parent->G_K[0][1]) / 2;
						ti_sum += ((*vec)[ind]->dinfo.parent->Ti[0][0] + (*vec)[ind]->dinfo.parent->Ti[0][1]) / 2;
					}
				}
				
				
			}
			values.push_back(linf_sum / vec->size());
			values.push_back(k_sum / vec->size());
			values.push_back(ti_sum / vec->size());
		}
		else {

		}
		return values;
	}
	//else if (param == "settlement") {
	else{
		float comp_sum = 0, sett_sum = 0, alphaS_sum = 0, betaS_sum = 0;
		for (int ind = 0; ind < vec->size(); ind++) {
			if (mod->repro == 0) {
				if (mod->size_or_time == 1) { //time depedent
					comp_sum += ((*vec)[ind]->dinfo.parent->comp_time[0][0]);
				}
				else {
					comp_sum += ((*vec)[ind]->dinfo.parent->comp_size[0][0]);
				}
				sett_sum += ((*vec)[ind]->dinfo.parent->S0[0][0]);
			}
			else {
				if (pop_stageinfo[0]->s_sett.sexdep == 1 && (*vec)[ind]->sex == 0) { //sex dependent
					if (mod->size_or_time == 1) { //time depedent
						comp_sum += ((*vec)[ind]->dinfo.parent->comp_time[1][0] + (*vec)[ind]->dinfo.parent->comp_time[1][1]) / 2;
					}
					else {
						comp_sum += ((*vec)[ind]->dinfo.parent->comp_size[1][0] + (*vec)[ind]->dinfo.parent->comp_size[1][1]) / 2;
					}
					sett_sum += ((*vec)[ind]->dinfo.parent->S0[1][0] + (*vec)[ind]->dinfo.parent->S0[1][1]) / 2;
				}
				else {
					if (mod->size_or_time == 1) { //time depedent
						comp_sum += ((*vec)[ind]->dinfo.parent->comp_time[0][0] + (*vec)[ind]->dinfo.parent->comp_time[0][1]) / 2;
					}
					else {
						comp_sum += ((*vec)[ind]->dinfo.parent->comp_size[0][0] + (*vec)[ind]->dinfo.parent->comp_size[0][1]) / 2;
					}
					sett_sum += ((*vec)[ind]->dinfo.parent->S0[0][0] + (*vec)[ind]->dinfo.parent->S0[0][1]) / 2;
				}
			}
		}
		values.push_back(comp_sum / vec->size());
		values.push_back(sett_sum / vec->size());

		if (pop_stageinfo[0]->s_sett.densdep == 1) {

			for (int ind = 0; ind < vec->size(); ind++) {
				if(mod->repro == 0) {
					alphaS_sum += ((*vec)[ind]->dinfo.parent->alphaS[0][0]);
					betaS_sum += ((*vec)[ind]->dinfo.parent->betaS[0][0]);
				}
				if (pop_stageinfo[0]->s_sett.sexdep == 1 && (*vec)[ind]->sex == 0) { //sex dependent
					alphaS_sum += ((*vec)[ind]->dinfo.parent->alphaS[1][0] + (*vec)[ind]->dinfo.parent->alphaS[1][1]) / 2;
					betaS_sum += ((*vec)[ind]->dinfo.parent->betaS[1][0] + (*vec)[ind]->dinfo.parent->betaS[1][1]) / 2;
				}
				else {
					alphaS_sum += ((*vec)[ind]->dinfo.parent->alphaS[0][0] + (*vec)[ind]->dinfo.parent->alphaS[0][1]) / 2;
					betaS_sum += ((*vec)[ind]->dinfo.parent->betaS[0][0] + (*vec)[ind]->dinfo.parent->betaS[0][1]) / 2;
				}
				
			}
			values.push_back(alphaS_sum / vec->size());
			values.push_back(betaS_sum / vec->size());
		}
		return values;
	}



}
//float Subpopulation::mean_evolution(string param, bool first) {
//	cout << "called mean evolution function, ";
//	float value = 0;
//	if (first == true) { //if this is year0 and therefore all individuals are in indivs and not surv_indivs!
//		cout << "in year0, there are " << indivs.size() << " surviving indivs in subpop " << patch_num->patch_ID;
//
//		for (int ind = 0; ind < indivs.size(); ind++) {
//			if (param == "emigration") {
//				//cout << " parameter was emigration, ";
//				if (pop_stageinfo[0]->s_emig.densdep == 1) { //density dependent
//					value += indivs[ind]->dinfo.D0;
//				}
//				value += indivs[ind]->dinfo.emig_prob;
//				//cout << "indiv " << indivs[ind]->ID << " has emig prob " << indivs[ind]->dinfo.emig_prob << endl;
//			}
//			else if (param == "alpha") {
//				value += indivs[ind]->dinfo.alpha;
//			}
//			else if (param == "beta") {
//				value += indivs[ind]->dinfo.beta;
//			}
//			else if (param == "min_active") {
//				if (mod->size_or_time == 1) { //time depedent
//					value += indivs[ind]->dinfo.min_active_time;
//				}
//				else {
//					value += indivs[ind]->dinfo.min_active_size;
//				}
//			}
//			else if (param == "growth") {
//				if (pop_stageinfo[0]->s_trans.pop_grow.method == 0) {
//					value += indivs[ind]->dinfo.grow_info.m;
//				}
//				else if (pop_stageinfo[0]->s_trans.pop_grow.method == 1) {
//
//				}
//				else {
//
//				}
//			}
//			else if (param == "comp") {
//				if (mod->size_or_time == 1) { //time depedent
//					value += indivs[ind]->dinfo.comp_time;
//				}
//				else {
//					value += indivs[ind]->dinfo.comp_size;
//				}
//			}
//			else if (param == "settlement") {
//				value += indivs[ind]->dinfo.S0;
//			}
//			else if (param == "alphaS") {
//				value += indivs[ind]->dinfo.alphaS;
//			}
//			else if (param == "betaS") {
//				value += indivs[ind]->dinfo.betaS;
//			}
//
//		}
//		value /= indivs.size();
//	}
//	else {
//		cout << "there are " << surv_indivs.size() << " surviving indivs in subpop " << patch_num->patch_ID;
//		for (int ind = 0; ind < surv_indivs.size(); ind++) {
//			if (param == "emigration") {
//				//cout << " parameter was emigration, ";
//				if (pop_stageinfo[0]->s_emig.densdep == 1) { //density dependent
//					value += surv_indivs[ind]->dinfo.D0;
//				}
//				value += surv_indivs[ind]->dinfo.emig_prob;
//				//cout << "indiv " << surv_indivs[ind]->ID << " has emig prob " << surv_indivs[ind]->dinfo.emig_prob << endl;
//			}
//			else if (param == "alpha") {
//				value += surv_indivs[ind]->dinfo.alpha;
//			}
//			else if (param == "beta") {
//				value += surv_indivs[ind]->dinfo.beta;
//			}
//			else if (param == "min_active") {
//				if (mod->size_or_time == 1) { //time depedent
//					value += surv_indivs[ind]->dinfo.min_active_time;
//				}
//				else {
//					value += surv_indivs[ind]->dinfo.min_active_size;
//				}
//			}
//			else if (param == "growth") {
//				value += surv_indivs[ind]->dinfo.grow_info.m;
//			}
//			else if (param == "comp") {
//				if (mod->size_or_time == 1) { //time depedent
//					value += surv_indivs[ind]->dinfo.comp_time;
//				}
//				else {
//					value += surv_indivs[ind]->dinfo.comp_size;
//				}
//			}
//			else if (param == "settlement") {
//				value += surv_indivs[ind]->dinfo.S0;
//			}
//			else if (param == "alphaS") {
//				value += surv_indivs[ind]->dinfo.alphaS;
//			}
//			else if (param == "betaS") {
//				value += surv_indivs[ind]->dinfo.betaS;
//			}
//
//		}
//		value /= surv_indivs.size();
//	}
//
//	
//	cout << "mean evo value in subpop" << patch_num->patch_ID << " is " << value << endl;
//
//	return value;
//}

//deleting individuals

void Subpopulation::delete_indivs() {
	for (int i = 0; i < indivs.size(); i++) {
		delete indivs[i];
	}
	indivs.clear();
}