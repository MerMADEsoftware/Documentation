#include "classes.h"

//random number generator definition, corresponding to extern declaration in classes.h
mt19937 eng(time(0));
//mt19937 eng(666); //for debugging, to get the same result every time
uniform_real_distribution<double> unif_dist(0, 1);
double M_PI = 3.141592653589793;
double M2_PI = 6.283185307179586;
//fill_mat

//void fill_mat333(mat333* mat, float base, float phi) {
//	
//	int drow, dcol;
//
//	int i0 = 1, i1 = base, i2 = pow(base, 2), i3 = pow(base, 3), i4 = pow(base, 4);
//
//
//	if (fabs(phi) > (7 * M_PI / 8)) { drow = 0; dcol = -1; }
//	else if (fabs(phi) > (5 * M_PI / 8)) {
//		dcol = -1;
//		if (phi > 0) { drow = 1; }
//		else { drow = -1; }
//	}
//	else if (fabs(phi) > (3 * M_PI / 8)) {
//		dcol = 0;
//		if (phi > 0) { drow = -1; }
//		else { drow = 1; }
//	}
//	else if (fabs(phi) > M_PI / 8) {
//		dcol = 1;
//		if (phi > 0) { drow = -1; }
//		else { drow = 1; }
//	}
//	else {
//		dcol = 1; drow = 0;
//	}
//
//	for (int lay = 0; lay < 3; lay++) {
//		mat->neigh[drow + 1][dcol + 1][lay] = i0;
//		mat->neigh[-drow + 1][-dcol + 1][lay] = i4;
//		if (drow == 0 || dcol == 0) { //if it points to a cardinal direction
//			mat->neigh[dcol + 1][drow + 1][lay] = i2; //need to switch them around because i am looking for the perpendicular coordinates
//			mat->neigh[-dcol + 1][-drow + 1][lay] = i2;
//
//			if (drow == 0) { //points east west
//				mat->neigh[drow][dcol + 1][lay] = i1;
//				mat->neigh[drow + 2][dcol + 1][lay] = i1;
//				mat->neigh[drow][-dcol + 1][lay] = i3;
//				mat->neigh[drow + 2][-dcol + 1][lay] = i3;
//			}
//			else { //points north south
//				mat->neigh[drow + 1][dcol][lay] = i1;
//				mat->neigh[drow + 1][dcol + 2][lay] = i1;
//				mat->neigh[-drow + 1][dcol][lay] = i3;
//				mat->neigh[-drow + 1][dcol + 2][lay] = i3;
//			}
//		}
//		else { //points to an ordinal direction
//			mat->neigh[drow + 1][-dcol + 1][lay] = i2;
//			mat->neigh[-drow + 1][dcol + 1][lay] = i2;
//			mat->neigh[drow + 1][dcol*dcol][lay] = i1;
//			mat->neigh[-drow + 1][dcol*dcol][lay] = i3;
//
//			//top diagonals
//			if (drow < 0) {
//				mat->neigh[-drow][dcol + 1][lay] = i1;
//				mat->neigh[-drow][-dcol + 1][lay] = i3;
//			}
//			else { //bottom diagonals
//				mat->neigh[drow][dcol + 1][lay] = i1;
//				mat->neigh[drow*drow][-dcol + 1][lay] = i3;
//			}
//
//
//		}
//		if (lay != 1) { //need to fill in the middle cells for the above and below layers
//			mat->neigh[1][1][lay] = i1;
//		}
//		else { mat->neigh[1][1][lay] = 0; }
//	}
//}

int main() {
	//set up model
	Model m1;
	m1.read_control(); //read control file. information in this file is the same for all batches
	
	//initialise landscape and population;
	Landscape ls1; Population pop1;
	pop1.lands = &ls1; //set the pointer to the created landscape
	ls1.mod = &m1; pop1.mod = &m1; //assign model pointers
	
	//when i start doing sims
	//for(int s=0; s<m1.nsims; s++){
		//set up landscape and population
		//for (int r = 0; r < m1.nreps; r++) {
			//run model
		//}
	//}

	cout << "created model, landscape and population" << endl;
	//set up landscape
	ls1.get_info(); //set parameters
	//cout << "number of cells in spdist is " << ls1.spdist_cells.size() << endl;

	cout << "set up landscape" << endl;
	//cout << "land dimensions: xmin " << ls1.get_land_att(9) << " xmax " << ls1.get_land_att(10) << ", ymin " << ls1.get_land_att(11) << ", ymax " << ls1.get_land_att(12) << ", zmin " << ls1.get_land_att(13) << ", zmax " << ls1.get_land_att(14) << endl;
	//cout << "there are " << size(ls1.suitable_cells) << " suitable cells" << endl;
	//ls1.fill_neigh(); //this is to test vector overheads vs 3D array in terms of memory
	//cout << "done with filling" << endl;
	//ls1.check_landscape();
	//set up Population
	pop1.set_pop(); //set up the population parameters
	
	cout << "set up population" << endl;
	
	//initialisation rules
	m1.read_init(); //reads initialisation file

	////open the output files
	ofstream population_file, indivs_file, indmov_file, range_file;
	stringstream ss, tt, qq;
	if (m1.Outpop == 1) {
		ss << m1.wrkdir + "/Outputs/Population.txt";
		population_file.open(ss.str().c_str());
		//add headers
		if (m1.patchmodel == 0) { population_file << "Rep" << "\t" << "Year" << "\t" << "Rep_Season" << "\t" << "Cell_ID" << "\t" << "Row" << "\t" << "Col" << "\t" << "Layer" << "\t" << "Patch_ID" << "\t" << "NInd" << "\t" << "NFemales" << "\t" << "NMales"; }
		else { population_file << "Rep" << "\t" << "Year" << "\t" << "Rep_Season" << "\t" << "Cell_ID" << "\t" << "Patch_ID" << "\t" << "NInd" << "\t" << "NFemales" << "\t" << "NMales"; }
		if (pop1.pop_dyn.pop_evo.ep_evo == 1) { //if emigration rate evolves
			population_file << "\t" << "Mean_EP";
			if (pop1.stages_b[0]->s_emig.densdep == 1) { //if it's density dependent
				population_file << "\t" << "Mean_alpha" << "\t" << "Mean_beta";
			}
		}
		if (pop1.pop_dyn.pop_evo.active_evo == 1) { //if active time/size evolves
			if (m1.size_or_time == 1) { //time-dependent
				population_file << "\t" << "Mean_active_time";
			}
			else { //size-dependent
				population_file << "\t" << "Mean_active_size";
			}
		}
		if (pop1.pop_dyn.pop_evo.dv_evo == 1) { //if active time/size evolves
			if (m1.size_or_time == 1) { //time-dependent
				population_file << "\t" << "Mean_dv_time";
			}
			else { //size-dependent
				population_file << "\t" << "Mean_dv_size";
			}
		}
		if (pop1.pop_dyn.pop_evo.growth_evo == 1) { //if growth rate evolves
			if (m1.grow_method == "linear") {
				population_file << "\t" << "Mean_growth_rate";
			}
			else if (m1.grow_method == "gompertz") {
				population_file << "\t" << "Mean_Linf" << "\t" << "Mean_K" << "\t" << "Mean_Ti";
			}
			else {

			}
			
		}
		if (pop1.pop_dyn.pop_evo.comp_evo == 1) { //if minimum competency time/size evolves
			if (m1.size_or_time == 1) { //time dependent
				population_file << "\t" << "Mean_comp_time";
			}
			else { //size dependent
				population_file << "\t" << "Mean_comp_size";
			}
		}
		if (pop1.pop_dyn.pop_evo.S0_evo == 1) { //if settlement rate evolves
			population_file << "\t" << "Mean_S0";
			if (pop1.stages_b[0]->s_sett.densdep == 1) { //if settlement is density dependent
				population_file << "\t" << "Mean_alphaS" << "\t" << "Mean_betaS";
			}
		}
		population_file << endl;
		population_file.close();
	}
	if (m1.Outind == 1) {
		tt << m1.wrkdir + "/Outputs/Individuals.txt";
		indivs_file.open(tt.str().c_str());
		if (m1.patchmodel == 0) { indivs_file << "Rep" << "\t" << "Year" << "\t" << "Rep_Season" << "\t" << "Ind_ID" << "\t" << "Status" << "\t" << "Natal_X" << "\t" << "Natal_Y" << "\t" << "Natal_Z" << "\t" << "X" << "\t" << "Y" << "\t" << "Z" << "\t" << "DistMoved"; }
		else { indivs_file << "Rep" << "\t" << "Year" << "\t" << "Rep_Season" << "\t" << "Ind_ID" << "\t" <<"Sex" << "\t" << "Status" << "\t" << "Starting_patch" << "\t" << "Natal_X" << "\t" << "Natal_Y" << "\t" << "Natal_Z" << "\t" << "Ending_patch" << "\t" << "X" << "\t" << "Y" << "\t" << "Z" << "\t" << "DistMoved"; }
		if (m1.size_or_time == 0) {
			indivs_file <<"\t" << "Size" << "\t";
		}
		if (pop1.pop_dyn.pop_evo.ep_evo == 1) { //if emigration rate evolves
			indivs_file << "\t" << "EP";
			if (pop1.stages_b[0]->s_emig.densdep == 1) { //if it's density dependent
				indivs_file << "\t" << "alpha" << "\t" << "beta";
			}
		}
		if (pop1.pop_dyn.pop_evo.active_evo == 1) { //if active time/size evolves
			if (m1.size_or_time == 1) { //time-dependent
				indivs_file << "\t" << "Min_active_time";
			}
			else { //size-dependent
				indivs_file << "\t" << "Min_active_size";
			}
		}
		if (pop1.pop_dyn.pop_evo.dv_evo == 1) {
			if (m1.size_or_time == 1) { //time-dependent
				indivs_file << "\t" << "Min_dv_time";
			}
			else { //size-dependent
				indivs_file << "\t" << "Min_dv_size";
			}
		}
		if (pop1.pop_dyn.pop_evo.growth_evo == 1) { //if growth rate evolves
			if (m1.grow_method == "linear") {
				indivs_file << "\t" << "Growth_rate_m";
			}
			else if (m1.grow_method == "gompertz") {
				indivs_file << "\t" << "Linf" << "\t" << "Growth_K" << "\t" << "Ti";
			}
			else {

			}

		}
		if (pop1.pop_dyn.pop_evo.comp_evo == 1) { //if minimum competency time/size evolves
			if (m1.size_or_time == 1) { //time dependent
				indivs_file << "\t" << "Min_comp_time";
			}
			else { //size dependent
				indivs_file << "\t" << "Min_comp_size";
			}
		}
		if (pop1.pop_dyn.pop_evo.S0_evo == 1) { //if settlement rate evolves
			indivs_file << "\t" << "S0";
			if (pop1.stages_b[0]->s_sett.densdep == 1) { //if settlement is density dependent
				indivs_file << "\t" << "alphaS" << "\t" << "betaS";
			}
		}
		indivs_file << endl;
		indivs_file.close();
	}
	
	/*if (m1.Outindmov == 1) {
		ii << m1.wrkdir + "/Outputs/Indiv_Mov.txt";
		cout << "opening indiv mov file " << m1.wrkdir + "/Outputs/Indiv_Mov.txt" << endl;
		indmov_file.open(ii.str().c_str());
		indmov_file << "Rep" << "\t" << "Year" << "\t" << "Rep_Season" << "\t" << "Ind_ID" << "\t" << "Sex" << "\t" << "Status" << "\t" << "XYZ_coords" << endl;
	}*/
	if (m1.Outrange == 1) {
		qq << m1.wrkdir + "/Outputs/Range.txt";
		range_file.open(qq.str().c_str());
		range_file << "Rep" << "\t" << "Year" << "\t" << "Rep_Season" << "\t" << "NOcc" << "\t" << "OccupSuit" << "\t" << "Min_X" << "\t" << "Max_X" << "\t" << "Min_Y" << "\t" << "Max_Y" << "\t" << "Min_Z" << "\t" << "Max_Z" << endl;
		range_file.close();
	}

	//set up a for loop for how many years , how many rep events, when does dispersal happen etc
	for (int r = 0; r < m1.nreps; r++) { //for each rep
		cout << "rep " << r << " ";
		stringstream ii;
		if (m1.Outindmov == 1) {
			ii << m1.wrkdir + "/Outputs/Indiv_Mov_rep" << r << ".txt";
			cout << "opening indiv mov file " << m1.wrkdir + "/Outputs/Indiv_Mov.txt" << endl;
			indmov_file.open(ii.str().c_str());
			indmov_file << "Rep" << "\t" << "Year" << "\t" << "Rep_Season" << "\t" << "Ind_ID" << "\t" << "Sex" << "\t" << "Status" << "\t" << "XYZ_coords" << endl;
			indmov_file.close();
		}
	//initialise subpopulations and individuals
		pop1.initialise(); //takes initialisation information from model, creates subpopulations, calls init_indivs for each subpop
		cout << "starting pop is this size: " << pop1.pop_dyn.total_births << endl;
		cout << "initialised population" << endl;
		/*pop1.write_popfile(population_file, r, 0, 0);
		pop1.write_indivfile(indivs_file, r, 0, 0);
		pop1.write_rangefile(range_file, r, 0, 0);*/
		int nchanges = 0; //this is to keep track of how many dynamic changes have been made so far
		for (int y = 0; y < m1.nyears; y++) { //for each year
			cout << "year " << y << "...";
			if (m1.dynamic == 1 && nchanges < ls1.dynamic_info.num) {
				//if (y == ls1.dynamic_info.burn_in) {

				//}
				//else if (y > ls1.dynamic_info.burn_in && (y - ls1.dynamic_info.burn_in) % ls1.dynamic_info.dyn_interval == 0) {

				//}
				if (y == ls1.dynamic_info.burn_in || y>ls1.dynamic_info.burn_in && (y-ls1.dynamic_info.burn_in)%ls1.dynamic_info.dyn_interval == 0) { //if we've past the burn-in period and are in the right interval
					ls1.do_dyn(nchanges); //then read in the new dynamic seascape hab/patch/hydrodynamics
					nchanges++; //now increase because we've made a change
				}
				
			}
			pop1.one_year(population_file, indivs_file, indmov_file, range_file, r, y); //reproduction, dispersal (emigration, transfer, settlement), survival, aging and development
			cout << "done" << endl;
			//if (ls1.patch_extinct==1 && y%ls1.p_ext_int == 0) {
			//	cout << " local patch extinction" << endl;
			//	vector<Individual*> proxy;
			//	float ran_i_ext = unif_dist(eng);
			//	if (ls1.p_ext_method == 0) { //random proportion
			//		float ran_p_ext = unif_dist(eng);
			//		for (int p = 0; p < pop1.subpops.size(); p++) {						
			//			if (pop1.subpops[p]->size > 0 && ran_p_ext < ls1.p_ext_p_prop) { //if the subpop isn't empty, and it is going to go extinct							
			//				for (int i = 0; i < pop1.subpops[p]->indivs.size(); i++) {
			//					if (ran_i_ext < ls1.p_ext_i_prop) {
			//						delete pop1.subpops[p]->indivs[i];
			//					}
			//					else {
			//						proxy.push_back(pop1.subpops[p]->indivs[i]);
			//					}
			//				}
			//				pop1.subpops[p]->indivs.clear();
			//				pop1.subpops[p]->indivs = proxy;
			//				proxy.clear();
			//			}
			//		}
			//	}
			//	else { //specific patch					
			//		if (pop1.subpops[ls1.p_ext_patch]->size > 0) { //if the subpop isn't empty, and it is going to go extinct
			//			cout << "clearing patch " << ls1.p_ext_patch << endl;
			//			for (int i = 0; i < pop1.subpops[ls1.p_ext_patch]->indivs.size(); i++) {
			//				if (ran_i_ext < ls1.p_ext_i_prop) {
			//					delete pop1.subpops[ls1.p_ext_patch]->indivs[i];
			//				}
			//				else {
			//					proxy.push_back(pop1.subpops[ls1.p_ext_patch]->indivs[i]);
			//				}
			//			}
			//			pop1.subpops[ls1.p_ext_patch]->indivs.clear();
			//			pop1.subpops[ls1.p_ext_patch]->indivs = proxy;
			//			proxy.clear();
			//		}
			//	}
			//}
		}

		pop1.reset();
		//cout << "after pop reset, total_births is " << pop1.pop_dyn.total_births << endl;
		indmov_file.close();
	}

	population_file.close(); indivs_file.close(); range_file.close();

	//delete everything that was created with new
	//subpopulations
	//pop1.delete_subpops(); //within each subpop, indivs are deleted. indivs in the matrix are also deleted within this function
	//cells
	ls1.delete_raster();
	//patches
	ls1.delete_patch_vector();
	cout << "done" << endl;
	int done;
	cin >> done;
	cin.get();
	return 0;
}