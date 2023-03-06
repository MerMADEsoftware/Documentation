#include "classes.h"

Patch::Patch() {
	//current_density = 0;
	in_spdist = false;
	patch_info.max_inds = 0;
	patch_info.K = 0;
	indepth_info.K = 0;
	indepth_info.max_inds = 0;
}

Patch::~Patch()
{
}

void Patch::set_limits(int row, int column,  int depth, int dint) { //remember that x=cols, y=rows
	if (included_cells.size() == 1) { //if this was the first cell to be added to that patch
		patchlimits.xmin = column; patchlimits.xmax = column; 
		patchlimits.ymin = row; patchlimits.ymax = row;
		patchlimits.zmin = dint*depth; patchlimits.zmax = dint*depth + dint; //need to add the depth interval, because otherwise this is just showing the top slice
	}
	else {
		if (column < patchlimits.xmin) { patchlimits.xmin = column; }
		if (column > patchlimits.xmax) { patchlimits.xmax = column; }
		if (row < patchlimits.ymin) { patchlimits.ymin = row; }
		if (row > patchlimits.ymax) { patchlimits.ymax = row; }
		if (dint*depth < patchlimits.zmin) { patchlimits.zmin = dint*depth; }
		if ((dint*depth +dint) > patchlimits.zmax) { patchlimits.zmax = dint*depth + dint; } //again, need it here 
	}
	
}

bool Patch::check_limits(int minX, int maxX, int minY, int maxY, int minZ, int maxZ) { //used for initialisation
	//the patch needs to be completely within the input limits in order for this function to return true
	if (minX <= patchlimits.xmin && patchlimits.xmax <= maxX && 
		minY <= patchlimits.ymin && patchlimits.ymax <= maxY /*&&
		minZ <= patchlimits.zmin && patchlimits.zmax <= maxZ*/) { //if the depth is within the patch limits

		return true;
	}
	else {
		return false;
	}
}

void Patch::calc_K(string whichstruct){
	//cout << "patch " << p_subpop << " has " << indepth_cells.size() << " in depth cells, and max inds is ";
	if (whichstruct == "patch_info") {
		float temp_K = 0;
		float temp_size_ha = 0, temp_max_inds = 0;
		for (int i = 0; i < included_cells.size(); i++) {
			temp_K += included_cells[i]->K; //add up all the Ks (K is in inds/ha)
			temp_size_ha += included_cells[i]->ha; //add up all the settle-able area in the patch
			temp_max_inds += included_cells[i]->max_inds; //add up all the max inds
		}
		patch_info.size_ha = temp_size_ha;
		patch_info.max_inds = temp_max_inds;
		patch_info.K = temp_K / included_cells.size(); //take the mean inds/ha carrying capacity
		//patch_info.max_inds = floor(patch_info.K * patch_info.size_ha); //how many has are in the patch? rounded down to nearest whole individual
	}
	else if (whichstruct == "indepth_info") {
		float temp_K = 0;
		float temp_size_ha = 0, temp_max_inds = 0;
		//cout << "there are " << indepth_cells.size() << " cells in the right depth" << endl;
		for (int i = 0; i < indepth_cells.size(); i++) {
			temp_K += indepth_cells[i]->K; //add up all the Ks (K is in inds/ha)
			temp_size_ha += included_cells[i]->ha; //add up all the settle-able area in the patch
			temp_max_inds += included_cells[i]->max_inds; //max inds will be the sum 
		}
		indepth_info.size_ha = temp_size_ha;
		indepth_info.max_inds = temp_max_inds;
		indepth_info.K = temp_K / indepth_cells.size(); //take the mean inds/ha carrying capacity
		//indepth_info.max_inds = floor(indepth_info.K * indepth_info.size_ha); //how many has are in the patch? rounded down to nearest whole individual
	}
	//cout << indepth_info.max_inds << endl;
	
}


