#include "classes.h"

Cell::Cell() {
	K = -9;
	habtype = -9;
}

Cell::Cell(Model* modp) {
	K = -9;
	habtype = -9;
	mod = modp;

	//initialise structs to -9
	/*for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				costs.neigh[i][j][k] = -9999;
				current.neigh[i][k][k] = -9999;
			}
		}
	}*/

	//initialise this to null
	//for (int r = 0; r < 3; r++) {
	//	for (int c = 0; c < 3; c++) {
	//		for (int l = 0; l < 3; l++) {
	//			ow_neigh[r][c][l] = nullptr;
	//		}
	//	}
	//}
}

Cell::~Cell()
{
}

bool Cell::check_limits(int minX, int maxX, int minY, int maxY, int minZ, int maxZ) {
	//if the cell is within the input range
	//remember x,y,z are the col, row, layer values 
	if (c >= minX && c <= maxX &&
		r >= minY && r <= maxY &&
		l >= minZ && l <= maxZ) {
		//cout << " cell is in init range" << endl;
		return true;
	}
	else {
		return false;
	}
}

void Cell::calc_midpoint() {

	float mid_x = cell_corners.llc[0] + (res / 2);
	float mid_y = cell_corners.ulc[1] - (res / 2);

	midpoint = { mid_x, mid_y, min_depth };
}


void Cell::calc_maxinds() { //K is the carrying capacity of the habitat type
	max_inds = floor(K * ha);//I also round down to make it a whole individual
}

void Cell::find_ow_neigh(Landscape* lands) {

	/*struct mat_index { int row, col, layer; };

	std::deque<mat_index> possibilities;*/

	if ((r - 1) >= 0 && lands->raster_3d[r - 1][c][l] != nullptr && lands->raster_3d[r - 1][c][l]->habtype == 0) { //north
		mat_index ind;
		ind.row = r - 1; ind.col = c; ind.layer = l; //record that location
		possibilities.push_back(ind); //and add it to the list of possibilities
	}

	if ((r + 1) < lands->get_land_att(1) && lands->raster_3d[r + 1][c][l] != nullptr && lands->raster_3d[r + 1][c][l]->habtype == 0) { //south
		mat_index ind;
		ind.row =r + 1; ind.col = c; ind.layer = l; //record that location
		possibilities.push_back(ind); //and add it to the list of possibilities
	}

	if ((c - 1) >= 0 && lands->raster_3d[r][c - 1][l] != nullptr && lands->raster_3d[r][c - 1][l]->habtype == 0) { //west
		mat_index ind;
		ind.row = r; ind.col = c - 1; ind.layer = l; //record that location
		possibilities.push_back(ind); //and add it to the list of possibilities
	}

	if ((c + 1) < lands->get_land_att(2) && lands->raster_3d[r][c + 1][l] != nullptr && lands->raster_3d[r][c + 1][l]->habtype == 0) { //east
		mat_index ind;
		ind.row = r; ind.col = c + 1; ind.layer = l; //record that location
		possibilities.push_back(ind); //and add it to the list of possibilities
	}

	if ((l - 1) >= 0 && lands->raster_3d[r][c][l - 1] != nullptr && lands->raster_3d[r][c][l - 1]->habtype == 0) { //up
		mat_index ind;
		ind.row = r; ind.col = c; ind.layer = l - 1; //record that location
		possibilities.push_back(ind); //and add it to the list of possibilities
	}
	if (possibilities.size() == 0) { cout << "ran find_ow_neighbours but found none at "<< r << ", " << c << ", " <<l << endl; }
}

void Cell::assess_pos_buffer(int detect_dist, int det_layers, Landscape* lands) {
	//because the distance to the suitable cell matters, i have to check the inner circle first before moving outwards
	//inner circle
	bufferpos = 0;
	int min_row = r - 1;
	int max_row = r + 1;
	int min_col = c - 1;
	int max_col = c + 1;
	int min_layer = l - det_layers;
	int max_layer = l + det_layers; 
	for (int layer = l+1; layer > (l-2); layer--) {
		for (int row = min_row; row < (max_row+1); row++) {
			for (int col = min_col; col < (max_col+1); col++) {
				if (row == r && col == c && layer == l) { continue; } //skip the current cell
				if (row >= lands->get_land_att(1) || row<0 || col>= lands->get_land_att(2) || col < 0 || layer<0 || layer >= lands->get_land_att(5)) {continue; } //if it's outside the landscape boundaries, just skip it, don't assess it
				//cout << "row " << row << ", col" << col << ", layer " << layer << endl;
				if (lands->raster_3d[row][col][layer] == nullptr) {  continue; }

				if (lands->raster_3d[row][col][layer]->K > 0) {
					if (mod->patchmodel == 1) {
						if (lands->raster_3d[row][col][layer]->patch==nullptr) { //if this cell is not part of any patch
							continue; //skip it. in patch-based models, only suitable cells within patches can be settled on
						}
						
					}
					bufferpos = 1; 
					buffer_focus.push_back(lands->raster_3d[row][col][layer]);

				}
				//else if (lands->raster_3d[row][col][layer]->habtype == 0) {
				//	if (layer < l) {
				//		open_water_neigh.up_neigh.push_back(lands->raster_3d[row][col][layer]);
				//	}
				//	else if (layer == l) {
				//		open_water_neigh.forward_neigh.push_back(lands->raster_3d[row][col][layer]);
				//	}
				//	else {
				//		open_water_neigh.down_neigh.push_back(lands->raster_3d[row][col][layer]); //open water neighbours are now organised by layer below to layer above

				//	}
				//}
			}
		}
	}
	if (detect_dist == 1) { return; } //even if there is a goal bias, dont continue calculating

	//further out circles
	else {
		for (int d = 2; d < detect_dist; d++) {
			for (int layer = min_layer; layer < (max_layer+1); layer++) {
				if (layer == lands->get_land_att(5) || layer<0) { continue; } //if the layer is too deep or above the surface, skip;
				min_row = r - d; max_row = r + d;
				min_col = c - d; max_col = c + d;
				for (int row = min_row; row < (max_row + 1); row++) {
					if (row >= lands->get_land_att(1) || row < 0) { continue; } //if the row is outside the landscape, just skip it
					
					if (row == min_row || row==max_row) { //if it's the top or bottom row, i want to do all the columns
						for (int col = min_col; col < (max_col + 1); col++) {
							if (col>= lands->get_land_att(2) || col < 0) { continue; } //if it's outside the landscape boundaries, just skip it, don't assess it
							if (lands->raster_3d[row][col][layer] == nullptr) { continue; }
							if (lands->raster_3d[row][col][layer]->K > 0) {
								if (mod->patchmodel == 1) {
									if (!lands->raster_3d[row][col][layer]->patch) { //if this cell is not part of any patch
										continue; //skip it. in patch-based models, only suitable cells within patches can be settled on
									}
								}
								bufferpos = 1; buffer_focus.push_back(lands->raster_3d[row][col][layer]);
							}
							
						}
					}
					else {
						//check the left-side neighbour
						if ( (c - d) < lands->get_land_att(2) && (c - d) > 0 ) { //if it's inside the landscape boundaries, assess the cells
							if (lands->raster_3d[row][c-d][layer] == nullptr) { continue; }

							if (lands->raster_3d[row][c-d][layer]->K > 0) {
								if (mod->patchmodel == 1) {
									if (!lands->raster_3d[row][c-d][layer]->patch) { //if this cell is not part of any patch
										continue; //skip it. in patch-based models, only suitable cells within patches can be settled on
									}
								}
								bufferpos = 1; buffer_focus.push_back(lands->raster_3d[row][c-d][layer]);
							}
						}
						//and the right-side neighbour
						if ( (c + d) < lands->get_land_att(2) && (c + d) > 0) { //if it's inside the landscape boundaries, assess the cells
							if (lands->raster_3d[row][c+d][layer] == nullptr) { continue; }

							if (lands->raster_3d[row][c + d][layer]->K > 0) {
								if (mod->patchmodel == 1) {
									if (!lands->raster_3d[row][c+d][layer]->patch) { //if this cell is not part of any patch
										continue; //skip it. in patch-based models, only suitable cells within patches can be settled on
									}
								}
								bufferpos = 1; buffer_focus.push_back(lands->raster_3d[row][c + d][layer]);
							}
						}
					}
				}
				
			}
			//if (bufferpos == 1) { return; } // if i have now found suitable cells in the outer circles, dont continue
		}

	}
}

//void Cell::assess_neg_buffer(int detect_dist, int det_layers, int neg, Landscape* lands) {
//	//because the distance to the suitable cell matters, i have to check the inner circle first before moving outwards
//	//inner circle
//	bufferneg = 0;
//	int min_row = r - 1;
//	int max_row = r + 1;
//	int min_col = c - 1;
//	int max_col = c + 1;
//	int min_layer = l - det_layers;
//	int max_layer = l + det_layers;
//
//	for (int layer = min_layer; layer < (max_layer + 1); layer++) {
//		for (int row = min_row; row < (max_row+1); row++) {
//			for (int col = min_col; col < (max_col + 1); col++) {
//				if (row == r && col == c && layer == l) { continue; } //skip the current cell
//				if (row >= lands->get_land_att(1) || row < 0 || col >= lands->get_land_att(2) || col < 0 || layer<0 || layer >= lands->get_land_att(5)) { continue; } //if it's outside the landscape boundaries, just skip it, don't assess it
//				if (lands->raster_3d[row][col][layer] == nullptr) { continue; }
//				cout << "looking at cell " << row << ", " << col << ", " << layer; // cout << "which has habtype " << lands->raster_3d[row][col][layer]->habtype << " and is part of patch " << lands->raster_3d[row][col][layer]->patch->patch_ID << endl;
//				if (mod->patchmodel == 0 && lands->raster_3d[row][col][layer]->cell_ID == neg) {
//					bufferneg = 1; buffer_natal.push_back(lands->raster_3d[row][col][layer]);
//				}
//				else if (mod->patchmodel == 1 && lands->raster_3d[row][col][layer]->patch!=nullptr && lands->raster_3d[row][col][layer]->patch->patch_ID == neg) {
//					bufferneg = 1; buffer_natal.push_back(lands->raster_3d[row][col][layer]);
//				}
//
//			}
//		}
//	}
//	if (detect_dist == 1) { return; }
//	
//	else { //further out circles
//		for (int d = 2; d < detect_dist; d++) {
//			for (int layer = min_layer; layer < max_layer; layer++) {
//				if (layer >= lands->get_land_att(5) || layer<0) { continue; } //if the layer is too deep, or above surface, skip;
//				min_row = r - d; max_row = r + d;
//				min_col = c - d; max_col = c + d;
//				for (int row = min_row; row < max_row + 1; row++) {
//					if (row >= lands->get_land_att(1) || row < 0) { continue; } //if the row is outside the landscape, just skip it
//
//					if (row == min_row || row == max_row) { //if it's the top or bottom row, i want to do all the columns
//						for (int col = min_col; col < max_col + 1; col++) {
//							if (col >= lands->get_land_att(2) || col < 0) { continue; } //if it's outside the landscape boundaries, just skip it, don't assess it
//							if (lands->raster_3d[row][col][layer] == nullptr) { continue; }
//
//							if (mod->patchmodel == 0 && lands->raster_3d[row][col][layer]->cell_ID == neg) {
//								bufferneg = 1; buffer_natal.push_back(lands->raster_3d[row][col][layer]);
//							}
//							else if (mod->patchmodel == 1 && lands->raster_3d[row][col][layer]->patch && lands->raster_3d[row][col][layer]->patch->patch_ID == neg) {
//								bufferneg = 1; buffer_natal.push_back(lands->raster_3d[row][col][layer]);
//							}
//
//						}
//					}
//					else {
//						//check the left-side neighbour
//						if (c - d < lands->get_land_att(2) && c - d > 0) { //if it's inside the landscape boundaries, assess the cells
//							if (lands->raster_3d[row][c - d][layer] == nullptr) { continue; }
//
//							if (mod->patchmodel == 0 && lands->raster_3d[row][c - d][layer]->cell_ID == neg) {
//								bufferneg = 1; buffer_natal.push_back(lands->raster_3d[row][c - d][layer]);
//							}
//							else if (mod->patchmodel == 1 && lands->raster_3d[row][c - d][layer]->patch && lands->raster_3d[row][c - d][layer]->patch->patch_ID == neg) {
//								bufferneg = 1; buffer_natal.push_back(lands->raster_3d[row][c - d][layer]);
//							}
//						}
//						//and the right-side neighbour
//						if (c + d < lands->get_land_att(2) && c + d > 0) { //if it's inside the landscape boundaries, assess the cells
//							if (lands->raster_3d[row][c + d][layer] == nullptr) { continue; }
//
//							if (mod->patchmodel == 0 && lands->raster_3d[row][c + d][layer]->cell_ID == neg) {
//								bufferneg = 1; buffer_natal.push_back(lands->raster_3d[row][c + d][layer]);
//							}
//							else if (mod->patchmodel == 1 && lands->raster_3d[row][c + d][layer]->patch && lands->raster_3d[row][c + d][layer]->patch->patch_ID == neg) {
//								bufferneg = 1; buffer_natal.push_back(lands->raster_3d[row][c + d][layer]);
//							}
//						}
//					}
//				}
//
//			}
//			if (bufferneg == 1) { return; } // if i have now found suitable cells in the outer circles, dont continue
//		}
//
//	}
//}



//float Cell::find_mean_cost(int minrow, int maxrow, int mincol, int maxcol, int minlay, int maxlay, Landscape* lands) {
//	float cost_sum = 0;
//	int cost_counter = 0;
//	float cost_mean=-9;
//	
//	//mean calculations are only done one layer at a time, the mean is not taking into acount all layers within z_PR!
//	for (int c_l = minlay; c_l < (maxlay + 1); c_l++) {
//
//		if (c_l < 0) { cout << "layer above surface, skipping" << endl; continue; }
//		else if (c_l >= lands->get_land_att(5)) { continue; }
//		for (int c_r = minrow; c_r < (maxrow + 1); c_r++) {
//
//			for (int c_c = mincol; c_c < (maxcol + 1); c_c++) {
//
//				if (c_r < 0 || c_r == lands->get_land_att(1) || c_c < 0 || c_c == lands->get_land_att(2)) { 
//					cout << "out of bounds" << endl; cost_sum += 20;  //if it's outside of the xy landscape boundaries, treat it like open water
//					cost_counter++; continue; 
//				}
//				else if (c_l < 0 || c_l == lands->get_land_att(5)) {
//					cost_sum += 50; //if it's the depth boundaries, treat it like unsuitable habitat
//					cost_counter++; continue;
//				}
//				if (lands->raster_3d[c_r][c_c][c_l] == nullptr) { cout << "null" << endl; continue; } //if it's a no data cell
//
//				cost_sum += lands->raster_3d[c_r][c_c][c_l]->cost; 
//				if (lands->raster_3d[c_r][c_c][c_l]->K > 0) {
//					if (mod->patchmodel == 1 && lands->raster_3d[c_r][c_c][c_l]->patch == nullptr) { //if it's suitable but not in a patch for a patch-based model, this doesn't count!
//						cost_sum += 50; //treat it like unsuitable habitat
//						cost_counter++;
//						continue;
//					}
//
//				}
//				
//				cost_counter++;
//				
//			}
//		}
//		if (cost_counter != 0) { cost_mean = cost_sum / cost_counter; }
//		/*cout << "cost mean is " << cost_mean << endl;
//		if (cost_mean == 0) { cost_mean = -9; }*/ //reset it to 1 if there was really no cell in that direction at all (i know this is the same as having suitable habitat everywhere but i need it to be 1 so it makes no diff as a multiplier!
//		//note that this is different than it being the surface! i need the surface to have probability=0 while indivs can still travel out of bounds
//	}
//	return cost_mean;
//}



//need to remember to check whether this has already been calculated when i call this function!
//vector<mat333> Cell::calc_weightings(int xy_pr, int z_pr, float rho, Landscape* lands) { //i need both the horizontal and the vertical perceptual ranges because they will most likely be different resolutions!
//	//fill cost matrix using 
//	//Remember: the matrix is 3x3x3, layer0 is top, layer 1 is current layer, layer 2 is layer below
//	vector<mat333> weightings;
//
//	int rowmin, rowmax, colmin, colmax, laymin, laymax;
//	//calculate which cells are needed to estimate effective cost and goal bias
//	for (int c_l = 0; c_l < 3; c_l++) {
//		for (int c_r = -1; c_r < 2; c_r++) {
//			for (int c_c = -1; c_c < 2; c_c++) {
//				if (c_l == 1 && c_r == 0 && c_c == 0) { continue; } //skip calculations for current cell
//				else if(l == 0 && c_l == 0) { //if the current cell is in the top layer, dont even calculate the layer above it
//					costs.neigh[c_r+1][c_c+1][c_l] = -9;
//					continue;
//				}
//				if (xy_pr % 2 != 0) {//for an odd PR
//					if (c_r < 0 || c_r > 0) { //upper or lower row
//						if (c_r < 0) { rowmin = r - xy_pr; rowmax = r - 1; } //upper
//						else if (c_r > 0) { rowmin = r + 1; rowmax = r + xy_pr; } //lower
//						if (c_c < 0) { colmin = c - xy_pr; colmax = c - 1; } //left
//						else if (c_c == 0) { colmin = c - (floor(xy_pr / 2)); colmax = c + (floor(xy_pr / 2)); } //right
//						else { colmin = c + 1; colmax = c + xy_pr; } //middle
//					}
//					else { //middle row
//						//this is only c-1 and c+1 since c=0 is the current cell!
//						rowmin = r - floor(xy_pr / 2); rowmax = r + floor(xy_pr / 2);
//						if (c_c < 0) { colmin = c - xy_pr; colmax = c - 1; } //left
//						if (c_c > 0) { colmin = c + 1; colmax = c + xy_pr; } //right
//					}
//				}
//				else { //for an even PR
//					if (c_r == 0) {
//						rowmin = r - floor(xy_pr / 2); rowmax = r + floor(xy_pr / 2);
//						if (c_c < 0) { colmin = c - xy_pr; colmax = c - 1; }
//						else if (c_c > 0) { colmin = c + 1; colmax = c + xy_pr; }
//						//dont need c==0 because that would be current cell
//					}
//					else if (c_r < 0) {
//						rowmin = r - xy_pr;
//						if (c_c < 0 || c_c>0) {
//							rowmax = r;
//							if (c_c < 0) { colmin = c - xy_pr; colmax = c - 1; }
//							else { colmin = c + 1; colmax = c + xy_pr; }
//						}
//						else {
//							rowmax = r - 1;
//							colmin = c - floor(xy_pr / 2); colmax = c + floor(xy_pr / 2);
//						}
//
//					}
//					else {
//						rowmax = r + xy_pr;
//						if (c_c < 0 || c_c > 0) {
//							rowmin = r;
//							if (c_c < 0) { colmin = c - xy_pr; colmax = c - 1; }
//							else { colmin = c + 1; colmax = c + xy_pr; }
//						}
//						else {
//							rowmin = r + 1;
//							colmin = c - floor(xy_pr / 2); colmax = c + floor(xy_pr / 2);
//						}
//					}
//				}
//
//				if (c_l == 1) { laymin = l; laymax = l;
//				} //if it's the current layer, means wont be multilayer
//				else if(c_l==0) { //layer above
//					laymin = l - z_pr;
//					laymax = l - 1;
//				}
//				else {
//					laymin = l + 1;
//					laymax = l + z_pr;
//				}
//
//				if (laymin < 0) { laymin = 0; }
//				if (laymax < 0) { laymax = 0; }
//				else if (laymax == lands->get_land_att(5)) { laymax = lands->get_land_att(5) - 1; } //if laymax has gone beyond the max layer
//
//				costs.neigh[c_r+1][c_c+1][c_l] = find_mean_cost(rowmin, rowmax, colmin, colmax, laymin, laymax, lands);
//			}
//		}
//		
//	}
//
//	
//	//calculate current weights
//	double M_PI = M2_PI / 2;
//	int drow, dcol;
//	int scaled_rho = rho * 10;
//	fill_mat333(&current, scaled_rho, cell_angle);
//	
//
//	weightings.push_back(costs);
//	weightings.push_back(current);
//
//	return weightings;
//}