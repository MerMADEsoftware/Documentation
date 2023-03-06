#include "classes.h"


Landscape::Landscape() {
	cell_based = false;
	patch_based = false;
	sp_dist = true;
	npatches = 0;
	nlayers = 1;
	
}

//destructor
Landscape::~Landscape() {
}

//get functions

float Landscape::get_land_att(int x) {
	float value;
	switch (x) {

	case 1:
		value = nrows;
		break;
	case 2:
		value = ncols;
		break;
	case 3:
		value = ncells;
		break;
	case 4:
		value = npatches;
		break;
	case 5:
		value = nlayers;
		break;
	case 6:
		value = nhabitats;
		break;
	case 7:
		value = resolution;
		break;
	case 8:
		value = dint;
		break;
	case 9: 
		value = x_min;
		break;
	case 10: 
		value = y_min;
		break;
	case 11:
		value = x_max;
		break;
	case 12:
		value = y_max;
		break;
	case 13:
		value = z_min;
		break;
	case 14:
		value = z_max;
		break;

	default:
		cout << "that is not a valid atttribute number" << endl;
		break;
	}
	return value;
}

//create landscape functions
void Landscape::create_raster() { 
	cout << "landscape has " << nrows << " rows and " << ncols << " columns " << endl;
	raster_3d = new Cell***[nrows];
	for (int i = 0; i < nrows; i++) {
		
		raster_3d[i] = new Cell **[ncols];
		for (int j = 0; j < ncols; j++) {
			raster_3d[i][j] = new Cell *[nlayers]; //don't fill them with an object because i don't know which cells have data and which dont
		}
	}
	cout << "created landscape" << endl;
}

void Landscape::create_patch_vector() {
	for (int i = 0; i < npatches; i++) {
		Patch* pPatch = new Patch();
		patch_vector.push_back(pPatch);
		patch_vector[i]->patch_ID = i;
	}
}

void Landscape::delete_patch_vector() {
	for (int p = 0; p < npatches; p++) {
		delete patch_vector[p];
	}
}

void Landscape::delete_raster() {

	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++) {
			for (int k = 0; k < nlayers; k++) {
				delete raster_3d[i][j][k];
			}
			delete raster_3d[i][j]; //deleting what they are pointing to
		}
		delete raster_3d[i];
	}

	delete[]raster_3d; //delete dynamic array itself

}


//fill information into the landscape 

void Landscape::fill_md() {
	float starting_x = x_min; //this information came from the info.txt file
	float starting_y = y_max; //using y_max means you start upper left corner

	vector<float> xes; //this gives all the different xes
	for (int i = 0; i < ncols; i++) {
		if (i == 0) { xes.push_back(starting_x); }
		else { xes.push_back((xes[i - 1]) + resolution); } //adds the size of a cell to the previous coordinate to get the next coordinate
	}

	vector<float> ys; //this gives all the different ys
	for (int i = 0; i < nrows; i++) {
		if (i == 0) { ys.push_back(starting_y); }
		else { ys.push_back((ys[i - 1]) - resolution); } //because you started at max_y you need to subtract to get to min_y
	}

	vector<float> ds; //this gives all the different MINIMUM ds for each layer
	for (int i = 0; i < nlayers; i++) { // i want minimum depths to start at landscape's min depth and then that + interval
		if (i == 0) { ds.push_back(z_min); }
		else(ds.push_back(z_min + i * dint));
	}


	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++) {
			for (int k = 0; k < nlayers; k++) {
				if (raster_3d[i][j][k]) {
					//raster_3d[i][j][k]->m_x = xes[j];
					//raster_3d[i][j][k]->m_y = ys[i];
					raster_3d[i][j][k]->min_depth = ds[k];
					raster_3d[i][j][k]->max_depth = ds[k] + dint;
					

					//update four_corners information
					//because these are coordinates and not row,col,layer numbers, x goes first!
					raster_3d[i][j][k]->cell_corners.llc.push_back(xes[j]); //low left corner
					raster_3d[i][j][k]->cell_corners.lrc.push_back(xes[j] + resolution); //low right corner
					raster_3d[i][j][k]->cell_corners.ulc.push_back(xes[j]);
					raster_3d[i][j][k]->cell_corners.urc.push_back(xes[j] + resolution);

					raster_3d[i][j][k]->cell_corners.llc.push_back(ys[i] - resolution);
					raster_3d[i][j][k]->cell_corners.lrc.push_back(ys[i] - resolution);
					raster_3d[i][j][k]->cell_corners.ulc.push_back(ys[i]);
					raster_3d[i][j][k]->cell_corners.urc.push_back(ys[i]);
					raster_3d[i][j][k]->calc_midpoint();
				}
			}

		}
	}
	
}


void Landscape::set_NoData(int layer, string input_loc, string filename) { //this is only used on habitat cells! matrix cells are always going to be separated from a nodata cell by a habitat cell
	string fpar = filename;
	ifstream hfile(input_loc + fpar);
	int h; bool create;
	for (int r = 0; r < nrows; r++) {

		for (int c = 0; c < ncols; c++) {

			hfile >> h;

			if (h == no_data) { cout << "hab was nodata"; create = false; }
			else if (layer == 0) { create = true; } //if it's the top layer, create all cells because indivs can settle everywhere
			else if (r == 0 || c == 0) { //if it's the first row or column of the layer,
				if (raster_3d[r][c][layer - 1] == nullptr || raster_3d[r][c][layer - 1]->habtype != 0) { //if the layer above is nullptr, dont need to create, and if it's solid habitat, don't create because it can be update later if next cell is water
					create = false;
				}
				else {
					create = true;
				}
			}
			else { //not in row 0 or col 0, not in layer 0 and not nodata
				if (h == 0) { //if it's water
					create = true; //create the cell
					if (raster_3d[r][c - 1][layer] == nullptr) { //check western neighbour, if it's ND but current cell is water, then create that cell! check the it wasnt actually no data from input
						raster_3d[r][c - 1][layer] = new Cell(mod); //create the cell object
						raster_3d[r][c - 1][layer]->habtype = -9;
						raster_3d[r][c - 1][layer]->cell_ID = -9;

					}
					if (raster_3d[r - 1][c][layer] == nullptr) { //check northern neighbour, if it's ND but current cell is water, then create that cell!
						raster_3d[r - 1][c][layer] = new Cell(mod); //create the cell object
						raster_3d[r - 1][c][layer]->cell_ID = -9;
						raster_3d[r - 1][c][layer]->habtype = -9;
					}
				}
				else { //if it's habitat
					if (raster_3d[r - 1][c][layer] != nullptr && raster_3d[r - 1][c][layer]->habtype == 0) { //so if western neighbour is water
						create = true;
					}
					if (raster_3d[r][c - 1][layer] != nullptr && raster_3d[r][c - 1][layer]->habtype == 0) { //northern neighbour
						create = true;
					}
					else {
						if (raster_3d[r][c][layer - 1] == nullptr) { create = false; } //if the cell above it is a nullptr, don't create
						else if (raster_3d[r][c][layer - 1]->habtype == 0) { create = true; } //if its upper neighbour is water, then still create the cell because indivs can settle on top of it 
						else { create = false; } //if it's not got water on any of these three edges, don't create
					}
				}
			}

			if (create == true) {
				raster_3d[r][c][layer] = new Cell(mod); //create the cell object
				raster_3d[r][c][layer]->cell_ID = r * ncols*(layer + 1) + c; //give it an ID number
				raster_3d[r][c][layer]->habtype = h;
				raster_3d[r][c][layer]->K = car_caps[h];

			}
			else {
				//cout << "trying to make nullptr";
				raster_3d[r][c][layer] = nullptr;
			}
		}
	}
	hfile.close(); //need to close it and then re-open it to start back at the beginning
	ifstream hfile2(input_loc + fpar);

	//now need to go through a fill in the information for the new cells
	for (int r = 0; r < nrows; r++) {
		for (int c = 0; c < ncols; c++) {

			hfile2 >> h;
			if (h == no_data) {
				delete(raster_3d[r][c][layer]); //delete the object created with new()
				raster_3d[r][c][layer] = nullptr;//if it's actually a nodata cell, turn it back to a nullptr
				continue;
			}
			else if (raster_3d[r][c][layer] != nullptr && raster_3d[r][c][layer]->cell_ID == -9) {
				raster_3d[r][c][layer]->cell_ID = r * ncols*(layer + 1) + c; //give it an ID number
				raster_3d[r][c][layer]->habtype = h;
				raster_3d[r][c][layer]->K = car_caps[h];
				//cout << "changing nullptr to new cell at " << r << ", " << c << ", " << layer << "with habtype" << h<< endl;
			}

			//do some checks in case i missed something
			if (h > 0 && layer > 0) { //if it's some kind of habitat (and not layer 0 because everything should have been created in layer 0
				if (raster_3d[r][c][layer] != nullptr) { //if the cell has been created,
					//check whether it has any water neighbours at all
					if (c - 1 >= 0) { //western
						if (raster_3d[r][c - 1][layer] != nullptr && raster_3d[r][c - 1][layer]->habtype == 0) { //if it's water
							continue; //that means at least one neighbour is water, so it's fine that the neighbour has been created 
						}
					}
					if (c + 1 < get_land_att(2)) {
						if (raster_3d[r][c + 1][layer] != nullptr && raster_3d[r][c + 1][layer]->habtype == 0) {
							continue;
						}
					}
					if (r - 1 >= 0) {
						if (raster_3d[r - 1][c][layer] != nullptr && raster_3d[r - 1][c][layer]->habtype == 0) {
							continue;
						}
					}
					if (r + 1 < get_land_att(1)) {
						if (raster_3d[r + 1][c][layer] != nullptr && raster_3d[r + 1][c][layer] == 0) {
							continue;
						}
					}
					if (raster_3d[r][c][layer - 1] != nullptr && raster_3d[r][c][layer - 1]->habtype == 0) {
						continue;
					}
					else { //if it has no water neighbours or all surrounding cells are null
						//cout << "a cell was created at " << r << ", " << c << ", " << layer << " that shouldn't have been. making it null" << endl;
						delete(raster_3d[r][c][layer]); //delete the object created with new()
						raster_3d[r][c][layer] = nullptr; //make it nullptr again, it should not have been created
					}

				}
			}
		}
	}
	hfile2.close();
}

//void Landscape::set_NoData(int layer, string input_loc, string filename) { //this is only used on habitat cells! matrix cells are always going to be separated from a nodata cell by a habitat cell
//	string fpar = filename;
//	ifstream hfile(input_loc + fpar);
//	int h; bool create;
//	for (int r = 0; r < nrows; r++) {
//		for (int c = 0;  c < ncols; c++) {
//			
//			hfile >> h;
//			if (layer == 0) { create = true; } //if it's the top layer, create all cells because indivs can settle everywhere
//			else if (h == no_data) { create = false; } //if there is a value here that is no data, dont create the cell
//			else if (r == 0 || c == 0) { create = true; } //if it's the first row or column of the layer, create cell
//			else { //not in row 0 or col 0
//				if (h == 0) { //if it's water
//					create = true; //create the cell
//					if (raster_3d[r][c - 1][layer] == nullptr && h!=no_data) { //check western neighbour, if it's ND but current cell is water, then create that cell! check the it wasnt actually no data from input
//						raster_3d[r][c - 1][layer] = new Cell(mod); //create the cell object
//
//						raster_3d[r][c - 1][layer]->cell_ID = -9;
//
//					}
//					else if (raster_3d[r-1][c][layer] == nullptr && h != no_data) { //check northern neighbour, if it's ND but current cell is water, then create that cell!
//						raster_3d[r - 1][c][layer] = new Cell(mod); //create the cell object
//						raster_3d[r - 1][c][layer]->cell_ID = -9;
//					}
//				} 
//				else { //if it's habitat
//					if (raster_3d[r - 1][c][layer] && raster_3d[r - 1][c][layer]->habtype == 0 || raster_3d[r][c - 1][layer] && raster_3d[r][c - 1][layer]->habtype == 0) { //so if western, or northern neighbour is water
//						create = true;
//					}
//					else {
//						if (layer != 0 && raster_3d[r][c][layer - 1] == nullptr) { create = false; } //if the cell above it is a nullptr, don't create
//						else if (layer!=0 && raster_3d[r][c][layer - 1]->habtype == 0) { create = true; } //if its upper neighbour is water, then still create the cell because indivs can settle on top of it 
//						else { create = false; } //if it's not got water on any of these three edges, don't create
//					}
//				}
//			}
//
//			if (create == true) {
//				raster_3d[r][c][layer] = new Cell (mod); //create the cell object
//				raster_3d[r][c][layer]->cell_ID = r * ncols*(layer + 1) + c; //give it an ID number
//				raster_3d[r][c][layer]->habtype = h;
//				raster_3d[r][c][layer]->K = car_caps[h];
//
//			}
//			else {
//				raster_3d[r][c][layer] = nullptr;
//			}
//		}
//	}
//	hfile.close(); //need to close it and then re-open it to start back at the beginning
//	ifstream hfile2(input_loc + fpar);
//
//	//now need to go through a fill in the information for the new cells
//	for (int r = 0; r < nrows; r++) {
//		for (int c = 0; c < ncols; c++) {
//
//			hfile2 >> h;
//			if (raster_3d[r][c][layer] != nullptr && raster_3d[r][c][layer]->cell_ID == -9) {
//				raster_3d[r][c][layer]->cell_ID = r * ncols*(layer + 1) + c; //give it an ID number
//				raster_3d[r][c][layer]->habtype = h;
//				raster_3d[r][c][layer]->K = car_caps[h];
//				//cout << "changing nullptr to new cell at " << r << ", " << c << ", " << layer << "with habtype" << h<< endl;
//
//			}
//		}
//	}
//	hfile2.close();
//}

void Landscape::fill_cells(string input_loc, vector<string> files_to_fill) {
	create_raster();
	ifstream ulayers(input_loc + files_to_fill[0]), vlayers(input_loc + files_to_fill[1]), wlayers(input_loc + files_to_fill[2]), hlayers(input_loc + files_to_fill[3]), 
		players(input_loc + files_to_fill[4]), splayers(input_loc + files_to_fill[5]), costlayers(input_loc+files_to_fill[6]), templayers(input_loc+files_to_fill[7]);
	

	string fpar;
	float upar, vpar, wpar, tpar, dpar; //u v and w are floats; temp is also a float, as i absolute depth of seafloor
	int hpar, ppar, spar, cpar; //habitat type, patch number, spdist and cost are int
	for (int j = 0; j < nlayers; j++) {
		ulayers >> fpar; ifstream ufile(input_loc + fpar);
		vlayers >> fpar; ifstream vfile(input_loc + fpar);
		wlayers >> fpar; ifstream wfile(input_loc + fpar);
		players >> fpar; ifstream pfile(input_loc + fpar); 
		splayers >> fpar; ifstream spfile(input_loc + fpar);
		costlayers >> fpar; ifstream cfile(input_loc + fpar);
		templayers >> fpar; ifstream tfile(input_loc + fpar);
		ifstream dfile(input_loc + files_to_fill[8]); //there is only one layer of this but i need to read it every time

		hlayers >> fpar;
		//cout << "calling set_NoData" << endl;
		set_NoData(j, input_loc, fpar); //this assigns habitat type and decides NoData cells. from here on, we can just check for nullptr
		//cout << "set no data for file " << endl;


		for (int k = 0; k < nrows; k++) {

			for (int m = 0; m < ncols; m++) {

				ufile >> upar; vfile >> vpar; wfile >> wpar; dfile >> dpar; if (patch_based == true) { pfile >> ppar; } if (sp_dist == true) { spfile >> spar; };
				if (files_to_fill[6] != "NULL") { cfile >> cpar; }
				if (files_to_fill[7] != "NULL") { tfile >> tpar; }

				if (raster_3d[k][m][j]) { //check if it's valid, if yes:
					
					raster_3d[k][m][j]->c = m; //give it its col
					raster_3d[k][m][j]->r = k; //row
					raster_3d[k][m][j]->l = j; //and layer coordinates
					raster_3d[k][m][j]->res = resolution;
					raster_3d[k][m][j]->seafloor_depth = dpar;
					//if (files_to_fill[6] != "NULL") { raster_3d[k][m][j]->cost = cpar; }
					if (files_to_fill[7] != "NULL") { raster_3d[k][m][j]->cell_temp = tpar; }
					
					//current speeds and angle calculations
					//if (upar != no_data && vpar != no_data && wpar != no_data) {
					if (raster_3d[k][m][j]->habtype == 0) {

						raster_3d[k][m][j]->u = upar * 60 * 60; //these velocities are in m/s so make them m/hr
						raster_3d[k][m][j]->v = vpar * 60 * 60;
						raster_3d[k][m][j]->w = wpar * 60 * 60;
						raster_3d[k][m][j]->cell_angle = atan2(raster_3d[k][m][j]->v, raster_3d[k][m][j]->u);
						float sum;
						if (wpar != -9) { sum = pow(raster_3d[k][m][j]->v, 2) + pow(raster_3d[k][m][j]->u, 2) + pow(raster_3d[k][m][j]->w, 2); } //if there is vertical velocity, 3D movement
						else { sum = pow(raster_3d[k][m][j]->v, 2) + pow(raster_3d[k][m][j]->u, 2); } //otherwise 2D movement
						raster_3d[k][m][j]->cell_speed = sqrt(sum);
						//cout << " cell " << k << ", " << m << ", " << j <<" has speed " << raster_3d[k][m][j]->cell_speed << "from u " << raster_3d[k][m][j]->u << ", v " << raster_3d[k][m][j]->v << ", w " << raster_3d[k][m][j]->w << endl;
					}

					else { //either its a nodata cell or it's habitat, which doesnt have these parameters
						raster_3d[k][m][j]->u = -9999;
						raster_3d[k][m][j]->v = -9999;
						raster_3d[k][m][j]->w = -9999;
						raster_3d[k][m][j]->cell_angle = -9999;
						raster_3d[k][m][j]->cell_speed = -9999;
					}
					
					//cell area
					
					if ( raster_3d[k][m][j] && raster_3d[k][m][j]->K>0) { //if it is suitable habitat and not a NODATA cell

						//calculate area based on neighbours
						raster_3d[k][m][j]->ha = (4 * (dint*resolution) + pow(resolution, 2)) / 10000; //maximum area possible (can't settle underneath a cell. no floating habitat
						if (k == 0) { //if it's the first row, i won't check previous rows
							raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //but still subtract one side because it's inaccessible
							if (m == 0) { //if it's the first column, i won't check previous columns
								raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //but still subtract one side because it's inaccessible
								if (j > 0) { //if not the top layer
									if (raster_3d[k][m][j - 1] && raster_3d[k][m][j - 1]->habtype > 0) {//cell in layer above is valid and not water
										raster_3d[k][m][j]->ha -= pow(resolution, 2) / 10000; // layers will cover the x,y plane, so resolution^2
									}
								}
							}
							else { //not the first column
								if (raster_3d[k][m - 1][j] && raster_3d[k][m - 1][j]->habtype > 0) { //cell in previous column valid and not water
									raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //depth interval * resolution, divide by 10,000 to make it into hectares, subtract from total possible area
									raster_3d[k][m - 1][j]->ha -= (dint*resolution) / 10000; //also subtract from previous cell
									raster_3d[k][m - 1][j]->calc_maxinds(); //update previous cell's max inds with it's hab type K
								} 
								if (j > 0) { //not the top layer
									if (raster_3d[k][m][j - 1] && raster_3d[k][m][j - 1]->habtype > 0) {//cell in layer above is valid and not water
										raster_3d[k][m][j]->ha -= pow(resolution, 2) / 10000; // layers will cover the x,y plane, so resolution^2
									}
								}
							}
						}
						
						else { //not first row
							if (raster_3d[k - 1][m][j] && raster_3d[k - 1][m][j]->habtype > 0) {//cell in previous row is valid and not water
								raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //subtract one side
								raster_3d[k - 1][m][j]->ha -= (dint*resolution) / 10000; //subtract from cell in previous row too
								raster_3d[k - 1][m][j]->calc_maxinds(); //update previous cell's max inds with it's hab type K
							}
							if (m == 0) { //if first column
								raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //subtract one side because it's inaccessible
							}
							else { //not first column
								if (raster_3d[k][m - 1][j] && raster_3d[k][m - 1][j]->habtype > 0) { //cell in previous column valid and not water
									raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //depth interval * resolution, divide by 10,000 to make it into hectares, subtract from total possible area
									raster_3d[k][m - 1][j]->ha -= (dint*resolution) / 10000; //subtract from cell in previous column too
									raster_3d[k][m - 1][j]->calc_maxinds(); //update previous cell's max inds with it's hab type K
								}
							}
							if (j > 0) { //not top later
								if (raster_3d[k][m][j - 1] && raster_3d[k][m][j - 1]->habtype > 0) {//cell in layer above is valid and not water
									raster_3d[k][m][j]->ha -= pow(resolution, 2) / 10000; // layers will cover the x,y plane, so resolution^2
								}
							}
						}
						suitable_cells.push_back(raster_3d[k][m][j]);  //if the habitat type has a carrying capacity, add it to the vector of suitable cells
						raster_3d[k][m][j]->calc_maxinds();
						//don't have to do layer because there is no floating habitat! but overhangs arent handled...
					}
				}
				
				//patch numbers
				if (patch_based == true) {
					if (raster_3d[k][m][j] == nullptr) { continue; } //this would happen if the cell was in the middle of a patch, with cells above it in upper layers. this is a cell that an indiv will never reach!
					else if (ppar != -9999 ) {
						raster_3d[k][m][j]->patch = patch_vector[ppar];
						patch_vector[ppar]->patch_info.size_cells++;//everytime you fill a cell with this patch number, increase the size of the patch (unless it's 0, which is matrix)
						patch_vector[ppar]->included_cells.push_back(raster_3d[k][m][j]); //add this cell's address to the patch's vector of included cells
						patch_vector[ppar]->set_limits(k, m, j, dint); //check to see if the patch's x,y,z limits need to be updated

					}
					else {
						raster_3d[k][m][j]->patch = nullptr;
					}
					

				}

				//species distribution
				if (sp_dist == true) { //species distribution
					if (raster_3d[k][m][j] == nullptr) { continue; }
					raster_3d[k][m][j]->spdist = spar;
					if (spar == 1) { //if this cell is part of spdist
						spdist_cells.push_back(raster_3d[k][m][j]); //add it to the species distribution vector
						if (patch_based == true) {
							raster_3d[k][m][j]->patch->in_spdist = true; // patch is now "in species distribution"
							raster_3d[k][m][j]->patch->in_spdist_cells.push_back(raster_3d[k][m][j]); //add cell to spdist vector in patch

						}
					}
				}
				
			}
			
		}
		
		ufile.close(); vfile.close(); wfile.close(); dfile.close(); if (patch_based == true) { pfile.close(); } if (sp_dist == true) { spfile.close(); }; cfile.close(); tfile.close();
		cout << "included data for layer " << j << endl;
	}

	ulayers.close(); vlayers.close(); wlayers.close(); hlayers.close(); players.close(); splayers.close(); costlayers.close(); templayers.close();

	cout << "read in info, putting md" << endl;
	fill_md();//will fill in metres and depths. this also checks for nullptr
	if (patch_based) { //calculate K and max_inds per patch
		for (int p = 0; p < npatches; p++) {
			patch_vector[p]->calc_K("patch_info"); //need to tell it to calculate K on scale of whole patch
		}
	}
	cout << "there are " << suitable_cells.size() << " suitable cells in this landscape" << endl;
}


//old version of fill_cells()

//void Landscape::fill_cells(string input_loc, vector<string> files_to_fill) {
//	create_raster();
//
//	for (int i = 0; i < 6; i++) { //0:u, 1:v, 2: w, 3:habtype, 4:patches, 5:spdist
//		ifstream parfile(input_loc + files_to_fill[i]); //open the file containing the names of each layer's files
//		string p;
//		float s;
//		int x;
//		int cell_counter = 0;
//
//		for (int j = 0; j < nlayers; j++) {
//			parfile >> p; //read each depth's file name
//			ifstream pardepth(input_loc + p); //open the file
//
//			for (int k = 0; k < nrows; k++) {
//				
//				for (int m = 0; m < ncols; m++) {
//					
//					if (i == 0) {//first thing we need is to see if the cell is no data. THIS ASSUMES NO velocity measure=NODATA!
//
//						pardepth >> s; //read the file
//						
//						if (s != no_data) { //if it's got data in it
//							raster_3d[k][m][j] = new Cell; //create the cell object
//							raster_3d[k][m][j]->cell_ID = cell_counter; //give it an ID number
//							raster_3d[k][m][j]->u = s; //give it the easterly flow velocity
//							raster_3d[k][m][j]->c = m; //give it its col
//							raster_3d[k][m][j]->r = k; //row
//							raster_3d[k][m][j]->l = j; //and layer coordinates
//							raster_3d[k][m][j]->res = resolution;
//							cell_counter++; //only increment for valid cells
//						}
//						else { //if it is no data
//							cout << "cell" << k << "," << m << "," << j << "is null because s is " << s << endl;
//							raster_3d[k][m][j] = nullptr; //set it to a null pointer
//						}
//
//					}
//					else if (i == 1) { //v
//						pardepth >> s;
//						if (raster_3d[k][m][j]) {
//							raster_3d[k][m][j]->v = s;
//						}
//					}
//					else if (i == 2) { //w
//						pardepth >> s;
//						if (raster_3d[k][m][j]) {
//							raster_3d[k][m][j]->w = s;
//							raster_3d[k][m][j]->cell_angle = atan2(raster_3d[k][m][j]->v, raster_3d[k][m][j]->u);
//							float sum;
//							if (s != -9) { sum = pow(raster_3d[k][m][j]->v, 2) + pow(raster_3d[k][m][j]->u, 2) + pow(raster_3d[k][m][j]->w, 2); } //if there is vertical velocity, 3D movement
//							else { sum = pow(raster_3d[k][m][j]->v, 2) + pow(raster_3d[k][m][j]->u, 2); } //otherwise 2D movement
//							raster_3d[k][m][j]->cell_speed = sqrt(sum);
//						}
//					}
//					
//					else if (i == 3) { //habitat type
//						
//						pardepth >> x;
//						
//						if (raster_3d[k][m][j]) {
//							raster_3d[k][m][j]->habtype = x;
//							raster_3d[k][m][j]->K = car_caps[x];
//
//
//							if (car_caps[x] > 0) { //if it is suitable habitat
//								
//								//calculate area based on neighbours
//								raster_3d[k][m][j]->ha = (4 * (dint*resolution) + pow(resolution, 2)) / 10000; //maximum area possible (can't settle underneath a cell. no floating habitat
//
//								if (k == 0) { //if it's the first row, i won't check previous rows
//									
//									raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //but still subtract one side because it's inaccessible
//									if (m == 0) { //if it's the first column, i won't check previous columns
//										
//										raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //but still subtract one side because it's inaccessible
//										if (j > 0) { //if not the top layer
//											
//											if (raster_3d[k][m][j-1] && raster_3d[k][m][j-1]->habtype > 0) {//cell in layer above is valid and not water
//												raster_3d[k][m][j]->ha -= pow(resolution, 2) / 10000; // layers will cover the x,y plane, so resolution^2
//											}
//										}
//
//									}
//									else { //not the first column
//										if (raster_3d[k][m-1][j] && raster_3d[k][m-1][j]->habtype > 0) { //cell in previous column valid and not water
//											raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //depth interval * resolution, divide by 10,000 to make it into hectares, subtract from total possible area
//											raster_3d[k][m-1][j]->ha -= (dint*resolution) / 10000; //also subtract from previous cell
//											raster_3d[k][m-1][j]->calc_maxinds(); //update previous cell's max inds with it's hab type K
//										}
//										if (j > 0) { //not the top layer
//											if (raster_3d[k][m][j-1] && raster_3d[k][m][j-1]->habtype > 0) {//cell in layer above is valid and not water
//												raster_3d[k][m][j]->ha -= pow(resolution, 2) / 10000; // layers will cover the x,y plane, so resolution^2
//											}
//										}
//									}
//
//								}
//								else { //not first row
//									if (raster_3d[k-1][m][j] && raster_3d[k-1][m][j]->habtype > 0) {//cell in previous row is valid and not water
//										raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //subtract one side
//										raster_3d[k-1][m][j]->ha -= (dint*resolution) / 10000; //subtract from cell in previous row too
//										raster_3d[k-1][m][j]->calc_maxinds(); //update previous cell's max inds with it's hab type K
//
//									}
//									if (m == 0) { //if first column
//										raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //subtract one side because it's inaccessible
//									}
//									else { //not first column
//										if (raster_3d[k][m-1][j] && raster_3d[k][m-1][j]->habtype > 0) { //cell in previous column valid and not water
//											raster_3d[k][m][j]->ha -= (dint*resolution) / 10000; //depth interval * resolution, divide by 10,000 to make it into hectares, subtract from total possible area
//											raster_3d[k][m-1][j]->ha -= (dint*resolution) / 10000; //subtract from cell in previous column too
//											raster_3d[k][m-1][j]->calc_maxinds(); //update previous cell's max inds with it's hab type K
//
//										}
//									}
//									if (j > 0) { //not top later
//										if (raster_3d[k][m][j-1] && raster_3d[k][m][j-1]->habtype > 0) {//cell in layer above is valid and not water
//											raster_3d[k][m][j]->ha -= pow(resolution, 2) / 10000; // layers will cover the x,y plane, so resolution^2
//										}
//									}
//								}
//
//								suitable_cells.push_back(raster_3d[k][m][j]);  //if the habitat type has a carrying capacity, add it to the vector of suitable cells
//								raster_3d[k][m][j]->calc_maxinds();
//
//								//don't have to do layer because there is no floating habitat! but overhangs arent handled...
//							}
//						}
//					}
//
//					else if (i == 4 && patch_based == true) { //patch number
//					pardepth >> x;
//					if (raster_3d[k][m][j]) {
//						raster_3d[k][m][j]->patch = patch_vector[x];
//						if (x != 0) {
//							patch_vector[x]->patch_info.size_cells++;//everytime you fill a cell with this patch number, increase the size of the patch (unless it's 0, which is matrix)
//							patch_vector[x]->included_cells.push_back(raster_3d[k][m][j]); //add this cell's address to the patch's vector of included cells
//							patch_vector[x]->set_limits(k, m, j); //check to see if the patch's x,y,z limits need to be updated
//						}
//					}
//					}
//					else if (i == 5 && sp_dist==true) { //species distribution
//						pardepth >> x;
//						if (raster_3d[k][m][j]) {
//							raster_3d[k][m][j]->spdist = x;
//							if (x == 1) { //if this cell is part of spdist
//								spdist_cells.push_back(raster_3d[k][m][j]); //add it to the species distribution vector
//								raster_3d[k][m][j]->patch->in_spdist = true; // patch is now "in species distribution"
//								raster_3d[k][m][j]->patch->in_spdist_cells.push_back(raster_3d[k][m][j]); //add cell to spdist vector in patch
//							}
//						}
//					}
//				}
//				
//			}
//
//
//
//			//the rest can be done using set_value_3 because that checks for whether it's still a nullptr or not
//			//else if (i == 1) { set_value(pardepth, j, 2); } //v
//			//else if (i == 2) { set_value(pardepth, j, 3); } //w
//			//else if (i == 3) { set_value(pardepth, j, 5); } //habtype
//			//else if (i == 4 && patch_based == true) { set_value(pardepth, j, 4); } //patches
//			//else if (i == 5 && sp_dist == true) { set_value(pardepth, j, 6); } //species distribution
//			
//			pardepth.close();
//		}
//
//		parfile.close();
//
//	}
//
//	cout << "read in info, putting md" << endl;
//	fill_md();//will fill in metres and depths. this also checks for nullptr
//	if (patch_based) { //calculate K and max_inds per patch
//		for (int p = 0; p < npatches; p++) {
//			patch_vector[p]->calc_K("patch_info"); //need to tell it to calculate K on scale of whole patch
//		}
//	}
//}

Landscape::dynland Landscape::declare_dynland(string u, string v, string w, string h, string p) {
	dynland c;
	return c;
}

void Landscape::read_dyn() {
	ifstream dynland_file(mod->wrkdir + "/Inputs/" + mod->dynfile);
	string header;
	int nheaders = 8;
	int sim;


	for (int h = 0; h < nheaders; h++) { dynland_file >> header; }
	dynland_file >> sim >> dynamic_info.num >> dynamic_info.dyn_files >> dynamic_info.burn_in >> dynamic_info.dyn_interval >> dynamic_info.dyn_hab >>
		dynamic_info.dyn_hydro >> dynamic_info.dyn_patches;

	dynland_file.close();

	for (int n = 0; n < dynamic_info.num; n++) {
		dynland c;
		dynamic_info.dyn_land_changes.push_back(c); //now there are enough structs initialised for all the changes
	}
	
	//now read in the actual data and create the structs for each new input set
	nheaders = 6;
	string h_file, u_file, v_file, w_file, p_file;
	ifstream dynfiles(mod->wrkdir + "/Inputs/" + dynamic_info.dyn_files);
	for (int h = 0; h < nheaders; h++) { dynfiles >> header; }
	dynfiles >> sim >> h_file >> p_file >> u_file >> v_file >> w_file;

	if (dynamic_info.dyn_hab == 1) {
		//save the habs file name
	}
	else { dynfiles >> header; }
	if (dynamic_info.dyn_hydro == 1) {
		ifstream ufile(mod->wrkdir + "/Inputs/" + u_file);
		for (int h = 0; h < dynamic_info.num; h++) {
			for (int d = 0; d < nlayers; d++) {
				string filename;
				ufile >> filename;
				dynamic_info.dyn_land_changes[h].u_files.push_back(filename);
			}
		}
		ifstream vfile(mod->wrkdir + "/Inputs/" + v_file);
		for (int h = 0; h < dynamic_info.num; h++) {
			for (int d = 0; d < nlayers; d++) {
				string filename;
				vfile >> filename;
				dynamic_info.dyn_land_changes[h].v_files.push_back(filename);
			}
		}
		ifstream wfile(mod->wrkdir + "/Inputs/" + w_file);
		for (int h = 0; h < dynamic_info.num; h++) {
			for (int d = 0; d < nlayers; d++) {
				string filename;
				wfile >> filename;
				dynamic_info.dyn_land_changes[h].w_files.push_back(filename);
			}
		}
		
	}
	else {
		dynfiles >> header;
	}
	if (dynamic_info.dyn_patches == 1) {
		//save the patchfile name
	}
	else { dynfiles >> header; }

	dynfiles.close();
}

void Landscape::do_dyn(int change_num) { //right now this is only convenient for hydro, i need to figure out an elegant way to potentially do all three dynamic types without repeating loads!
	cout << "doing dynamic seascape number " << change_num << endl;
	for (int lay = 0; lay < nlayers; lay++) {
		ifstream ufile(mod->wrkdir + "/Inputs/" + dynamic_info.dyn_land_changes[change_num].u_files[lay]); //accesses the struct that is the right change number and then the element in the vector that is the right layer
		ifstream vfile(mod->wrkdir + "/Inputs/" + dynamic_info.dyn_land_changes[change_num].v_files[lay]); 
		ifstream wfile(mod->wrkdir + "/Inputs/" + dynamic_info.dyn_land_changes[change_num].w_files[lay]);
		
		float north, east, up;
		for (int r = 0; r < nrows; r++) {
			for (int c = 0; c < ncols; c++) {

				ufile >> east; vfile >> north; wfile >> up;
				if (raster_3d[r][c][lay] == nullptr) { continue; } //if it's null then it doesnt need new hydro data
				if (raster_3d[r][c][lay]->habtype != 0) { continue; } //if it's not water, it doesn't need new hydro data
				raster_3d[r][c][lay]->u = east * 60 * 60; //to bring m/s to m/hr
				raster_3d[r][c][lay]->v = north * 60 * 60;
				raster_3d[r][c][lay]->w = up * 60 * 60;
			}
		}
	}

	//cout << "checking: u at 0,45,0" << raster_3d[0][45][0]->u << endl;
}

//read files and set parameters. calls fill_cells which calls fill_md and set_value
void Landscape::get_info() { 
	if (mod->speciesdist == 0) { sp_dist = false; }

	//read in the LandFile
	ifstream landfile(mod->wrkdir + "/Inputs/" + mod->landfile);
	string header;

	
	int nheaders = 20 + mod->nhabs; //when i start using simulations then this would be nheaders=(9 + mod->nhabs) +(sim_number*9)
	
	for (int h = 0; h < nheaders; h++) {landfile >> header;}

	int sim; //not doing anything with this yet!

	landfile >> sim >> nhabitats >> npatches;
	
	if (mod->patchmodel == 1) {
		patch_based = true;
		create_patch_vector();
	}
	else { cell_based = true; }
	
	for (int k = 0; k < nhabitats; k++) {
		float temp;
		landfile >> temp;
		car_caps.push_back(temp);
	}

	string dimfile;
	landfile >> dimfile;

	//read in the dimensions file
	ifstream dimensionfile(mod->wrkdir + "/Inputs/" + dimfile);
	dimensionfile >> header >> ncols >> header >> nrows >> header >> xllcorner >> header >> yllcorner >> header >> resolution >>
		header >> no_data >> header >> nlayers >> header >> dint >> header >> z_min;

	ncells = nrows * ncols;
	x_min = ceil(xllcorner); //to start off with a whole number
	x_max = x_min + resolution * ncols;
	y_min = floor(yllcorner);
	y_max = y_min + resolution * nrows;
	z_max = z_min + dint * (nlayers);


	dimensionfile.close();

	string us, vs, ws, landscape, abs_Depth, patches, spdist, costs, temps;
	landfile >> us >> vs >> ws >> landscape >> abs_Depth >> patches >> spdist >> costs >> temps;
	vector<string> files_to_fill = { us, vs, ws, landscape, patches, spdist, costs, temps, abs_Depth };
	//cout << "absdepth file is " << abs_Depth << endl;

	landfile >> patch_extinct >> p_ext_method >> p_ext_p_prop >> p_ext_patch >> p_ext_i_prop >> p_ext_int >> p_ext_burnin;
	cout << "proportion of indivs to go extinct is " << p_ext_i_prop<< endl;
	landfile.close();

	cout << "read landfile" << endl;
	//cout << "nodata value is " <<no_data << endl;
	//fill cells with information
	fill_cells(mod->wrkdir + "/Inputs/", files_to_fill);

	if (mod->dynamic == 1) {
		dynamic = 1;
		//cout << "checking baseline: u at 0,45,0" << raster_3d[0][45][0]->u << endl;
		read_dyn();
	}
}


//void Landscape::dyn_land_change(string lf) {
//	//read in the LandFile
//	ifstream landfile(mod->wrkdir + "/Inputs/" + lf);
//	string header;
//
//	int nheaders = 12 + mod->nhabs; //when i start using simulations then this would be nheaders=(9 + mod->nhabs) +(sim_number*9)
//
//	for (int h = 0; h < nheaders; h++) { landfile >> header; }
//
//	int sim; //not doing anything with this yet!
//
//	landfile >> sim >> nhabitats >> npatches;
//
//	if (mod->patchmodel == 1) {
//		patch_based = true;
//		create_patch_vector();
//	}
//	else { cell_based = true; }
//
//	for (int k = 0; k < nhabitats; k++) {
//		float temp;
//		landfile >> temp;
//		car_caps.push_back(temp);
//	}
//
//	string us, vs, ws, landscape, patches, spdist, costs, temps;
//	landfile >> us >> vs >> ws >> landscape >> patches >> spdist >> costs >> temps;
//	vector<string> files_to_fill = { us, vs, ws, landscape, patches, spdist, costs, temps };
//
//	landfile.close();
//
//	cout << "read landfile" << endl;
//
//	//fill cells with information
//	fill_cells(mod->wrkdir + "/Inputs/", files_to_fill);
//}

//this is the function to test memory use of vectors vs 3D array
//void Landscape::fill_neigh() {
//	
//	for (int r = 0; r < get_land_att(1); r++) {
//		for (int c = 0; c < get_land_att(2); c++) {
//			for (int la = 0; la < get_land_att(5); la++) {
//				if (raster_3d[r][c][la] == nullptr) {
//					continue;
//				}
//				else {
//					for (int nr = -1; nr < 2; nr++) {
//						for (int nc = -1; nc < 2; nc++) {
//							for (int nl = -1; nl < 2; nl++) {
//								int row = r + nr, col = c + nc, lay = la + nl;
//								if (row<0 || row>= get_land_att(1) || col<0 || col>= get_land_att(2) || lay<0 || lay>= get_land_att(5)) { //outside boundaries
//									continue; //since the matrix is already initialised to nullptr, i can just skip it
//								}
//								else if (raster_3d[row][col][lay] == nullptr) {
//									continue;
//								}
//								else {
//									if (raster_3d[row][col][lay]->habtype == 0) { //if the neighbour is water
//										//for vectors
//										//if (lay < la) { //layer above
//										//	raster_3d[r][c][la]->open_water_neigh.up_neigh.push_back(raster_3d[row][col][lay]);
//										//}
//										//else if (lay == la) { //same layer
//										//	raster_3d[r][c][la]->open_water_neigh.forward_neigh.push_back(raster_3d[row][col][lay]);
//										//}
//										//else { //layer below
//										//	raster_3d[r][c][la]->open_water_neigh.down_neigh.push_back(raster_3d[row][col][lay]);
//										//}
//
//										//for 3D array
//										raster_3d[r][c][la]->ow_neigh[nr + 1][nc + 1][nl + 1] = raster_3d[row][col][lay];
//									}
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	
//	
//}

void Landscape::check_landscape() {
	for (int lay = 0; lay < get_land_att(5); lay++) {
		for (int r = 0; r < get_land_att(1); r++) {
			for (int c = 0; c < get_land_att(2); c++) {
				//cout << "row " << r << ", col " << c << ", layer " << lay;
				if (raster_3d[r][c][lay] == nullptr) {
					//cout << "was null ptr.";
					continue;
				}
				else if (raster_3d[r][c][lay]->habtype == 0 && raster_3d[r][c][lay]->cell_speed == -9999) {
					cout << "row " << r << ", col " << c << ", layer " << lay << " is water but has no cell speed!"<< endl;;
				}
			}
			
		}
	}

}