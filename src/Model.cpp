#include "classes.h"

Model::Model() {

}

Model::~Model()
{
}

void Model::read_control() {
	cout << "what is the filepath to your inputs?"; //remember this needs to have the / after it because I only add the file names in subsequent functions!
	//ie it needs to be ../tests/test1/
	cin >> wrkdir;

	string header;

	ifstream controlfile(wrkdir + "/Inputs/CONTROL.txt");

	controlfile >> header >> batchnum >> header >> nsims >> header >> nreps >> header >> nyears >> header >>
		patchmodel >> header >> nhabs >> header >> speciesdist >> header >> distres >> header >>
		stagestruct >> header >> stages >> header >> repro >> header >> temp_dep >> header >> size_or_time >> header >> 
		growthfile >> header >> grow_method >> header >> dynamic >> header >> landfile >> header >> dynfile >> header >>
		initfile >> header >> popfile >> header >> stagestructfile >> header >> emigfile >> header >>
		transfile >> header >> setfile >> header >> evofile >> header >> Outpop >> header >> Outind >> header >> Outindmov >> header >>
		Outrange >> header >> outpop_int >> header >> outind_int >> header >> outrange_int;
	cout << "growth method is " << grow_method << endl;
	
	cout << "dynamic? " << dynamic << ", dynfile: " << dynfile << endl;
}

void Model::read_init() {

	ifstream init_file(wrkdir + "/Inputs/" + initfile);
	string header; int nheaders;

	if (stagestruct==1) {nheaders = 17 + (stages - 1);} //when i have simulations this will be +sim_num*nheaders
	else { nheaders = 17; }

	for (int h = 0; h < nheaders; h++) {
		init_file >> header;
	}
	int sim;
	init_file >> sim >> seedtype >> freetype >> freeprop >> minX >> maxX >> minY >> maxY >> minZ >> maxZ >>
		sptype >> spprop >> initdens >> indsha >> whichpatch >> propmales;
	//cout << "propmales is " << propmales << endl;
	cout << "init min x is " << minX << ", min y is " << minY << ", min z is " << minZ << endl;
	if (stagestruct == 1) {
		init_file >> initage; //only applies to stage structured because otherwise you have nonoverlapping generations
		for (int s = 1; s < stages; s++) {
			float props;
			init_file >> props;
			stage_props.push_back(props);
		}
	}
}