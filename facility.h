#ifndef Facility_H
#define Facility_H

//Facility definition
struct Facility {
	int p_idx;				//Facility center p, index of point
	double fc;				//Facility cost, non-uniform clusterization
	TVector <int> U_idxs;			//Connected clients ui, indexes of points
	bool del;				//Facility marked as deleted
	Facility(const unsigned int& p_idx_, const double& fc_) : p_idx(p_idx_), fc(fc_), del(false) {}
};

#endif