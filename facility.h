#ifndef Facility_H
#define Facility_H

//Facility definition
struct Facility {
	int u_id;				//Facility center, ID of point
	double fc;				//Facility cost, non-uniform clusterization
	TVector <int> Cl_id;			//Connected clients, ID of points
	bool del;				//Facility marked as deleted

	Facility(const unsigned int& u_id_, const double& fc_) : u_id(u_id_), fc(fc_), del(false) {}
};

#endif