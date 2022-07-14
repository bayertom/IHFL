#ifndef RegressionPlane_H
#define RegressionPlane_H

#include "tvector.h"
#include "point3d.h"

//Parameters of the regression plane given by P
struct RegressionPlane
{
	double  a,	//Parametric eaquation
		b,	//Parametric eaquation
		c,	//Parametric eaquation
		d,	//Parametric eaquation
		sigma,	//Cummulated distance of points from the plane: S[min]
		abn,    //Maximum ABN, point and its neighbors
		var;	//Variability of normals S[min] / (S[min] + S[mid] + S[max])

	RegressionPlane() : a(1.0), b(0.0), c(0.0), d(0.0), sigma(0.0), abn(0.0), var(0.0) {}

	void computeRegressionPlane(const TVector <Point3D>& P);
};

#endif