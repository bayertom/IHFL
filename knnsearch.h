#ifndef KNNSearch_H
#define KNNSearch_H

#include "point3d.h"
#include "tvector.h"
#include "tvector2d.h"

//KNN search using Nano-Flann library
class KNNSearch
{
	private:
		TVector <Point3D> points;

	public:
		KNNSearch(const TVector <Point3D>& points_) : points(points_) {}
		void findAllKNN(const TVector <Point3D>& qpoints, const TVector<int>& K, TVector2D <size_t>& knn_id, TVector2D <float>& knn_dist);
};

#endif