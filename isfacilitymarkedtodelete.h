#ifndef IsFacilityMarkedToDelete_H
#define IsFacilityMarkedToDelete_H

#include <string>
#include "point3d.h"
#include "tvector.h"

//Subset of the input point cloud divided by KD-tree
struct KDPointTile
{
	static int kd_tile_counter;						//Counter of created KD point tile
	int kd_tile_id;								//Id of the KD point tile
	int depth;								//Depth of the split
	int dir;								//Split direction
	const TVector <Point3D>& points;						//Points of the subset

	KDPointTile(const int depth_, const int dir_, const TVector<Point3D>& points_) : depth(depth_), dir(dir_), points(points_), kd_tile_id(kd_tile_counter++) {}

};

#endif
