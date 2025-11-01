// Description: 3D point with r, g, b components and initial facility cost

// Copyright (c) 2021 - 2023
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.


#ifndef Point3D_H
#define Point3D_H

//3D point
struct Point3D
{
	int id;					//Point ID in the file
	double x, y, z;				//Cartesian coordinates
	double fc;				//Facility cost
	short r, g, b;				//R, G, B components
	bool deleted;				//Flag, if a point needs to be deleted

	Point3D() : id(0), x(0.0), y(0.0), z(0.0), fc(1.0), r(0), g(0), b(0), deleted(false) {}
	Point3D(const int& ide_, const double& x_, const double& y_, const double& z_) : id(ide_), x(x_), y(y_), z(z_), fc(1.0), r(0), g(0), b(0), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_, const double& fc_) : id(id_), x(x_), y(y_), z(z_), fc(fc_), r(0), g(0), b(0), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_, const short& r_, const short& g_, const short& b_) : id(id_), x(x_), y(y_), z(z_), fc(1.0), r(r_), g(g_), b(b_), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_, const double& fc_, const short& r_, const short& g_, const short& b_) : id(id_), x(x_), y(y_), z(z_), fc(fc_), r(r_), g(g_), b(b_), deleted(false) {}

	bool operator == (const Point3D& p) const
	{
		return (this->x == p.x) && (this->y == p.y) && (this->z == p.z);
	}
};

#endif
