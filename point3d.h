#ifndef Point3D_H
#define Point3D_H

//3D point
struct Point3D
{
	int id;					//Point ID
	double x, y, z;				//Cartesian coordinates
	double fc;				//Cost
	short r, g, b;				//R, G, B components
	bool deleted;				//Flag, if a point needs to be deleted

	Point3D() : id(0), x(0.0), y(0.0), z(0.0), fc(1.0), r(0), g(0), b(0), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_) : id(id_), x(x_), y(y_), z(z_), fc(1.0), r(0), g(0), b(0), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_, const double& c_) : id(id_), x(x_), y(y_), z(z_), fc(c_), r(0), g(0), b(0), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_, const short& r_, const short& g_, const short& b_) : id(id_), x(x_), y(y_), z(z_), fc(1.0), r(r_), g(g_), b(b_), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_, const double& fc_, const short& r_, const short& g_, const short& b_) : id(id_), x(x_), y(y_), z(z_), fc(fc_), r(r_), g(g_), b(b_), deleted(false) {}

	bool operator == (const Point3D& p) const
	{
		return this->x == p.x && this->y == p.y && this->z == p.z;
	}
};

#endif
