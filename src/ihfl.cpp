#define _CRT_SECURE_NO_WARNINGS
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <random>
#include <stack>
#include <filesystem>
#include <nanoflann.hpp>
#include <limits>
#include "utils.h"

#include <Core>
#include <SVD>

#include "libalgo/source/types/TVector.h"
#include "libalgo/source/types/TVector2D.h"

#include "libalgo/source/io2/DXFExport.h"

#include "libalgo/source/exceptions/FileReadException.h"

using namespace std;
using namespace nanoflann;

//Initialize static variables: clusterization parameters
static double l = 1.0;
static double mju = 0.9;

struct Point3D
{
	int id;					//Point ID
	double x, y, z;				//Cartesian coordinates
	short r, g, b;				//R, G, B components
	bool deleted;				//Flag, if a point needs to be deleted

	Point3D() : id(0), x(0), y(0), z(0), r(0), g(0), b(0), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_) : id(id_), x(x_), y(y_), z(z_), r(0), g(0), b(0), deleted(false) {}
	Point3D(const int& id_, const double& x_, const double& y_, const double& z_, const short &r_, const short &g_, const short &b_) : id(id_), x(x_), y(y_), z(z_), r(r_), g(g_), b(b_), deleted(false) {}
};


//Check, whether a point needs to be deleted
bool isPointMarkedToDelete (const Point3D &p)
{
	return p.deleted;
}


struct sortPoints3DByX
{
	bool operator() (const Point3D& p1, const Point3D& p2) const
	{
		return p1.x < p2.x;
	}
};

struct sortPoints3DByY
{
	bool operator() (const Point3D& p1, const Point3D& p2) const
	{
		return p1.y < p2.y;
	}
};

struct sortPoints3DByZ
{
	bool operator() (const Point3D& p1, const Point3D& p2) const
	{
		return p1.z < p2.z;
	}
};

//Pointer to the norm
typedef double (*pfnorm) (const Point3D &, const Point3D &, const TVector <double>&, const TVector <double>&, const double);

//Facility
typedef struct Facility {
	Point3D u;			//Facility center
	double m;			//Radius of the ball (for FFL/FNFL algorithm)
	double fc;			//Facility cost: non-uniform clusterization
	TVector <Point3D> C;		//Connected vertices of the cluster
	bool deleted;			//Facility markes as deleted (for IFL algorithm)

	Facility(const Point3D &u_, const double &m_, const double &fc_) : u(u_), m(m_), fc(fc_), deleted(false) {}
};

//Check, whether a facility needs to be deleted
bool isFacilityMarkedToDelete (const Facility & f)
{
	return f.deleted;
}

//Subset of the input point cloud
typedef struct KDSubset
{
	static int kd_subsets_id_counter;					//Counter of created KD subsets
	int kd_subset_id;							//Id of the KD subset
	int depth;								//Depth of the split
	int dir;								//Split direction
	TVector<Point3D> points;						//Points of the subset

	KDSubset(const int depth_, const int dir_, const TVector<Point3D> &points_) : depth(depth_), dir(dir_), points(points_), kd_subset_id(kd_subsets_id_counter++) {}
};


//Initialize static variables: general parameters
int KDSubset::kd_subsets_id_counter = 0;
static double max_float = std::numeric_limits<float>::max();
static double max_int = std::numeric_limits<int>::max();


//WGS84 => JTSK
void WGS84ToJTSK(const double lat, const double lon, const double h, double &X_JTSK, double &Y_JTSK, double &Z_JTSK)
{
	//Deg => rad
	const double PI = 4.0 * atan(1.0), Ro = 57.29577951;

	//WGS-84
	const double A_WGS = 6378137.0000, B_WGS = 6356752.3142;
	const double E2_WGS = (A_WGS * A_WGS - B_WGS * B_WGS) / (A_WGS * A_WGS);

	//Bessel
	const double A_Bes = 6377397.1550, B_Bes = 6356078.9633;
	const double E2_Bes = (A_Bes * A_Bes - B_Bes * B_Bes) / (A_Bes * A_Bes), E_Bes = sqrt(E2_Bes);

	//Scale, Translation, Rotation
	const double m = -3.5623e-6;
	const double X0 = -570.8285, Y0 = -85.6769, Z0 = -462.8420;
	const double OMX = 4.9984 / 3600 / Ro, OMY = 1.5867 / 3600 / Ro, OMZ = 5.2611 / 3600 / Ro;

	//JTSK
	const double  FI0 = 49.5 / Ro;
	const double U0 = (49.0 + 27.0 / 60 + 35.84625 / 3600) / Ro;
	const double ALFA = 1.000597498372;
	const double LA_FERRO = (17.0 + 40.0 / 60) / Ro;
	const double K = 0.9965924869;
	const double R = 6380703.6105;
	const double UK = (59.0 + 42.0 / 60 + 42.6969 / 3600) / Ro;
	const double VK = (42.0 + 31.0 / 60 + 31.41725 / 3600) / Ro;
	const double Ro0 = 1298039.0046;
	const double S0 = 78.5 / Ro;

	//Point parameters
	const double W = sqrt(1 - E2_WGS * sin(lat / Ro) * sin(lat / Ro));
	const double M = A_WGS * (1 - E2_WGS) / (W * W * W);
	const double N = A_WGS / W;

	//Transformation (B,L,H)WGS => (X,Y,Z)WGS
	const double X_WGS = (N + h) * cos(lat / Ro) * cos(lon / Ro);
	const double Y_WGS = (N + h) * cos(lat / Ro) * sin(lon / Ro);
	const double Z_WGS = (N * (1 - E2_WGS) + h) * sin(lat / Ro);

	//Transformation (X,Y,Z)WGS =>(X,Y,Z)Bes
	const double X_Bes = X0 + (m + 1) * (X_WGS + Y_WGS * OMZ - Z_WGS * OMY);
	const double Y_Bes = Y0 + (m + 1) * (-X_WGS * OMZ + Y_WGS + Z_WGS * OMX);
	const double Z_Bes = Z0 + (m + 1) * (X_WGS * OMY - Y_WGS * OMX + Z_WGS);

	//Transformation (X,Y,Z)Bes => (BLH)Bess
	const double rad = sqrt(X_Bes * X_Bes + Y_Bes * Y_Bes);
	double la = 2 * atan(Y_Bes / (rad + X_Bes));
	const double p = atan((A_Bes * Z_Bes) / (B_Bes * rad));
	const double t = (Z_Bes + E2_Bes * A_Bes * A_Bes / B_Bes * pow(sin(p), 3)) / (rad - E2_Bes * A_Bes * pow(cos(p), 3));
	const double fi = atan(t);
	const double H = sqrt(1 + t * t) * (rad - A_Bes / sqrt(1 + (1 - E2_Bes) * t * t));

	//Transformation (fi,la)Bes => (u,v)sphere  (Gauss conformal projection)
	la = la + LA_FERRO;
	const double ro = 1 / K * pow(tan(fi / 2 + PI / 4) * pow((1 - E_Bes * sin(fi)) / (1 + E_Bes * sin(fi)), E_Bes / 2), ALFA);
	const double u = 2 * atan(ro) - PI / 2;
	const double v = ALFA * la;

	//Transformation (u,v)sphere => (s,d)sphere
	const double s = asin(sin(UK) * sin(u) + cos(UK) * cos(u) * cos(VK - v));
	const double d = asin(sin(VK - v) * cos(u) / cos(s));

	//Transformation (s,d)sphere => (Ro,Eps)plane (Lambert conformal projection)
	const double n = sin(S0);
	const double Ro_JTSK = Ro0 * pow((tan((S0 / 2 + PI / 4)) / (tan(s / 2 + PI / 4))), n);
	const double eps_JTSK = n * d;

	//(Ro, eps) => (x,y)
	X_JTSK = Ro_JTSK * cos(eps_JTSK);
	Y_JTSK = Ro_JTSK * sin(eps_JTSK);
	Z_JTSK = H;
}


TVector <Point3D> WG84PointsToJTSK(const TVector <Point3D> & points)
{
	//Convert point in WGS-84 to JTSK
	TVector <Point3D> points_output;

	for (int i = 0; i < points.size(); i++)
	{
		//Conversion
		double x_jtsk, y_jtsk, z_jtsk;
		WGS84ToJTSK(points[i].y, points[i].x, points[i].z, x_jtsk, y_jtsk, z_jtsk);
		
		//Create converted point
		Point3D p(i, x_jtsk, y_jtsk, z_jtsk, points[i].r, points[i].g, points[i].b);
		
		//Add point to the list
		points_output.push_back(p);
	}

	return points_output;
}


void loadPoints(const std::string& file_name, TVector <Point3D> & nl)
{
	//Load points from the list
	int index = 0;
	std::string line;
	std::ifstream file;

	try
	{
		//Open file
		file.open(file_name);

		//Read line by line
		while (std::getline(file, line))
		{
			char* cline = &line[0u];
			const char* item = strtok(cline, " \t");

			//Delimit the row
			std::vector<std::string> row;
			while (item != NULL)
			{
				row.push_back(item);
				item = strtok(NULL, " \t");
			}

			//Add point P = [x, y, z] to the list
			if (row.size() == 3)
				nl.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2])));

			//Add point P = [x, y, z, r, g, b] to the list
			else if (row.size() >= 6)
				nl.push_back(Point3D(index, std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stoi(row[3]), std::stoi(row[4]), std::stoi(row[5])));

			index++;
		}
	}

	//Throw exception
	catch (std::ifstream::failure e)
	{
		std::cerr << "Can not open file> :" << file_name << '\n';
	}
}


void loadKDSubsets(const std::string &file_name, TVector<std::string> &file_name_subsets)
{
	//Load points from the list
	std::string line;
	std::ifstream file;

	try
	{
		//Open file
		file.open(file_name);

		//Read line by line
		while (std::getline(file, line))
		{
			//Skip empty line
			if (line.size() > 0)
				file_name_subsets.push_back(line);
		}
	}

	//Throw exception
	catch (std::ifstream::failure e)
	{
		std::cerr << "Can not open file> :" << file_name << '\n';
	}
}


void savePoints(const std::string& file_name, const TVector <Point3D> & nl)
{
	//Save points to the file
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		file << std::fixed;

		//Process all points
		for (auto point : nl)
		{
			//Write coordinates
			file << std::setprecision(8);
			file << point.x << "  " << point.y << "  " << point.z << "  ";
			
			//Write r, g, b components
			file << std::setprecision(0);
			file << point.r << "  " << point.g << "  " << point.b << '\n';
		}
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}


void saveKDSubsets(const std::string& file_name, const TVector<std::string> &file_name_subsets)
{
	//Save points to the file
	std::string line;
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name);

		//Process all kd subsets
		for (auto f : file_name_subsets)
			file << f << '\n';
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}


void splitPoints(const TVector <Point3D> & U, const double median, const int dir, TVector <Point3D> & UL, TVector <Point3D> & UU)
{
	//Split set of points in median according to the direction
	for (auto p : U)
	{
		//Split according to X
		if (dir == 0)
		{
			if (p.x <= median)
				UL.push_back(p);
			else
				UU.push_back(p);
		}

		//Split according to Y
		else if (dir == 1)
		{
			if (p.y <= median)
				UL.push_back(p);
			else
				UU.push_back(p);
		}

		//Split according to Z
		else
		{
			if (p.z <= median)
				UL.push_back(p);
			else
				UU.push_back(p);
		}
	}
}


void createKDSubsets(const std::string& file_name, const int n_max, TVector <std::string>& file_names)
{
	//Create subsets from the input point cloud divided in median
	TVector <Point3D> U;

	std::cout << ">> Loading point cloud: ";
	loadPoints(file_name, U);
	std::cout << "OK \n";

	std::cout << ">> Partition to subsets: ";
	KDSubset kds(0, 0, U);

	//Create stack
	std::stack <KDSubset> S;
	S.push(kds);

	//Divide in median until S is empty
	while (!S.empty())
	{
		//Get element on the top of the stack
		KDSubset kdss = S.top();

		//Subset too large, split
		double median;
		if (kdss.points.size() > n_max)
		{
			//Find median in x
			if (kdss.dir == 0)
			{
				std::nth_element(kdss.points.begin(), kdss.points.begin() + kdss.points.size() / 2, kdss.points.end(), sortPoints3DByX());
				median = kdss.points[kdss.points.size() / 2].x;
			}

			//Find median in y
			else if (kdss.dir == 1)
			{
				std::nth_element(kdss.points.begin(), kdss.points.begin() + kdss.points.size() / 2, kdss.points.end(), sortPoints3DByY());
				median = kdss.points[kdss.points.size() / 2].y;
			}

			//Find median in z
			else
			{
				std::nth_element(kdss.points.begin(), kdss.points.begin() + kdss.points.size() / 2, kdss.points.end(), sortPoints3DByZ());
				median = kdss.points[kdss.points.size() / 2].z;
			}

			//Split in two data sets
			TVector <Point3D> UL, UU;
			splitPoints(kdss.points, median, kdss.dir, UL, UU);

			//Change direction
			const int kdss_dir = (kdss.dir + 1) % 3;

			//Increment depth
			const int kdss_depth = kdss.depth + 1;

			//Create 2 KD subsets
			KDSubset kdsl(kdss_depth, kdss_dir, UL);
			KDSubset kdsu(kdss_depth, kdss_dir, UU);

			//Remove kds
			S.pop();

			//Add subsets to the stack
			S.push(kdsl);
			S.push(kdsu);

			std::cout << ".";
		}

		//Point cloud of the acceptable size: store to the file
		else
		{
			//Remove kds
			S.pop();

			//Create file name
			std::string file_name_kdss = file_name + "_" + std::to_string(kdss.kd_subset_id) + "_" + std::to_string(kdss.depth) + "_" + std::to_string(kdss.dir);

			//Save file
			savePoints(file_name_kdss, kdss.points);

			//Add file name to the list
			file_names.push_back(file_name_kdss);
		}
	}

	std::cout << "OK \n\n";
}


void generateClusters(const double w, const double h, const double rad, const int nc, const int n, TVector <Point3D>& U)
{
	//Generate nc clusters, each is represented by n points
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu1(0.0, 1.0);

	//Generate points
	for (int i = 0; i < nc; i++)
	{
		//Generate random center of the cluster
		const double xc = disu1(gen) * w - rad;
		const double yc = disu1(gen) * h - rad;

		//Generate points of the cluster
		std::uniform_real_distribution<double> disu2(xc, xc + rad);
		std::uniform_real_distribution<double> disu3(yc, rad);

		for (int j = 0; j < n; j++)
			U.push_back(Point3D(0, disu2(gen), disu3(gen), 0));
	}
}


void sideConePoint(const double a, const double b, double& h, double& r, double& t)
{
	//Points on the side of cone
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	h = a * sqrt(disu(gen));
	r = (b / a) * h;
	t = 2.0 * M_PI * disu(gen);
}


void baseConePoint(const double a, const double b, double& h, double& r, double& t)
{
	//Points on the base of cone
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	h = a;
	r = b * sqrt(disu(gen));
	t = 2.0 * M_PI * disu(gen);
}


void sideCylinderPoint(const double a, const double b, double& h, double& r, double& t)
{
	//Points on the side of cylinder
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	h = a * disu(gen);
	r = b;
	t = 2.0 * M_PI * disu(gen);
}


void baseCylinderPoint(const double a, const double b, double& h, double& r, double& t)
{
	//Points on the base of cone
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> disi(0, 1);
	std::uniform_real_distribution<double> disu(0, 1);

	h = disi(gen) * a;
	r = b * sqrt(disu(gen));
	t = 2.0 * M_PI * disu(gen);
}


void generateCone(const double a, const double b, const int n, TVector <Point3D>& U)
{
	//Generate points on the cone surface
	double t = 0;

	//Initialize random number generator
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	//Generate points
	for (int i = 0; i < n; i++)
	{
		double h, r, t;

		//Point on the base
		if (disu(gen) < b / (b + sqrt(a * a + b * b)))
			baseConePoint(a, b, h, r, t);

		//Point on the side
		else
			sideConePoint(a, b, h, r, t);

		//Cone coordinates
		const double x = r * cos(t);
		const double y = h;
		const double z = r * sin(t);

		//Apend to the list
		U.push_back(Point3D(0, x, z, -y));
	}
}


void generateCylinder(const double a, const double b, const int n, TVector <Point3D> & U)
{
	//Generate cone
	double t = 0;

	//Initialize random number generator
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	//Generate points
	for (int i = 0; i < n; i++)
	{
		double h, r, t;

		//Point on the base
		if (disu(gen) < b / (2 * (a + b)))
			baseCylinderPoint(a, b, h, r, t);

		//Point on the side
		else
			sideCylinderPoint(a, b, h, r, t);

		//Cone coordinates
		const double x = r * cos(t);
		const double y = h;
		const double z = r * sin(t);

		//Apend to the list
		U.push_back(Point3D(0, x, z, -y));
	}
}


void generateCube(const double a, const int n, TVector <Point3D>& U)
{
	//Generate points on the surface of cube
	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_int_distribution <int> disi(0, 5);
	std::uniform_real_distribution<double> disu(0, 1);

	//Generate points
	for (int i = 0; i < n; i++)
	{
		//Select random side of the cube
		int side = disi(gen);

		//Create empty row
		double xyz[] = { 0, 0, 0 };

		//Coordinate index
		int c = side % 3;

		//Create a point
		if (side > 2)
			xyz[c] = a;
		else
			xyz[c] = 0;

		xyz[(c + 1) % 3] = a * disu(gen);
		xyz[(c + 2) % 3] = a * disu(gen);

		//Apend to the list
		U.push_back(Point3D(0, xyz[0], xyz[1], xyz[2]));
	}
}


double dist(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	//Compute distance
	const double dx = x2 - x1;
	const double dy = y2 - y1;
	const double dz = z2 - z1;

	return sqrt(dx * dx + dy * dy + dz * dz);
}


double getPointPlaneDist(const double qx, const double qy, const double qz, const double a, const double b, const double c, const double d)
{
	//Distance of the point q = [qz, qy, qz] from the plane given by (a, b, c, d)
	return fabs(a * qx + b * qy + c * qz + d) / sqrt(a * a + b * b + c * c);
}


double getPointLineDistance(const Point3D& q, const Point3D& p1, const Point3D& p2)
{
	//Get perpendicular distance between point q and the line (p1, p2)
	
	//Direction vector (p2 - p1)
	const double d12x = p2.x - p1.x;
	const double d12y = p2.y - p1.y;
	const double d12z = p2.z - p1.z;

	//Direction vector (q - p1)
	const double dq1x = q.x - p1.x;
	const double dq1y = q.y - p1.y;
	const double dq1z = q.z - p1.z;

	//Compute the cross product
	const double cx = d12y * dq1z - dq1y * d12z;
	const double cy = d12z * dq1x - dq1z * d12x;
	const double cz = d12x * dq1y - dq1x * d12y;

	//Distance (p1, p2)
	const double d12 = sqrt(d12x * d12x + d12y * d12y + d12z * d12z);

	//Perpendicular distance
	const double numer = sqrt(cx * cx + cy * cy + cz * cz);
	return numer / d12;
}


double getPointLineSegmentDistance(const Point3D& q, const Point3D& p1, const Point3D& p2)
{
	//Get distance between point q and the line segment (p1, p2)

	//Direction vector (p2 - p1)
	const double d12x = p2.x - p1.x;
	const double d12y = p2.y - p1.y;
	const double d12z = p2.z - p1.z;

	//Direction vector (q - p1)
	const double dq1x = q.x - p1.x;
	const double dq1y = q.y - p1.y;
	const double dq1z = q.z - p1.z;

	//Point behind the start point p1 of the segment
	const double dot_s = d12x * dq1x + d12y * dq1y + d12z * dq1z;

	//Return distance q, p1
	if (dot_s <= 0.0)
		return sqrt(dq1x * dq1x + dq1y * dq1y + dq1z * dq1z);

	//Direction vector (q - p2)
	const double dq2x = q.x - p2.x;
	const double dq2y = q.y - p2.y;
	const double dq2z = q.z - p2.z;

	//Point after the end point p2 of the segment
	const double dot_e = d12x * dq2x + d12y * dq2y + d12z * dq2z;

	//Return distance q, p2
	if (dot_e >= 0.0)
		return sqrt(dq2x * dq2x + dq2y * dq2y + dq2z * dq2z);

	//Otherwise, perpendicular distance point and line
	return getPointLineDistance(q, p1, p2);
}


double abn(const TVector <double>& na, const TVector <double>& nb)
{
	//Compute abn criterion
	const double norma = sqrt(na[0] * na[0] + na[1] * na[1] + na[2] * na[2]);
	const double normb = sqrt(nb[0] * nb[0] + nb[1] * nb[1] + nb[2] * nb[2]);
	const double dot_ab = fabs(na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2]);

	return acos(std::min(1.0, dot_ab / (norma * normb)));
}


double ablp(const Point3D& a, const Point3D& b, const TVector <double>& na, const TVector <double>& nb)
{
	//Compute average angle between line (a, b) and plane nb and line (b, a) and plane na
	const double ux = b.x - a.x;
	const double uy = b.y - a.y;
	const double uz = b.z - a.z;

	const double norma = sqrt(na[0] * na[0] + na[1] * na[1] + na[2] * na[2]);
	const double normb = sqrt(nb[0] * nb[0] + nb[1] * nb[1] + nb[2] * nb[2]);
	const double normu = sqrt(ux * ux + uy * uy + uz * uz);

	const double dot_nau = fabs(ux * na[0] + uy * na[1] + uz * na[2]);
	const double dot_nbu = fabs(-ux * nb[0] - uy * nb[1] - uz * nb[2]);

	const double ablp1 = asin(std::min(1.0, dot_nau / (norma * normu)));
	const double ablp2 = asin(std::min(1.0, dot_nbu / (normb * normu)));

	return 0.5 * (ablp1 + ablp2);
}


double dfp(const Point3D& a, const Point3D& b, const TVector <double>& na, const TVector <double>& nb)
{
	//Compute dfp criterion
	const double dfa = getPointPlaneDist(a.x, a.y, a.z, nb[0], nb[1], nb[2], nb[3]);
	const double dfb = getPointPlaneDist(b.x, b.y, b.z, na[0], na[1], na[2], na[3]);

	///std::cout << dfa << " : " << dfb << " : " << sqrt(dfa * dfa + dfb * dfb) << "   ";
	return sqrt(dfa * dfa + dfb * dfb);
}


double nL2(const Point3D& a, const Point3D& b, const TVector <double>& an, const TVector <double>& bn, const double radius)
{
	//L2 norm
	return dist(a.x, a.y, a.z, b.x, b.y, b.z);
}


double nDIS(const Point3D& a, const Point3D& b, const TVector <double>& na, const TVector <double>& nb, const double radius)
{
	//Norm: combined DIS criterion, tangent model
	const double dab = nL2(a, b, na, nb, radius);

	//Identical points
	if (dab < MIN_FLOAT)
		return 0;

	//Compute abn criterion
	const double dabn = abn(na, nb);
	const double dh = sin(0.5 * dabn) * dab;

	//Constrained (hybrid) metric
	const double cmax = sqrt(2.0) / 2.0;
	
	/*
	if (dab < radius)
		return dh;
	else
		return dh + pow(dab - radius, l);
	
	//std::cout << radius << " ";
	*/

	if (dab < radius)
		return (mju * dh / cmax + (1 - mju) * dab);
	else
		return (mju * dh / cmax + (1 - mju) * dab) + mju * pow(dab - radius, l);
}


double nARC(const Point3D& a, const Point3D& b, const TVector <double>& na, const TVector <double>& nb, const double radius)
{
	//Norm: combined arc criterion, tangent model
	const double dab = nL2(a, b, na, nb, radius);

	//Identical points
	if (dab < MIN_FLOAT)
		return 0;

	//Compute abn criterion
	const double dabn = abn(na, nb);
	const double dh = 0.5 * dabn * dab;

	//Constrained (hybrid) metric
	const double cmax = M_PI / 4.0;;

	if (dab < radius)
		return (mju * dh / cmax + (1 - mju) * dab);
	else
		return (mju * dh / cmax + (1 - mju) * dab) + mju * pow(dab - radius, l);
}


double nARC2(const Point3D& a, const Point3D& b, const TVector <double>& na, const TVector <double>& nb, const double radius)
{
	//Norm: combined ARC2 criterion, secant model
	const double dab = nL2(a, b, na, nb, radius);

	//Identical points
	if (dab < MIN_FLOAT)
		return 0;

	//Compute ablp criterion
	const double dablp = ablp(a, b, na, nb);
	const double dh = sin(0.5 * dablp) * dab;

	//Constrained (hybrid) metric
	const double cmax = M_PI / 4.0;

	if (dab < radius)
		return (mju * dh / cmax + (1 - mju) * dab);
	else
		return (mju * dh / cmax + (1 - mju) * dab) + mju * pow(dab - radius, l);
}


double nDFP(const Point3D & a, const Point3D & b, const TVector <double>& na, const TVector <double>& nb, const double radius)
{
	//Norm: distance from planes (DFP)
	const double dab = nL2(a, b, na, nb, radius);

	//Identical points
	if (dab < MIN_FLOAT)
		return 0;

	//Compute dfp criterion
	double dh = dfp(a, b, na, nb);

	//Constrained (hybrid) metric
	const double cmax = sqrt(2.0);
	if (dab < radius)
		return (mju * dh / cmax + (1 - mju) * dab);
	else
		return (mju * dh / cmax + (1 - mju) * dab) + mju * pow(dab - radius, l);
}


void getNearestF(const Point3D & u, const TVector <Facility>& F, const TVector2D <double>& AN, const pfnorm& fnorm, const double radius, int& i_fmin, double& d_fmin)
{
	//Find closest facility in F closest to u
	const int n = F.size();

	//Empty list
	if (n == 0)
	{
		i_fmin = -1; d_fmin = -1;
		return;
	}

	//Initialize
	i_fmin = 0;
	d_fmin = fnorm(u, F[0].u, AN[u.id], AN[F[0].u.id], radius) + 2.0 * F[0].m;

	//Browse all facilities
	for (int i = 1; i < n; i++)
	{
		//Compute the given norm
		const double d_f = fnorm(u, F[i].u, AN[u.id], AN[F[i].u.id], radius) + 2.0 * F[i].m;

		//Actualize the nearest factory
		if (d_f < d_fmin)
		{
			d_fmin = d_f;
			i_fmin = i;
		}
	}

	//d_fmin = fnorm(u, F[i_fmin].u, AN[u.id], AN[F[i_fmin].u.id], radius) + 2.0 * F[i_fmin].m;
}


void getNearestPoint(const Point3D& u, const TVector <Point3D>& U, const TVector2D <double>& AN, const TVector <double>& FC, const TVector<int> &IF, const double rb, const int j, const pfnorm& fnorm, const double radius, int& i_fmin, double& fd_min)
{
	//Find closest facility in F closest to u
	const int n = U.size();

	//Empty list
	if (n == 0)
	{
		i_fmin = -1; fd_min = -1;
		return;
	}

	//Initialize
	i_fmin = -1;
	fd_min = MAX_FLOAT;

	//Browse all facilities
	for (int i = 0; i <= j; i++)
	{
		//Point is not a facility
		if (IF[i] == 0)
		{
			//Compute distance
			const double d_f = fnorm(u, U[i], AN[u.id], AN[U[i].id], radius);

			//Point inside the ball
			if (d_f <= rb)
			{
				//Compute the given norm
				const double fd = d_f + FC[i];

				//Actualize the nearest factory
				if (fd < fd_min)
				{
					fd_min = fd;
					i_fmin = i;
				}
			}
		}
	}
}


void findAllKNN(const TVector <Point3D>& points, const TVector <Point3D>& qpoints, const TVector<int>& K, TVector2D <size_t>& knn_id, TVector2D <float>& knn_dist)
{
	//Find all k nearest neighbors to any qpoint in points
	const int n = points.size(), nq = qpoints.size(),  nk = K.size();

	std::cout << ">> Convert points to the cloud: ";

	//Convert points to the cloud
	PointCloud<double> pointsc;
	pointsc.pts.resize(n);
	for (int i = 0; i < n; i++)
	{
		pointsc.pts[i].x = points[i].x;
		pointsc.pts[i].y = points[i].y;
		pointsc.pts[i].z = points[i].z;
	}

	std::cout << "OK \n";

	//Create kd tree
	std::cout << ">> Building KD-tree: ";
	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud<double> >, PointCloud<double>, 3> my_kd_tree_t;
	my_kd_tree_t  index(3, pointsc, KDTreeSingleIndexAdaptorParams(10));
	index.buildIndex();
	std::cout << "OK \n";

	//Perform knn search
	std::cout << ">> Perform knn-search: ";
	for (int i = 0; i < nq; i++)
	{
		//Fixed or variable amount of neighbors
		const int k = max((nk == 1 ? K[0] : K[i]), 1);

		//Auxilliary data structures
		TVector<size_t> ret_indices(k);
		TVector<float> out_dists_sqr(k);

		//Initialize result set
		nanoflann::KNNResultSet<float> resultSet(k);
		resultSet.init(&ret_indices[0], &out_dists_sqr[0]);

		//Perform knn search
		double query_pt[3] = { qpoints[i].x, qpoints[i].y, qpoints[i].z };
		index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

		//Add to the list
		knn_id.push_back(ret_indices);
		knn_dist.push_back(out_dists_sqr);
	}

	std::cout << "OK \n";
}


void regressionPlane(const TVector <Point3D>& XA, double& a, double& b, double& c, double& d, double &sigma, double& var)
{
	//Compute regression plane using SVD
	const int n = XA.size();

	//Not enough points for the regression plane
	if (n < 3)
		return;

	Eigen::MatrixXd X(n, 3);

	//Convert points to matrix
	for (int i = 0; i < n; i++)
	{
		X(i, 0) = XA[i].x;
		X(i, 1) = XA[i].y;
		X(i, 2) = XA[i].z;
	}

	//Compute centroid
	Eigen::RowVector3d C = Eigen::RowVector3d(X.col(0).mean(), X.col(1).mean(), X.col(2).mean()).transpose();

	//Subtract centroid
	Eigen::MatrixXd CX = X.rowwise() - C;

	//Compute SVD
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CX, Eigen::ComputeFullV | Eigen::ComputeFullU/*Eigen::ComputeThinU | Eigen::ComputeThinV*/);
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::MatrixXd S = svd.singularValues();
	Eigen::MatrixXd V = svd.matrixV();

	//Compute sigma: cummulated distance of points from the plane
	sigma = S(2, 0);
	
	//Compute variability of normals S[min] / (S[min] + S[mid] + S[max])
	var = S(2, 0) / (S(0, 0) + S(1, 0) + S(2, 0));
	
	//std::cout << "\n\n" << S;

	//Compute parameters a, b, c of the regression plane (last column of V)
	a = V(0, 2);
	b = V(1, 2);
	c = V(2, 2);

	//Compute parameter d of the regression plane
	d = -(a * C(0, 0) + b * C(0, 1) + c * C(0, 2));
}


TVector <Point3D> findPointsCloserThanMedian(const TVector <Point3D>& XA)
{
	//Find points closer than median distance from the plane
	TVector <Point3D> XAR;

	//Compute the regression plane
	double a = 0, b = 0, c = 0, d = 0, sigma = 0, var = 0;
	regressionPlane(XA, a, b, c, d, sigma, var);

	//Compute distances from the plane
	TVector <double> dists;
	for (auto p : XA)
		dists.push_back(getPointPlaneDist(p.x, p.y, p.z, a, b, c, d));

	//Find the median distance
	std::nth_element(dists.begin(), dists.begin() + dists.size() / 2, dists.end());
	const double dist_median = dists[dists.size() / 2];

	//Find points closer than median
	for (int i = 0; i < dists.size(); i++)
		if (dists[i] < dist_median)
			XAR.push_back(XA[i]);

	return XAR;
}


void regressionPlaneRobust(const TVector <Point3D>& XA, double& a, double& b, double& c, double& d, double& sigma, double& var)
{
	//Robust computation of regression plane using SVD and median
	const int n = XA.size(), iter_max = n / 2;
	int iter = 0;
	const double abn_max = 2.0 * M_PI / 180;
	double abnr = 0;

	//Initial computation of the regression plane
	regressionPlane(XA, a, b, c, d, sigma, var);

	//Not enough neighbors: stop the iteration process
	if (n < 5)
		return;

	//Iteration process
	TVector <Point3D> XAR;
	do
	{
		//Find points closer than median
		XAR = findPointsCloserThanMedian(XAR);

		//Recompute plane
		double an = 0, bn = 0, cn = 0, dn = 0, sigman = 0;
		regressionPlane(XAR, an, bn, cn, dn, sigma, var);

		//Not enough relevant points
		if (XAR.size() < 5)
			break;

		//Compute abn between new and old plane
		TVector <double> na = { a, b, c, d }, nb = { an, bn, cn, dn };
		abnr = abn(na, nb);

		//Assign normal vector
		a = an; b = bn; c = cn; d = dn; sigma = sigman;

		iter++;
	}while ((abnr > abn_max) && (iter < iter_max));

}


void assignPointsToFacilities(const TVector <Point3D>& U, const TVector <int>& IF, const TVector2D <double>& AN, const double radius, const pfnorm& fnorm, TVector <Facility>& F)
{
	//Assign point of the cloud to the nearest facility according to a pseudo norm
	//Assign vertices
	for (int i = 0; i < U.size(); i++)
	{
		//Point is not facility
		if (IF[i] == 0)
		{
			//Find nearest facility
			int imin;
			double gmin;
			getNearestF(U[i], F, AN, fnorm, radius, imin, gmin);

			//Get the nearest facility
			Facility& fn = F[imin];

			//Assign the point to the nerest facility
			fn.C.push_back(U[i]);
		}
	}
}


void assignPointsToFacilitiesL2(const TVector <Point3D>& U, const TVector <int> &IF, TVector <Facility>& F)
{
	//Assign point of the cloud to the nearest facility according to L2 norm
	//Fast and approximate solution using nano flann
	const int nf = F.size();

	//Convert facilities to the cloud
	PointCloud<double> pointsc;
	pointsc.pts.resize(nf);

	for (int i = 0; i < nf; i++)
	{
		pointsc.pts[i].x = F[i].u.x;
		pointsc.pts[i].y = F[i].u.y;
		pointsc.pts[i].z = F[i].u.z;
	}

	//Create KD tree
	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud<double> >, PointCloud<double>, 3> my_kd_tree_t;
	my_kd_tree_t  index(3, pointsc, KDTreeSingleIndexAdaptorParams(10));
	index.buildIndex();

	//Perform knn search
	int k = 1, n = U.size();
	for (int i = 0; i < n; i++)
	{
		//Point is not facility
		if (IF[i] == 0)
		{
			//Auxilliary data structures
			size_t ret_index;
			float out_dist_sqr;

			//Initialize result set
			nanoflann::KNNResultSet<float> resultSet(k);
			resultSet.init(&ret_index, &out_dist_sqr);

			//Perform knn search
			double query_pt[3] = { U[i].x, U[i].y, U[i].z };
			index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

			//Get the nearest facility
			Facility& fn = F[ret_index];

			//Assign the point to the nerest facility
			fn.C.push_back(U[i]);
		}
	}
}


void getAveragePointNormal(TVector <Point3D>& U, const TVector2D <size_t> &knn_id, TVector2D <double>& AN)
{
	//Compute average point normal using SVD decomposition, regression error and ABN
	const int n = U.size();

	std::cout << ">> SVD decomposition: ";
	for (int i = 0; i < n; i++)
	{
		//Add all neighbors to the list
		TVector <Point3D> KNN;
		for (size_t index : knn_id[i])
			KNN.push_back(U[index]);

		//Compute average normal
		double a = 0, b = 0, c = 0, d = 0, sigma = 0, var = 0;
		regressionPlane(KNN, a, b, c, d, sigma, var);
		//U[i].id = i;
		//Add average normal vector and regression error to the list
		//ABN will be computed later
		AN.push_back({ a, b, c, d, sigma, 0, var });
	}

	//Compute max ABN for points and their neighbors
	for (int i = 0; i < n; i++)
	{
		//Find max ABN
		const size_t index0 = knn_id[i][0];
		double abn_max = 0;
		for (size_t index : knn_id[i])
		{
			//const double abn_j = abn(AN[index0], AN[index]);
			const double abn_j = nL2(U[index0], U[index], AN[index0], AN[index], 0.0) * sin(abn(AN[index0], AN[index]) / 2);

			//Actualize maximum
			if (abn_j > abn_max)
				abn_max = abn_j;
		}

		//Store ABN max value to the list
		AN[i][5] = abn_max;
	}

	std::cout << "OK \n";
}


float genRandom(const float a, const float b)
{
	//Create random number inside interval
	const float diff = b - a;
	return (((float)rand() / RAND_MAX) * diff) + a;
}


void resamplePointCloud(const TVector <Point3D>& points, const int k, TVector <Point3D>& points_res)
{
	//Resample point cloud: each point is the average of k nearest neighbors
	TVector2D <size_t> knn_id;
	TVector2D <float> knn_dist;
	TVector <int> K(1, k);

	std::cout << ">> Resample point cloud: \n";

	//Find all k-nearest neighbors
	findAllKNN(points, points, K, knn_id, knn_dist);

	//Compute average of k points
	//srand(time(NULL));
	for (int i = 0; i < points.size(); i++)
	{
		//Add all neighbors to the list
		double x_mean = points[knn_id[i][0]].x, y_mean = points[knn_id[i][0]].y, z_mean = points[knn_id[i][0]].z;
		for (int j = 1; j < knn_id[i].size(); j++)
		{
			x_mean += points[knn_id[i][j]].x; //+ genRandom(-0.01, 0.01);
			y_mean += points[knn_id[i][j]].y; //+ genRandom(-0.01, 0.01);
			z_mean += points[knn_id[i][j]].z; //+ genRandom(-0.01, 0.01);
		}

		//Compute mean of coordinates
		x_mean /= k;
		y_mean /= k;
		z_mean /= k;

		//Create new point
		Point3D p(points[i].id, x_mean, y_mean, z_mean, points[i].r, points[i].g, points[i].b);
		points_res.push_back(p);

		if (i % 25000 == 0)
			std::cout << ".";
	}

}


void recomputeFacilityCostsFNFL(const double f, const double crit, double rat, const TVector2D <double>& AN, const pfnorm& fnorm, TVector <double>& FC)
{
	//Recompute facility cost for non-uniform clustering
	const double eps = 0.0001;

	for (auto an : AN)
	{
		//ABN
		if (fnorm == nDIS || fnorm == nARC || fnorm == nARC2)
		{
			//Compute new facility cost
			const double f_new = max(min(f * f / (an[5] + eps), rat * f), f / rat);
			FC.push_back(f_new);
			//FC.push_back(f);
		}

		//DFP
		else if (fnorm == nDFP)
		{
			const double f_new = max(min(f * f / (an[4] + eps), rat * f), f / rat);
			FC.push_back(f_new);
			//FC.push_back(f);
		}
		
		//L2-ABN
		else
		{
			const double f_new = max(min(f * crit / (an[5] + eps), rat * f), f / rat);
			FC.push_back(f_new);
			//FC.push_back(f);
		}

	}

	//for (auto fc : FC)
	//	std::cout << fc << " ";
}


double computeTotalCost(const TVector <Facility>& F, const pfnorm& fnorm, const TVector2D <double>& AN, const double radius)
{
	//Compute total cost of the clusterization
	double total_cost = 0;

	for (const auto f : F)
	{
		//Add the facility cost
		total_cost += f.fc;

		//Process all assigned points
		for (const auto u : f.C)
		{
			total_cost += fnorm(f.u, u, AN[f.u.id], AN[u.id], radius);
		}
	}

	return total_cost;
}

void findEdgePoints(const TVector2D <double>& AN, const double var_max, const TVector <Point3D>& P, TVector <Point3D>& PE)
{
	//Find points located on edges according to the normal variability
	for (int i = 0; i < P.size(); i++)
	{
		//std::cout << AN[i][6] << ' ';

		if (AN[i][6] > var_max)
		{
			PE.push_back(P[i]);
		}
	}
}


void printHeader(const std::string& file_name)
{
	//Print header with the local statistics
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name, std::ofstream::out | std::ofstream::app);

		//Write file header
		file << std::setw(8) << "N"
			<< std::setw(8) << "N_F"
			<< std::setw(10) << "N_F/N"
			<< std::setw(12) << "N_pts_min"
			<< std::setw(12) << "N_pts_max"
			<< std::setw(12) << "N_pts_aver"
			<< std::setw(12) << "Time"
			<< std::setw(8) << "N_del"
			<< std::setw(12) << "Cl_rad_min"
			<< std::setw(12) << "Cl_rad_max"
			<< std::setw(12) << "Cl_rad_aver"
			<< std::setw(12) << "F_cost_min"
			<< std::setw(12) << "F_cost_max"
			<< std::setw(12) << "F_cost_aver"
			<< std::setw(12) << "F_cost_total"
			<< std::setw(12) << "dABN_min"
			<< std::setw(12) << "dABN_max"
			<< std::setw(12) << "dABN_aver"
			<< std::setw(12) << "dDFP_min"
			<< std::setw(12) << "dDFP_max"
			<< std::setw(12) << "dDFP_aver"
			<< std::setw(12) << "d_edge_min"
			<< std::setw(12) << "d_edge_max"
			<< std::setw(12) << "d_edge_aver"
			<< std::setw(12) << "ovrlp_min"
			<< std::setw(12) << "ovrlp_max"
			<< std::setw(12) << "ovrlp_aver"
			<< std::endl;

		file.close();
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}


void printStatisticsInLine(const double vals [], const std::string& file_name)
{
	//Print statistics in line
	std::ofstream file;

	try
	{
		//Open file
		file.open(file_name, std::ofstream::out | std::ofstream::app);

		//Set parameters
		file << std::setw(4) << std::fixed << std::right;

		//Print results
		file << std::setw(8) << std::setprecision(0) << vals[0]
			<< std::setw(8) << std::setprecision(0) << vals[1]
			<< std::setw(10) << std::setprecision(2) << vals[2]
			<< std::setw(12) << std::setprecision(2) << vals[3]
			<< std::setw(12) << std::setprecision(2) << vals[4]
			<< std::setw(12) << std::setprecision(2) << vals[5]
			<< std::setw(12) << std::setprecision(2) << vals[6]
			<< std::setw(8) << std::setprecision(0) << vals[7]
			<< std::setw(12) << std::setprecision(3) << vals[8]
			<< std::setw(12) << std::setprecision(3) << vals[9]
			<< std::setw(12) << std::setprecision(3) << vals[10]
			<< std::setw(12) << std::setprecision(3) << vals[11]
			<< std::setw(12) << std::setprecision(3) << vals[12]
			<< std::setw(12) << std::setprecision(3) << vals[13]
			<< std::setw(12) << std::setprecision(3) << vals[14]
			<< std::setw(12) << std::setprecision(3) << vals[15]
			<< std::setw(12) << std::setprecision(3) << vals[16]
			<< std::setw(12) << std::setprecision(4) << vals[17]
			<< std::setw(12) << std::setprecision(4) << vals[18]
			<< std::setw(12) << std::setprecision(4) << vals[19]
			<< std::setw(12) << std::setprecision(4) << vals[20]
			<< std::setw(12) << std::setprecision(4) << vals[21]
			<< std::setw(12) << std::setprecision(4) << vals[22]
			<< std::setw(12) << std::setprecision(4) << vals[23]
			<< std::setw(12) << std::setprecision(4) << vals[24]
			<< std::setw(12) << std::setprecision(4) << vals[25]
			<< std::setw(12) << std::setprecision(4) << vals[26]
			<< std::endl;

		file.close();
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << file_name << '\n';
	}
}


void computeGlobalStatistics(const TVector<std::string> &file_name_subsets, const std::string& fnorm_text, const std::string &method_text, const double f, const double dc, const std::string &output_file_text)
{
	//Compute global statistics for all kd subsets
	TVector <int> n_points, n_deleted, n_facilities;
	TVector <double> time_diff, n_clust_points_min, n_clust_points_max, n_clust_points_aver, clust_radius_min, clust_radius_max, clust_radius_aver,
		fac_cost_min, fac_cost_max, fac_cost_aver, fac_cost_total, fac_abn_diff_min, fac_abn_diff_max, fac_abn_diff_aver, fac_dfp_diff_min, fac_dfp_diff_max, 
		fac_dfp_diff_aver, fac_dist_edge_min, fac_dist_edge_max, fac_dist_edge_aver, fac_overlap_min, fac_overlap_max, fac_overlap_aver;
	
	//Print header
	printHeader(output_file_text);

	unsigned n_subsets = file_name_subsets.size();
	for (auto f_name : file_name_subsets)
	{
		//List of results
		std::string result_file_subset = f_name + "_" + fnorm_text + "_l" + std::to_string(l) + "_" + method_text + "_f" + std::to_string(f) + "_dc" + std::to_string(dc) + "_results.log";

		//Output file exists
		if (filesystem::exists(result_file_subset))
		{
			//Process all facilities
			int i = 0, n_pts, n_facil, deleted;
			double  n_c_ratio, nc_point_min, nc_point_max, nc_point_aver, dtime, c_radius_aver, c_radius_min, c_radius_max, f_cost_aver,
				f_cost_min, f_cost_max, f_cost_total, f_abn_diff_aver, f_abn_diff_min, f_abn_diff_max, f_dfp_diff_aver, f_dfp_diff_min, 
				f_dfp_diff_max, f_dist_edge_min, f_dist_edge_max, f_dist_edge_aver, f_overlap_min, f_overlap_max, f_overlap_aver;

			//Read a single line from file
			ifstream infile(result_file_subset);
			
			try
			{
				infile >> n_pts >> n_facil >> n_c_ratio >> nc_point_min >> nc_point_max >> nc_point_aver >> dtime >> deleted
					>> c_radius_min >> c_radius_max >> c_radius_aver >> f_cost_min >> f_cost_max >> f_cost_aver >> f_cost_total
					>> f_abn_diff_min >> f_abn_diff_max >> f_abn_diff_aver >> f_dfp_diff_min >> f_dfp_diff_max >> f_dfp_diff_aver >> f_dist_edge_min
					>> f_dist_edge_max >> f_dist_edge_aver >> f_overlap_min >> f_overlap_max >> f_overlap_aver;
			}

			//Throw exception
			catch (std::ostream::failure e)
			{
				std::cerr << "Can not read the file> :" << result_file_subset << '\n';
			}

			//Results
			const double results[] =
			{
				n_pts,
				n_facil,
				n_c_ratio,
				nc_point_min,
				nc_point_max,
				nc_point_aver,
				dtime,
				deleted,
				c_radius_min,
				c_radius_max,
				c_radius_aver,
				f_cost_min,
				f_cost_max,
				f_cost_aver,
				f_cost_total,
				f_abn_diff_min,
				f_abn_diff_max,
				f_abn_diff_aver,
				f_dfp_diff_min,
				f_dfp_diff_max,
				f_dfp_diff_aver,
				f_dist_edge_min,
				f_dist_edge_max,
				f_dist_edge_aver,
				f_overlap_min, 
				f_overlap_max, 
				f_overlap_aver
			};

			//Print results
			printStatisticsInLine(results, output_file_text);

			//Add analyzed parameters to the list
			n_points.push_back(n_pts);
			n_facilities.push_back(n_facil);
			n_deleted.push_back(deleted);
			time_diff.push_back(dtime);

			n_clust_points_min.push_back(nc_point_min);
			n_clust_points_max.push_back(nc_point_max);
			n_clust_points_aver.push_back(nc_point_aver);

			fac_cost_min.push_back(f_cost_min);
			fac_cost_max.push_back(f_cost_max);
			fac_cost_aver.push_back(f_cost_aver);
			fac_cost_total.push_back(f_cost_total);

			clust_radius_min.push_back(c_radius_min);
			clust_radius_max.push_back(c_radius_max);
			clust_radius_aver.push_back(c_radius_aver);

			fac_abn_diff_min.push_back(f_abn_diff_min);
			fac_abn_diff_max.push_back(f_abn_diff_max);
			fac_abn_diff_aver.push_back(f_abn_diff_aver);

			fac_dfp_diff_min.push_back(f_dfp_diff_min);
			fac_dfp_diff_max.push_back(f_dfp_diff_max);
			fac_dfp_diff_aver.push_back(f_dfp_diff_aver);

			fac_dist_edge_min.push_back(f_dist_edge_min);
			fac_dist_edge_max.push_back(f_dist_edge_max);
			fac_dist_edge_aver.push_back(f_dist_edge_aver);

			fac_overlap_min.push_back(f_overlap_min);
			fac_overlap_max.push_back(f_overlap_max);
			fac_overlap_aver.push_back(f_overlap_aver);
		}
	}

	//Compute average results
	TVector <double> n_n_clust;
	std::transform(n_facilities.begin(), n_facilities.end(), n_points.begin(), back_inserter(n_n_clust), std::divides<double>());

	const double results_a[] =
	{
		accumulate(n_points.begin(), n_points.end(), 0.0) / n_subsets,
		accumulate(n_facilities.begin(), n_facilities.end(), 0.0) / n_subsets,
		accumulate(n_n_clust.begin(), n_n_clust.end(), 0.0) / n_subsets * 100,
		accumulate(n_clust_points_min.begin(), n_clust_points_min.end(), 0.0) / n_subsets,
		accumulate(n_clust_points_max.begin(), n_clust_points_max.end(), 0.0) / n_subsets,
		accumulate(n_clust_points_aver.begin(), n_clust_points_aver.end(), 0.0) / n_subsets,
		accumulate(time_diff.begin(), time_diff.end(), 0.0) / n_subsets,
		accumulate(n_deleted.begin(), n_deleted.end(), 0.0) / n_subsets,
		accumulate(clust_radius_min.begin(), clust_radius_min.end(), 0.0) / n_subsets,
		accumulate(clust_radius_max.begin(), clust_radius_max.end(), 0.0) / n_subsets,
		accumulate(clust_radius_aver.begin(), clust_radius_aver.end(), 0.0) / n_subsets,
		accumulate(fac_cost_min.begin(), fac_cost_min.end(), 0.0) / n_subsets,
		accumulate(fac_cost_max.begin(), fac_cost_max.end(), 0.0) / n_subsets,
		accumulate(fac_cost_aver.begin(), fac_cost_aver.end(), 0.0) / n_subsets,
		accumulate(fac_cost_total.begin(), fac_cost_total.end(), 0.0) / n_subsets,
		accumulate(fac_abn_diff_min.begin(), fac_abn_diff_min.end(), 0.0) / n_subsets,
		accumulate(fac_abn_diff_max.begin(), fac_abn_diff_max.end(), 0.0) / n_subsets,
		accumulate(fac_abn_diff_aver.begin(), fac_abn_diff_aver.end(), 0.0) / n_subsets,
		accumulate(fac_dfp_diff_min.begin(), fac_dfp_diff_min.end(), 0.0) / n_subsets,
		accumulate(fac_dfp_diff_max.begin(), fac_dfp_diff_max.end(), 0.0) / n_subsets,
		accumulate(fac_dfp_diff_aver.begin(), fac_dfp_diff_aver.end(), 0.0) / n_subsets,
		accumulate(fac_dist_edge_min.begin(), fac_dist_edge_min.end(), 0.0) / n_subsets,
		accumulate(fac_dist_edge_max.begin(), fac_dist_edge_max.end(), 0.0) / n_subsets,
		accumulate(fac_dist_edge_aver.begin(), fac_dist_edge_aver.end(), 0.0) / n_subsets,
		accumulate(fac_overlap_min.begin(), fac_overlap_min.end(), 0.0) / n_subsets,
		accumulate(fac_overlap_max.begin(), fac_overlap_max.end(), 0.0) / n_subsets,
		accumulate(fac_overlap_aver.begin(), fac_overlap_aver.end(), 0.0) / n_subsets
	};

	//Write empty line
	std::ofstream file;
	try
	{
		file.open(output_file_text, std::ofstream::out | std::ofstream::app);
		file << std::endl;
		file.close();
	}

	//Throw exception
	catch (std::ostream::failure e)
	{
		std::cerr << "Can not write the file> :" << output_file_text << '\n';
	}

	//Print average results
	printStatisticsInLine(results_a, output_file_text);
}


void computeLocalStatistics(TVector <Point3D>& U, TVector <Point3D>& US, const double f, const int k, const double radius, const pfnorm& fnorm, const int deleted, const double time, const TVector2D <double>& AN, const TVector <Facility>& F, const std::string& file_name)
{
	//Compute local statistics on the kd-subsets
	TVector2D <size_t> knn_id_f;
	TVector2D <float> knn_dist_f;
	TVector <int> K;

	//Add amount of points in cluster
	for (auto f : F)
		K.push_back(f.C.size());
	
	//Find all k-nearest neighbors
	findAllKNN(US, US, K, knn_id_f, knn_dist_f);

	//Compute average normals for clusters
	TVector <double> FCS;
	TVector2D <double> ANS;
	getAveragePointNormal(US, knn_id_f, ANS);

	//Process all facilities
	int i = 0, n_facil = F.size();
	double nc_point_aver = 0, nc_point_min = max_int, nc_point_max = 0, f_cost_aver = 0, f_cost_min = max_float, f_cost_max = 0,
		f_cost_total = 0, c_radius_aver = 0, c_radius_min = max_float, c_radius_max = 0, f_abn_diff_aver = 0, f_abn_diff_min = max_float,
		f_abn_diff_max = 0, f_dfp_diff_aver = 0, f_dfp_diff_min = max_float, f_dfp_diff_max = 0, f_dist_edge_min = max_float, 
		f_dist_edge_max = 0, f_dist_edge_aver = 0, f_overlap_min = max_float, f_overlap_max = 0, f_overlap_aver = 0;

	for (auto f : F)
	{
		//Add amount of points of the cluster
		nc_point_min = min(nc_point_min, static_cast<double>(f.C.size()) + 1);
		nc_point_max = max(nc_point_max, static_cast<double>(f.C.size()) + 1);
		nc_point_aver += f.C.size() + 1;

		//Add facility cost
		f_cost_min = min(f_cost_min, f.fc);
		f_cost_max = max(f_cost_max, f.fc);
		f_cost_aver += f.fc;
		f_cost_total += f.fc;

		//Total cost
		//Process all assigned points
		for (const auto u : f.C)
		{
			f_cost_total += fnorm(f.u, u, AN[f.u.id], AN[u.id], radius);
		}


		//Analyze radius, xmin, xmax, ymin, ymax of the cluster
		double c_radius = 0;
		double xmin = MAX_FLOAT, ymin = MAX_FLOAT, zmin = MAX_FLOAT, xmax = -MAX_FLOAT, ymax = -MAX_FLOAT, zmax = -MAX_FLOAT;
		for (auto& c : f.C)
		{
			//Find vertices of the min-max box
			if (c.x < xmin)
				xmin = c.x;
			if (c.y < ymin)
				ymin = c.y;
			if (c.z < zmin)
				zmin = c.z;
			if (c.x > xmax)
				xmax = c.x;
			if (c.y > ymax)
				ymax = c.y;
			if (c.z > zmax)
				zmax = c.z;

			//Find max radius
			c_radius = max(dist(f.u.x, f.u.y, f.u.z, c.x, c.y, c.z), c_radius);
		}

		c_radius_min = min(c_radius_min, c_radius);
		c_radius_max = max(c_radius_max, c_radius);
		c_radius_aver += c_radius;

		//Analyze DFP differencies in facilities
		//Use facilities with > 3 points
		double a, b, c, d, ddfp = 0, var = 0;
		if (f.C.size() > 3)
		{
			regressionPlane(f.C, a, b, c, d, ddfp, var);

			f_dfp_diff_min = min(f_dfp_diff_min, ddfp);
			f_dfp_diff_max = max(f_dfp_diff_max, ddfp);
			f_dfp_diff_aver += ddfp;
		}
		
		//Analyze normal variance in facilities
		double dabn = 0;
		for (auto & c : f.C)
			dabn = max(abn(AN[f.u.id], AN[c.id]), dabn);  //K3
		
		//std::cout << f.u.x << ' '  << f.u.y << ' ' << f.u.z << ' ' << dabn << '\n';
		f_abn_diff_min = min(f_abn_diff_min, dabn);
		f_abn_diff_max = max(f_abn_diff_max, dabn);
		f_abn_diff_aver += dabn;
		
		//Cluster overlaps
		double f_overlap = 0;
		for (auto ff : F)
		{
			//Point inside min-max box generated by cluster
			if ((ff.u.x >= xmin) && (ff.u.x <= xmax) && (ff.u.y >= ymin) && (ff.u.y <= ymax) && (ff.u.z >= zmin) && (ff.u.z <= zmax))
				f_overlap++;
		}

		f_overlap_min = min(f_overlap_min, f_overlap);
		f_overlap_max = max(f_overlap_max, f_overlap);
		f_overlap_aver += f_overlap;

		//Increment facility index
		i++;
	}

	//Compute new normals for edges
	TVector2D <size_t> knnf_id;
	TVector2D <float> knnf_dist;
	TVector <int> Kf(1, k);
	TVector2D <double> ANF;
	findAllKNN(US, US, Kf, knnf_id, knnf_dist);
	getAveragePointNormal(US, knnf_id, ANF);

	//Find edges, original/simplified point cloud: edge 0.10 (cube), 0.15 (cone)
	TVector <Point3D> PE, PFE;
	findEdgePoints(AN, 0.10, U, PE);
	findEdgePoints(ANF, 0.10, US, PFE);

	//Find point on the simplified edge closest to the original edge point
	TVector2D <size_t> knne_id;
	TVector2D <float> knne_dist;
	TVector <int> Ke(1, 1);
	findAllKNN(PFE, PE, Ke, knne_id, knne_dist);

	//At least one edge point has been found
	if (PFE.size() > 0)
	{
		for (auto d : knne_dist)
		{
			f_dist_edge_min = min(f_dist_edge_min, (double)sqrt(d[0]));
			f_dist_edge_max = max(f_dist_edge_max, (double)sqrt(d[0]));
			f_dist_edge_aver += sqrt(d[0]);

			//std::cout << sqrt(d[0]) << "  " << f_dist_edge_min << "  " << f_dist_edge_max << "  " << f_dist_edge_aver << '\n';
		}

		f_dist_edge_aver = f_dist_edge_aver / knne_dist.size();
	}

	//std::cout << "dist: " << f_dist_edge_aver << " " << knne_dist.size() << '\n';
	/*
	TVector <Point3D> start, end;
	double dd = 0;
	for (int i = 0; i < knne_dist.size(); i++)
	{
		start.push_back(PE[i]);
		end.push_back(PFE[knne_id[i][0]]);

		dd += dist(PE[i].x, PE[i].y, PE[i].z, PFE[knne_id[i][0]].x, PFE[knne_id[i][0]].y, PFE[knne_id[i][0]].z);
	}

	//std::cout << "dist2: " << dd/PFE.size() << '\n';
	*/
	//Export edges
	//DXFExport::exportPointsToDXF("edges_input.dxf", PE, 5, 1);
	//DXFExport::exportPointsToDXF("edges_output.dxf", PFE, 5, 1);
	//DXFExport::exportLinesToDXF("lines_output.dxf", start, end, 5, 1);
	
	//Create results: conversion to mm (multiplication by 1000)
	const double results[] =
	{
		U.size(),
		n_facil,
		(n_facil * 100.0) / U.size(),
		nc_point_min,
		nc_point_max,
		nc_point_aver / n_facil,
		time,
		deleted,
		1000 * c_radius_min,
		1000 * c_radius_max,
		1000 * c_radius_aver / n_facil,
		f_cost_min,
		f_cost_max,
		f_cost_aver / n_facil,
		f_cost_total,
		f_abn_diff_min * 180 / M_PI,
		f_abn_diff_max * 180 / M_PI,
		f_abn_diff_aver / n_facil * 180 / M_PI,
		1000 * f_dfp_diff_min,
		1000 * f_dfp_diff_max,
		1000 * f_dfp_diff_aver / n_facil,
		1000 * f_dist_edge_min,
		1000 * f_dist_edge_max,
		1000 * f_dist_edge_aver,
		f_overlap_min,
		f_overlap_max,
		f_overlap_aver / n_facil
	};

	//Print results
	printStatisticsInLine(results, file_name);
}



void updateClusters(const Point3D& p, const double radius, const pfnorm& fnorm, TVector <Facility>& F, TVector2D <double>& AN, const TVector <double> &FC)
{
//Incremental method, update of the clusterization takes one of four strategies:
//   1) Create new facility at p: S1
//   2) Add p to the nearest facility: S2
//   3) Assign parts of clusters closest to p: create cluster at p + partial reassignment: S3 a)
//   4) Assign all clusters closest to p: create cluster at p + fulll reassignment + deletition: S3 b)
	int i = 0, i_nearest = -1;

	//const double cb = computeTotalCost(F, fnorm, AN, radius);

	//Initial costs for different strategies
	double c_nearest = MAX_FLOAT;			   //Cost of the nearest facility, strategy S1
	double c_new = FC[p.id];			   //Cost for creation of the new facility, strategy S2
	double c_reassign_clusters = FC[p.id];		   //Cost for the cluster modification (full or partial reassignment to another facility), strategy S3 a) b)

	//Browse all facilities
	TVector2D <unsigned int> reassignment_type;        //Changed cluster id and type of reassignment operation (false = partial, true = full)
	for (const Facility& f : F)
	{
		const double dpf = fnorm(p, f.u, AN[p.id], AN[f.u.id], radius);                 //Distance between point p and facility
		double dc_all = FC[p.id] - f.fc;						//Cost diifference: reassignment of all cluster to  p - cost for the old facility deletition
		double dc_closer = FC[p.id] - dpf;						//Cost difference: reassignment of cluster points closer to p

		//Distance point and the current facility: Strategy S1
		if (dpf < c_nearest && dpf > 0)  //Actualize distance to the nearest facility
		{
			c_nearest = dpf;
			i_nearest = i;
		}

		//Browse all points of the facility
		unsigned int j = 0;
		TVector <unsigned int> reassigned_oper;					      //Reassign operation: 0 = partial to p, 1 = new cluster at p; r[0] = id cluster, r[1] = oper. type, r[>2] reassigned indices
		reassigned_oper.push_back(i); reassigned_oper.push_back(0);		      //Temporary set for operation 0 for ith cluster
		
		for (const Point3D& u : f.C)
		{
			const double dup = fnorm(u, p, AN[u.id], AN[p.id], radius);           //Distance of the cluster point u to the proposed new center p
			const double duf = fnorm(u, f.u, AN[u.id], AN[f.u.id], radius);       //Distance of the cluster point u to its cluster center f.u

			//Current cluster point u is closer to p: cost for the reassignment u to the new center p
			//Strategy S3 a)
			if (dup < duf)
			{
				//Compute new cost increment
				dc_closer += dup - duf;

				//Mark point u for the deletition: it will be reassigned
				reassigned_oper.push_back(j);
			}

			//Cost for the reassignment of any cluster point to p
			//Strategy S3 b)
			dc_all += dup - duf;

			//Increment j
			j++;
		}

		//Entire cluster is reassigned to p: does the newly created cluster at p decrese the cost?
		//Strategy S3 b)
		if ((dc_all < dc_closer) && (dc_all < min(dpf, c_new)))				     //Cost is better than assignment p to the facility
		{
			c_reassign_clusters = c_reassign_clusters + dc_all - FC[p.id] + dpf;	     //Expected cost increment
			reassigned_oper[1] = 1;					                     //Change type of operation to 1 (new cluster at p)
			reassignment_type.push_back(reassigned_oper);
		}

		//Reassign part of the cluster to p: is there any improvement?
		//Strategy S3 a)
		else if ((dc_closer < dc_all) && (dc_closer < min(dpf, c_new)))			      //Cost is better than assignment p to the facility
		{
			c_reassign_clusters = c_reassign_clusters + dc_closer - FC[p.id] + dpf;       //Expected increment
			reassignment_type.push_back(reassigned_oper);
		}

		//Increment i
		i++;
	}

	//Heuristic decision: find minimum cost increment and choose the optimal strategy (S1-S3 a) b))
	const double c_min = min(c_new, min(c_nearest, c_reassign_clusters));

	//S1: Create new facility at p is the cheapest
	if (c_min == c_new)
	{
		Facility fac(p, 0, FC[p.id]);
		F.push_back(fac);
		//std::cout << "C";
	}

	//S2: Assign p to the nearest facility is the cheapest
	else if (c_min == c_nearest)
	{
		F[i_nearest].C.push_back(p);
		//std::cout << "N";
	}

	//S3 a): Reasign entire clusters to p is the cheapest
	//S3 b): Reassign parts of clusters to p is the cheapest
	else
	{
		//Create new facility at p
		Facility fac_new(p, 0, FC[p.id]);

		//Process all cluster affected by the reassignment one by one
		for (const auto &r : reassignment_type)
		{
			//Get index of the changed cluster in the list
			const unsigned int k = r[0];

			//S3 a)
			//Reassign part of the cluster to the new facility at p: all closest points marked as deleted are deleted
			if (!r[1])
			{
				//Set flag for reassigned and deleted points
				for (unsigned int j = 2; j < r.size(); j++) F[k].C[r[j]].deleted = true;

				//std::cout << "RP ";

				//Copy (and reassign) marked points to the new facility at p
				std::copy_if(F[k].C.begin(), F[k].C.end(), std::back_inserter(fac_new.C), isPointMarkedToDelete);

				//Delete marked points from the the old cluster
				F[k].C.erase(std::remove_if(F[k].C.begin(), F[k].C.end(), isPointMarkedToDelete), F[k].C.end());

				//Reset flag for reassigned  points
				for (auto & p : fac_new.C) p.deleted = false;
			}

			//S3 b)
			//Assign entire cluster to the new facility at p
			else
			{
				//std::cout << "RA ";

				//Insert all points to the new facility at p
				fac_new.C.insert(fac_new.C.end(), F[k].C.begin(), F[k].C.end());

				//Add old facility center to the facility at p
				fac_new.C.push_back(F[k].u);

				//Mark old cluster for the deletetion
				F[k].deleted = true;
			}
		}

		//Add new facility at p to the list
		F.push_back(fac_new);
	}
	
	//Delete marked facilities of the old cluster (now it is reassigned)
	F.erase(std::remove_if(F.begin(), F.end(), isFacilityMarkedToDelete), F.end());

	/*
	const double ca = computeTotalCost(F, fnorm, AN, radius);
	const double diff = ca - cb;
	std::cout << "CB: " << cb << ", CA: " << ca << ", D: " << diff << ", DE: " << c_min << '\n';
	*/
}



void IFL(TVector <Point3D>& U, const double f, const double alpha, const double radius, const int k, const pfnorm& fnorm, const bool ifl_non_uniform, TVector <Facility>& F, TVector2D <double>& AN, int& deleted)
{
	//Perform uniform facility cost clustering by FFL according to a given norm
	//Proposed incremental algorithms selecting one of four strategies
	const int n = U.size();

	//Find all k-nearest neighbors
	TVector2D <size_t> knn_id;
	TVector2D <float> knn_dist;
	TVector <int> K(1, k);
	findAllKNN(U, U, K, knn_id, knn_dist);

	//Compute average normals, regression errors
	getAveragePointNormal(U, knn_id, AN);

	//Compute average normals, facility costs
	std::vector <double> FC;

	//Non-uniform version of IFL
	if (ifl_non_uniform)
	{
		if (fnorm != nL2)
			recomputeFacilityCostsFNFL(f, 0.0, 5.0, AN, fnorm, FC);
		else
			recomputeFacilityCostsFNFL(radius, f, 5.0, AN, nL2, FC);
	}

	//Uniform version of IFL
	else
		FC.insert(FC.end(), n, f);

	//Process all points
	std::cout << ">> Clusterization: ";
	TVector <int> IF(n);

	//First point is the facility
	Facility f0(U[0], 0, FC[U[0].id]);
	F.push_back(f0);

	for (int i = 1; i < n; i++)
	//for (int i = 0; i < 18; i++)       //hrensko_subset_small.txt
	{
		//Print
		if (i % (25000) == 0)
			std::cout << i << " ";

		//Update clusters using the heuristics
		updateClusters(U[i], radius, fnorm, F, AN, FC);
	}
}


void FFL(TVector <Point3D>& U, const double f, const double alpha, const double radius, const int k, const pfnorm& fnorm, TVector <Facility>& F, TVector2D <double>& AN, int& deleted)
{
	//Perform uniform facility cost clustering by FFL according to a given norm
	//Fast uniform facility location algorithm by Fotakis (non-deterministic)
	const int n = U.size();

	//Initialize random number generator
	std::mt19937_64 rd;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rd.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	//Find all k-nearest neighbors
	TVector2D <size_t> knn_id;
	TVector2D <float> knn_dist;
	TVector <int> K(1, k);
	findAllKNN(U, U, K, knn_id, knn_dist);

	//Compute average normals, regression error and ABN
	std::vector <double> FC;
	getAveragePointNormal(U, knn_id, AN);

	//Process all points
	std::cout << ">> Clusterization: ";
	TVector <int> IF(n);

	for (int i = 0; i < n; i++)
	{
		//Print
		if (i % (25000) == 0) std::cout << i << " ";

		//Find nearest facility
		int imin = -1;
		double gmin;
		getNearestF(U[i], F, AN, fnorm, radius, imin, gmin);

		//No facility exists
		if (imin == -1)
			gmin = f + 1;

		//Create new facility with the probability
		const double p = min(gmin / (alpha * f), 1.0);
		double r = unif(rd);
		if (r < p)
		{
			//Compute new facility radius
			const double m = min(gmin, alpha * f) / 6;

			//Create new facility
			Facility w(U[i], m, f);

			//Set node as a facility
			IF[i] = 1;

			auto it_z = F.begin();
			while (it_z != F.end())
			{
				//Factory is closer than threshold
				if (fnorm(it_z->u, w.u, AN[it_z->u.id], AN[w.u.id], radius) <= it_z->m)
				{
					//Reset facility flag
					IF[it_z->u.id] = 0;

					//std::cout << "   R:" << it_z->u.x << " " << it_z->u.y << " " << it_z->u.z << " " << it_z->m << '\n';

					//Remove factory
					it_z = F.erase(it_z);

					deleted++;
				}
				else
					++it_z;
			}

			//std::cout << i << " " << w.u.x << " " << w.u.y << " " << w.u.z << " " << w.m << '\n';

			//Add facility to the list
			F.push_back(w);
		}
	}

	//Assign the nearest facility according to a pseudonorm
	if (fnorm != nL2)
		assignPointsToFacilities(U, IF, AN, radius, fnorm, F);

	//Assign the nearest facility according to L2 norm  (accelerated)
	else
		assignPointsToFacilitiesL2(U, IF, F);
}


void FNFL(TVector <Point3D>& U, const double f, const double alpha, const double lambda, const double radius, const int k, const pfnorm& fnorm, TVector <Facility>& F, TVector2D <double>& AN, int &deleted)
{
	//Perform non-uniform facility cost clustering by FNFL according to a given norm
	//Fast non-uniform facility location algorithm by Fotakis (non-deterministic)
	const int n = U.size();

	//Initialize random number generator
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> disu(0.0, 1.0);

	//Find all k-nearest neighbors
	TVector2D <size_t> knn_id;
	TVector2D <float> knn_dist;
	TVector <int> K(1, k);
	findAllKNN(U, U, K, knn_id, knn_dist);

	//Compute average normals, facility cost, regression error and ABN
	std::vector <double> FC;
	getAveragePointNormal(U, knn_id, AN);
	if (fnorm != nL2)
		recomputeFacilityCostsFNFL(f, 0.0, 5.0, AN, fnorm, FC);
	else
		recomputeFacilityCostsFNFL(radius, f, 5.0, AN, nL2, FC);

	//Process all points
	std::cout << ">> Clusterization: ";
	TVector <int> IF(n);
	for (int i = 0; i < n; i++)
	{
		//Print
		if (i % (25000) == 0) std::cout << i << " ";

		//Find nearest facility
		int imin_fu = -1;
		double gmin_fu;
		getNearestF(U[i], F, AN, fnorm, radius, imin_fu, gmin_fu);

		//No facility exists
		if (imin_fu == -1)
			gmin_fu = FC[i];

		//Compute radius of the ball Bu(rb)
		const double rb = gmin_fu / (lambda - 1);

		//Find cheapest alternative in the neighborhood of u (point inside the ball Bu)
		int imin_p = -1;
		double gmin_p;
		getNearestPoint(U[i], U, AN, FC, IF, rb, i, fnorm, radius, imin_p, gmin_p);

		//Create new facility with the probability p
		const double p = std::min(gmin_fu / (alpha * FC[imin_p]), 1.0);
		if (disu(gen) < p)
		{
			//Get nearest facility
			int imin_fw = -1;
			double gmin_fw;

			//No cheaper point inside the ball of rb radius
			//The nearest facility did not change
			if (imin_p == i)
			{
				imin_fw = imin_fu;
				gmin_fw = gmin_fu;
			}

			//We found a cheaper point inside the ball
			else
				getNearestF(U[imin_p], F, AN, fnorm, radius, imin_fw, gmin_fw);

			//Compute new ball radius
			const double f_uw = (imin_fw == -1 ? 0 : fnorm(U[i], U[imin_p], AN[U[i].id], AN[U[imin_p].id], radius));
			const double m = (imin_fw == -1 ? 0.5 * FC[imin_p] : min(gmin_fw, 3.0 * FC[imin_p] + (lambda - 3.0) / 2.0 * f_uw) / 6.0);

			//Create new facility
			Facility w(U[imin_p], m, FC[imin_p]);

			//Set node as a facility
			IF[imin_p] = 1;

			auto it_z = F.begin();
			while (it_z != F.end())
			{
				//Factory is closer than threshold
				if (fnorm(it_z->u, w.u, AN[it_z->u.id], AN[w.u.id], radius) <= it_z->m)
				{
					//Reset facility flag
					IF[it_z->u.id] = 0;

					//Remove factory
					it_z = F.erase(it_z);

					deleted++;
				}

				else
					++it_z;
			}

			//Add facility to the list
			F.push_back(w);
		}
	}

	//Assign the nearest facility according to a pseudonorm
	if (fnorm != nL2)
		assignPointsToFacilities(U, IF, AN, radius, fnorm, F);

	//Assign the nearest facility according to L2 norm  (accelerated)
	else
		assignPointsToFacilitiesL2(U, IF, F);
}



int main(int argc, char* argv[])
{
	//Load points
	TVector < Point3D> input_points, output_points;
	//TVector < Point3D>  resampled_points;

	loadPoints("rauenstein_jtsk.txt", input_points);
	//resamplePointCloud(input_points, 300, resampled_points);
	savePoints("rauenstein_jtsk_300.txt", resampled_points);


	std::cout << "*** HYBRID CLUSTERIZATION OF THE POINT CLOUD *** \n";

	//Parameters of the clusterization algorithm
	bool statistics = true, ifl_non_uniform = false, export_dxf = false;
	int knn = 50, ns = 200000;
	double f = 0.01, dc = 0.25, alpha_ffl = 19.0 / 8, alpha_fnfl = 38.0/7, lambda = 20;
	pfnorm fnorm = &nDIS;
	std::string file_name, fnorm_text = "dis", method_text = "";

	//Process command line parameters
	while (--argc > 0)
	{
		//Get - (A parameter follows)
		if (*argv[argc] == '-')
		{
			//Process parameter after -
			for (char* parameter = argv[argc]; ; )
			{
				switch (*(++parameter))
				{

					//Set statistics
					case 's':
					{
						statistics = true;
						break;
					}

					//Non-uniform version of IFL
					case 'n':
					{
						ifl_non_uniform = true;
						break;
					}


					//Non-uniform version of IFL
					case 'e':
					{
						export_dxf = true;
						break;
					}

					//Terminate character \0 of the argument
					case '\0':
						break;

					//Throw exception
					default:
						throw std::exception("Exception: Invalid parameter in command line!");
				}

				//Stop processing of the argument
				break;
			}
		}

		//Set values
		else if (*argv[argc] == '+')
		{
			//Get new command line parameter
			char* attribute = const_cast <char*> (argv[argc] + 1), * value = NULL;
			char* command = attribute;

			//Find splitter: =
			for (; (*command != '=') && (*command != '\0'); command++);

			//We found splitter, trim command and copy to value
			if ((*command == '=') && (*command != '\0'))
			{
				*command = '\0';
				value = command + 1;
			}

			//Throw exception
			if (attribute == NULL || value == NULL)
				throw std::exception("Exception: Invalid value in command line!");

			//Set clusterization norm
			if (!strcmp("norm", attribute))
			{
				//ABN
				if (!strcmp("dis", value))
				{
					fnorm = &nDIS;
					fnorm_text = "dis";
				}

				//arc
				else if (!strcmp("arc", value))
				{
					fnorm = &nARC;
					fnorm_text = "arc";
				}

				//DFP
				else if (!strcmp("dfp", value))
				{
					fnorm = &nDFP;
					fnorm_text = "dfp";
				}

				//ABLP
				else if (!strcmp("arc2", value))
				{
					fnorm = &nARC2;
					fnorm_text = "arc2";
				}

				//L2
				else if (!strcmp("l2", value))
				{
					fnorm = &nL2;
					fnorm_text = "l2";
				}

				else 
					throw std::exception("Exception: Invalid clusterization norm type in command line!");
			}

			//Set method
			else if (!strcmp("met", attribute))
			{
				//IFL
				if (!strcmp("ifl", value))
					method_text = "ifl";

				//FFL
				else if (!strcmp("ffl", value))
					method_text = "ffl";

				//FNFL
				else if (!strcmp("fnfl", value))
					method_text = "fnfl";

				else
					throw std::exception("Exception: Invalid method type in command line!");
			}

			//Set the facility cost
			else if (!strcmp("f", attribute))
			{
				f = std::max(std::min(atof(value), 1000.0), 0.0001);
			}

			//Set max cluster size
			else if (!strcmp("dc", attribute))
			{
				dc = std::max(std::min(atof(value), 10.0), 0.01);
			}

			//Split cloud to subsets
			else if (!strcmp("ns", attribute))
			{
				ns = std::max(std::min(atoi(value), 500000), 100);
			}

			//Amount of knn
			else if (!strcmp("knn", attribute))
			{
				knn = std::max(std::min(atoi(value), 500), 5);
			}

			//Set isotropic ratio
			else if (!strcmp("mju", attribute))
			{
				mju = std::max(std::min(atof(value), 1.0), 0.0);
			}

			//Set order
			else if (!strcmp("l", attribute))
			{
				l = std::max(std::min(atof(value), 3.0), 0.5);
			}

			//Bad argument
			else
			{
				//std::cout << attribute << '\n';
				throw std::exception("Exception: Invalid attribute in command line!");
			}
		}

		//Process file
		else
		{
			file_name = argv[argc];
		}
	}

	method_text = "ffl";
	statistics = true;
	
	if (method_text == "fnfl" || method_text == "ffl")
	{
		fnorm = &nDFP;
		fnorm_text = "dfp";
	}

	else
	{
		fnorm = &nL2;
		fnorm_text = "l2";
		fnorm = &nABN;
		fnorm_text = "abn";
	}
	
	
	//List of files
	std::string file_list = file_name + ".list";

	//Print parameters
	std::cout << "\nParameters: norm = " << fnorm_text << ", method = " << method_text << ", f_cost = " << f << ", radius = " << dc << ", l = " << l << ", mju = " << mju << ", max_clust = " << dc << ", knn = " << knn << ", n_split = " << ns << ".\n\n";

	//Load list of kd subsets if they exist
	TVector<std::string> file_name_subsets;
	loadKDSubsets(file_list, file_name_subsets);

	//Otherwise, create list of kd subsets
	if (file_name_subsets.size() == 0)
	{
		createKDSubsets(file_name, ns, file_name_subsets);
		saveKDSubsets(file_list, file_name_subsets);
	}

	//KD subsets has been loaded
	else
	{
		std::cout << ">> KD subsets loaded (" << file_name_subsets.size() << ") ...";
	}

	//Supplementary variables
	TVector <int> n_points, n_deleted,  n_facilities;
	TVector <double> time_diff, n_clust_points_min, n_clust_points_max, n_clust_points_aver, clust_radius_min, clust_radius_max, clust_radius_aver, 
		fac_cost_min, fac_cost_max, fac_cost_aver, fac_abn_diff_min, fac_abn_diff_max, fac_abn_diff_aver, fac_dfp_diff_min, fac_dfp_diff_max, fac_dfp_diff_aver;

	//Process subsets one by one
	unsigned int i = 0, n_subsets = file_name_subsets.size();
	for (auto f_name : file_name_subsets)
	{
		//Load KD-subsets points
		TVector <Point3D> kd_subset, output_points_subset;
		loadPoints(f_name, kd_subset);

		std::cout << "\nFile: " << f_name << " (" << ++i << "/" << n_subsets << ", "  << kd_subset.size() << " points)" << '\n';

		//Output files
		std::string facil_file_subset = f_name +"_" + fnorm_text + "_l" + std::to_string(l) + "_" + method_text + "_f" + std::to_string(f) + "_dc" + std::to_string(dc) + "_facil.txt";
		std::string dxf_file_subset = f_name + "_" + fnorm_text + "_l" + std::to_string(l) + "_" + method_text + "_f" + std::to_string(f) + "_dc" + std::to_string(dc) + ".dxf";
		std::string result_file_subset = f_name + "_" + fnorm_text + "_l" + std::to_string(l) + "_" + method_text + "_f" + std::to_string(f) + "_dc" + std::to_string(dc) + "_results.log";

		//Output file does not exist
		if (!filesystem::exists(facil_file_subset))
		{
			//Apply facility location clusterization
			TVector <Facility> F;
			TVector2D <double> AN;

			//Auxilliary variables
			int deleted = 0;
			const clock_t begin_time = clock();

			//Incremental facility location (IFL)
			if (method_text.compare("ifl") == 0)
				IFL(kd_subset, f, alpha_ffl, dc, knn, fnorm, ifl_non_uniform, F, AN, deleted);

			//Uniform facility location, stochastic (FFL)
			else if (method_text.compare("ffl") == 0)
				FFL(kd_subset, f, alpha_ffl, dc, knn, fnorm, F, AN, deleted);

			//Non-uniform facility location, stochastic (FNFL)
			else
				FNFL(kd_subset, f , alpha_fnfl, lambda, dc, knn, fnorm, F, AN, deleted);

			//Statistics
			double time = (clock() - begin_time) / (CLOCKS_PER_SEC) / 1000.0;
			std::cout << "(time: " << time << "s, ";
			std::cout << F.size() << " facilities, ";
			std::cout << "deleted: " << deleted << "). ";
			std::cout << "OK \n\n";

			//Export facilities to the list
			for (Facility f : F)
				output_points_subset.push_back(f.u);

			//Compute local statistics
			if (statistics)
			{
				std::cout << ">> Computing statistics: \n";
				computeLocalStatistics(kd_subset, output_points_subset, f, knn, dc, fnorm, deleted, time, AN, F, result_file_subset);
			}

			//Export clusters to DXF
			if (export_dxf)
				DXFExport::exportClustersToDXF<double>(dxf_file_subset, F, AN);

			//Save facilities
			savePoints(facil_file_subset, output_points_subset);

			//Compute total cost
			const double total_cost = computeTotalCost(F, fnorm, AN, lambda);
			std::cout << "Total cost: " << total_cost << '\n';
		}
		
		//Output file with facilities exists, load facilities
		else
		{
			loadPoints(facil_file_subset, output_points_subset);
		}
		
		//Add all output facilities to the list
		output_points.insert(output_points.end(), output_points_subset.begin(), output_points_subset.end());
	}

	//Compute global statistics
	if (statistics)
	{
		std::string log_file = file_name + "_" + fnorm_text + "_l" + std::to_string(l) + "_" + method_text + "_f" + std::to_string(f) + "_dc" + std::to_string(dc) + ".res";
		computeGlobalStatistics(file_name_subsets, fnorm_text, method_text, f, dc, log_file);
	}

	//Save output points
	std::string facil_file = file_name + "_" + fnorm_text + "_l" + std::to_string(l) + "_" + method_text + "_f" + std::to_string(f) + "_dc" + std::to_string(dc) + "_facil_all.txt";
	savePoints(facil_file, output_points);

	std::cout << "Finished... \n";

	return 0;
}
