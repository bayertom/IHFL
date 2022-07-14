1. IHFL, short manual

Facility location clustering of the point cloud according to the hybrid constrained pseudonorm with additional penalty. 


1.1. Running the software
Open the command prompt and use the following combination of parameters, their value and options: 

ihfl file_name +parameter1=value1 +parameter2=value2 -option1 -option2

1.1 Input file
The input txt file contains the Cartesian coordinates X, Y, Z of points of the input point cloud. Each item is separated by space or TAB; see examples of input files.
Its nameis specified inside the quotes "".

Example:

	ihfl "test.txt"

1.2 List of parameters

1.2.1 Setting the pseudonorm
The pseudonorm can be set using the parameter "norm"

+norm=val

	abn: tangent model, pseudonorm G1
	dis: tangent model, pseudonorm G2
	ablp: secant model, pseudonorm G3
	dfp: secant model, pseudonorm G4
	l2: L2 norm
	
Example: Clusterization according to the dfp pseudonorm

	ihfl "test.txt" +norm=dfp
	
1.2.3 Setting the pseudonorm threshold
User-defined pseudonorm threshold refering to the maximum surface complexity (a maximum acceptable notch or protrusion) 
can be set using the parameter "f"

	+f=val

Typical value used for point clouds acquired by ALS is 1-5cm.

Example: Clusterization according to dfp pseudonorm using IHFL algorithm with the maximum surface complexity of 2 cm

	ihfl "test.txt" +norm=dfp +f=0.02

1.2.4 Setting the maximum ball radius
User-defined maximum value of the ball radius lambda can be set using the parameter "ball". This value represents the maximum
radius of the cluster.

	+lambda=value

Typical value used for point clouds acquired by ALS is 20-70 cm.

Example: Clusterization according to dfp pseudonorm using IHFL algorithm with the maximum surface complexity of 2 cm
and maximum ball radius of 50 cm

	ihfl "test.txt" +norm=dfp +f=0.02 +lambda=0.5

1.2.5 Setting the subset size
The input datasets can be recursively partitioned into subsets using kD-tree. The maximum amount points per a subset can be set
using the parameter "ns"

	+ns=value

Typical size of the subset used for point clouds acquired by ALS is 100000.

Example: Clusterization according to dfp pseudonorm using IHFL algorithm with the maximum surface complexity of 2 cm
and maximum ball radius of 50 cm. The point cloud is partitioned into subsets with the maximum size of 100 000 points.

	ihfl "test.txt" +norm=dfp +f=0.02 +lambda=0.5 +ns=100000

1.2.6 Setting the subset size
User-defined value of the k-nearest neighbors used for the estimation of the normal using PCA can be set using the parameter "knn"

	+knn=value

Typical amount of k-nearest neighbors for point clouds acquired by ALS is 30.

Example: Clusterization according to dfp pseudonorm using IHFL algorithm with the maximum surface complexity of 2 cm
and maximum ball radius of 50 cm. The point cloud is partitioned into subsets with the maximum size of 100 000 points, the normal
vector is estimated from 30 k-nearest neighbors.

	ihfl "test.txt" +norm=dfp +f=0.02 +lambda=0.5 +ns=100000 +knn=30

1.2.7 Setting the isotropic ratio
User defined isotropic factor mju, mju in (0,1), regulating the influence of the L2 metric and pseudometric.  Important parameter of 
the clusterization, significantly affects the behavior of the clusterization process: mju=0 -> L2 metric (fully isotropic), 
mju=1 ->pseudometric (fully anisotropic)

	+mju=value

Typical value of the isotropic factor is mju=0.95

Example: Clusterization according to dfp pseudonorm using IHFL algorithm with the maximum surface complexity of 2 cm
and maximum ball radius of 50 cm. The point cloud is partitioned into subsets with the maximum size of 100 000 points, the normal
vector is estimated from 30 k-nearest neighbors, the isotropic factor is set to 0.95.

	ihfl "test.txt" +norm=dfp +f=0.02 +lambda=0.5 +ns=100000 +knn=30 +mju=0.95

1.3 List of switches

1.3.1 Recompute values of facility costs
The costs of input points can be recomputed according to the ehavior of normal vectors using the switch "n"

	-n 	

Otherwise, the loaded or default costs are used.

Example: Clusterization according to dfp pseudonorm using IHFL algorithm with the maximum surface complexity of 2 cm
and maximum ball radius of 50 cm. The point cloud is partitioned into subsets with the maximum size of 100 000 points, the normal
vector is estimated from 30 k-nearest neighbors, the isotropic factor is set to 0.95.

	ihfl test.txt +norm=dfp +f=0.02 +lambda=0.5 +ns=100000 +knn=30 +mju=0.95 -n

1.3.3 Exporting clusters to DXF
The resulted facilities and connected clients can be exported into DXF file using the switch "e".

	-e

Example: non-uniform clusterization according to dfp pseudonorm using IHFL algorithm with the maximum surface complexity of 2 cm
and maximum ball radius of 50 cm. The point cloud is partitioned into subsets with the maximum size of 100 000 points, the normal
vector is estimated from 30 k-nearest neighbors, the isotropic factor is set to 0.95, the resulted clusters are exported into DXF file.

	ihfl test.txt +norm=dfp +met=ihfl +f=0.02 +lam=0.5 +ns=100000 +knn=30 +mju=0.95 -n -e

This option reduces the performance of clustering! 


1.4 Results of the clusterization
The output facilities as well as the output statistics are stored into *.txt files. They can easily be imported into external SW tool, for example the Cloud Compare.
