// RBF.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "RBF_2D.h"
#include "RBF_3D.h"
int _tmain(int argc, _TCHAR* argv[])
{
    RBF_2D rbf2d;
	rbf2d.Read2Dpoints("g:\\spring2016\\amcs252\\project\\actualproject\\data\\2d\\seamount.csv");
	std::cout << rbf2d._2dpts.size() << std::endl;
	rbf2d.thinPlateRFB();
	std::ofstream("matrix.txt") << rbf2d.A_rbf;
	rbf2d.SolveRFB();
	rbf2d.InterpolateRBF();


	

	

	//RBF_3D rbf3d;
	//rbf3d.rbf_type = Cubic;
	//rbf3d.Read3Dpoints("G:\\Spring2016\\AMCS252\\Project\\ActualProject\\Data\\PointClouds\\sphere926.xyz",true);
	//std::cout << rbf3d._3dpts.size() << std::endl;
	////std::cout << rbf3d._3dnormals.size() << std::endl;

	//rbf3d.getMinMax();
	////std::cout << rbf3d.max_p_ << std::endl;
	////std::cout << rbf3d.min_p_ << std::endl;
	//rbf3d.ComputeRFB();
	////std::ofstream("DataExt_matrix3D.txt") << rbf3d.DataExtended;
	////std::ofstream("A_matrix3D.txt") << rbf3d.A_rbf;
	////std::ofstream("b_Vector3D.txt") << rbf3d.b_rbf;
	//
	//rbf3d.SolveRFB();
	//std::ofstream("x_Vector3D.txt") << rbf3d.x_rbf;
	////Eigen::MatrixXd cx(4, 1);
	//int N = rbf3d.DataExtended.size();
	////cx << rbf3d.x_rbf(N), rbf3d.x_rbf(N + 1), rbf3d.x_rbf(N + 2), rbf3d.x_rbf(N + 3);

	////std::cout <<  cx << std::endl;
	////rbf3d.interpolateRBF();
	system("PAUSE");
	return 0;
}

