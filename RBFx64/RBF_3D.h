#pragma once


#include "stdafx.h"

const float displacement_eps = 1.0;
const int grid_res = 50;

enum RBF_TYPE
{
	Gaussian, Linear, ThinPlateSpline, Cubic
};



class RBF_3D
{
public:
	RBF_3D();
	~RBF_3D();
	// I Will go with float persion to save space 
	Eigen::MatrixXf A_rbf;
	Eigen::MatrixXf P_rbf;

	Eigen::MatrixXf DataExtended;
	Eigen::MatrixXf b_rbf;
	Eigen::MatrixXf x_rbf;
	float displacement_eps ;
	int grid_res ;

	RBF_TYPE rbf_type;
	Eigen::Vector3f min_p_;
	Eigen::Vector3f max_p_;
	void getMinMax( );
	int count_words(std::string input_text);
	void Read3Dpoints(std::string filename, bool withNormals);
	void ComputeRFB();
	void SolveRFB();
	std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > _3dpts;
	std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > _3dnormals;
	std::vector<Eigen::Vector4f, Eigen::aligned_allocator<Eigen::Vector4f> > grid4MC;



	void interpolateRBF();
};

