#pragma once

#include "stdafx.h"

const int _2dgridSize = 60;
class RBF_2D
{
public:
	std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > _2dpts;
	std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > _2grid;
	std::vector<std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > >_2gridM;
	Eigen::MatrixXf A_rbf;
	Eigen::MatrixXf b_rbf;
	Eigen::MatrixXf x_rbf;
	Eigen::Vector2f min_p_;
	Eigen::Vector2f max_p_;
	void getMinMax();

	RBF_2D();
	~RBF_2D();
	int count_words(std::string input_text);
	void Read2Dpoints(std::string filename);
	void thinPlateRFB();
	void SolveRFB();
	void InterpolateRBF();
};

