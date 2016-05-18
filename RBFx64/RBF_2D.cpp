#include "stdafx.h"
#include "RBF_2D.h"


RBF_2D::RBF_2D()
{
	
}


RBF_2D::~RBF_2D()
{
}

int RBF_2D::count_words(std::string input_text)
{
	int number_of_words = 1;
	for (int i = 0; i < input_text.length(); i++)
		if (input_text[i] == ' ')
			number_of_words++;
	return number_of_words;
}

void RBF_2D::Read2Dpoints(std::string filename)
{
	std::string line;
	std::ifstream pFile(filename.c_str());
	int rows, col;
	
	double x, y ,z;

	if (pFile.is_open()){

		//getline(pFile, line);
		//rows = atoi(line.c_str());
		while (!pFile.eof()){
			getline(pFile, line);
			// count of space in each line
			col = count_words(line);
			// if the line is x y height map z
			if (col == 3){

				std::stringstream ss(line);
				ss >> x;
				ss >> y;
				ss >> z;
				_2dpts.push_back(Eigen::Vector3f(x, y,  z));

			}
			
		}
		

	}
}




void RBF_2D::thinPlateRFB()
{
	int N = _2dpts.size();
	A_rbf.resize(N, N);
	b_rbf.resize(N,1);
	for (size_t i = 0; i < N; i++){

		Eigen::Vector2f Xi(_2dpts.at(i).x(), _2dpts.at(i).y());
		b_rbf(i, 0) = _2dpts.at(i).z();
		for (size_t j = 0; j < N; j++){
			if (i != j) {
				Eigen::Vector2f Xj(_2dpts.at(j).x(), _2dpts.at(j).y());
				double r = (Xi - Xj).norm();
				A_rbf(i, j) = r*r * log(r);
			}
			else {
				A_rbf(i, j) = 0;
			}
		}
	}
}


void RBF_2D::SolveRFB()
{
	boost::progress_timer t;  // start timing
	//boost::timer::auto_cpu_timer t;
	int N = _2dpts.size();

	Eigen::MatrixXf AInvrbf = A_rbf.inverse();
	x_rbf = AInvrbf * b_rbf;
	std::ofstream("solution.txt") <<  A_rbf *  x_rbf;


}


void RBF_2D::getMinMax(){

	int N = _2dpts.size();
	Eigen::Matrix2Xf min;
	Eigen::Matrix2Xf max;
	Eigen::Matrix2Xf minmax;
	minmax.resize(N, 2);
	float xmn = FLT_MAX, xmx = -FLT_MAX, ymn = FLT_MAX, ymx = -FLT_MAX;
	for (size_t i = 0; i < N; i++){

		minmax(i, 0) = _2dpts.at(i).x();
		minmax(i, 1) = _2dpts.at(i).y();
		if (minmax(i, 0) < xmn) xmn = minmax(i, 0);
		if (minmax(i, 0) > xmx) xmx = minmax(i, 0);
		if (minmax(i, 1) < ymn) ymn = minmax(i, 1);
		if (minmax(i, 1) > ymx) ymx = minmax(i, 1);
	}

	min_p_[0] = xmn;
	max_p_[0] = xmx;

	min_p_[1] = ymn;
	max_p_[1] = ymx;
	
	/*std::cout<<"min 2d: " <<min_p_ <<std::endl;
	std::cout << "max 2d: " << max_p_ << std::endl;*/

}

void RBF_2D::InterpolateRBF()
{

	//boost::timer::auto_cpu_timer t;
	int N = _2dgridSize;
	int M = _2dpts.size(); // this number of actual data points
	_2grid.resize(N*N);
	_2gridM.resize(N);
	std::ofstream outGrid("surface.csv");
	getMinMax();
	for (int x = 0; x < N; ++x){
		_2gridM[x].resize(N);
		for (int y = 0; y < N; ++y){
			Eigen::Vector2f point;
			point[0] = min_p_[0] + (max_p_[0] - min_p_[0]) * x / N;
			point[1] = min_p_[1] + (max_p_[1] - min_p_[1]) * y / N;
			double f = 0.0f;
			// do function val in a loop x * A
				for (int j = 0; j < M; ++j){
					Eigen::Vector2f Xj(_2dpts.at(j).x(), _2dpts.at(j).y());
					double r = (point - Xj).norm();
					f += (r*r * log(r)) * x_rbf(j);
				}
		
			_2grid[x * N + y] =Eigen::Vector3f (point[0], point[1],f);
			_2gridM[x][y] = Eigen::Vector3f(point[0], point[1], f);
			//std::cout << _2grid[x * N + y] << std::endl;
			outGrid << point[0] << ","<< point[1] << "," <<f <<std::endl;

		}
	}
	//
	//for (int j = 0; j < N; ++j){
	//	Eigen::Vector2d Xj(_2dpts.at(j).x(), _2dpts.at(j).y());
	//	double fi = 0;
	//	for (int i = 0; i < N; ++i) {
	//		Eigen::Vector2d pos = Eigen::Vector2d(i - N / 2, j - N / 2)*35.0 / N;
	//		double r = (pos - Xj).norm();
	//		fi +=  (r*r * log(r)) * x_rbf(i);  // x[0] * pos[0] * pos[0] + x[1] * pos[1] * pos[1] + x[2] * pos[0] * pos[1] + x[3] * pos[0] + x[4] * pos[1] + x[5];
	//		_2grid[j * N + i] = Eigen::Vector3d(pos[0], pos[1], fi);
	//	}
	//}
	/*for (int j = 0; j < N; ++j){

		for (int i = 0; i < N; ++i) {
			outGrid << _2gridM[j][i].x() << "," << _2gridM[j][i].y() << "," << _2gridM[j][i].z() << std::endl;
		}

	}*/

}
