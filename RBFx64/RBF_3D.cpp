#include "stdafx.h"
#include "RBF_3D.h"


RBF_3D::RBF_3D()
{
	 displacement_eps = 0.1;
	 grid_res = 10;
}


RBF_3D::~RBF_3D()
{
}
void RBF_3D::getMinMax(){

	int N = _3dpts.size();
	/*Eigen::Matrix3Xf min;
	Eigen::Matrix3Xf max;
	Eigen::Matrix3Xf minmax;
	minmax.resize(N, 3);*/
	
	float xmn = FLT_MAX, xmx = -FLT_MAX, ymn = FLT_MAX, ymx = -FLT_MAX, zmn = FLT_MAX, zmx = -FLT_MAX;
	for (size_t i = 0; i < N; i++){
		Eigen::Vector3f minmax;
		minmax(0) = _3dpts.at(i).x();
		minmax(1) = _3dpts.at(i).y();
		minmax(2) = _3dpts.at(i).z();
		if (minmax(0) < xmn) xmn = minmax( 0);
		if (minmax(0) > xmx) xmx = minmax( 0);
		if (minmax(1) < ymn) ymn = minmax( 1);
		if (minmax(1) > ymx) ymx = minmax( 1);
		if (minmax(2) < zmn) zmn = minmax( 2);
		if (minmax(2) > zmx) zmx = minmax( 2);
	}

	min_p_[0] = xmn;
	max_p_[0] = xmx;

	min_p_[1] = ymn;
	max_p_[1] = ymx;

	min_p_[2] = zmn;
	max_p_[2] = zmx;
	
}
	

int RBF_3D::count_words(std::string input_text)
{
	int number_of_words = 1;
	for (int i = 0; i < input_text.length(); i++)
		if (input_text[i] == ' ')
			number_of_words++;
	return number_of_words;
}
void RBF_3D::Read3Dpoints(std::string filename, bool withNormals)
{
	std::string line;
	std::ifstream pFile(filename.c_str());
	int rows, col;

	double x, y, z;
	double nx, ny, nz;

	if (pFile.is_open()){

		//getline(pFile, line);
		//rows = atoi(line.c_str());
		while (!pFile.eof()){
			getline(pFile, line);
			// count of space in each line
			col = count_words(line);
			// if the line is x y z
			if (/*col == 6 &&*/ withNormals== true){

				std::stringstream ss(line);
				ss >> x;
				ss >> y;
				ss >> z;
				ss >> nx;
				ss >> ny;
				ss >> nz;
				_3dpts.push_back(Eigen::Vector3f(x, y, z));//.normalized());
				_3dnormals.push_back(Eigen::Vector3f(nx, ny, nz));//.normalized());

			}
			if (/*col == 3 &&*/ withNormals == false){

				std::stringstream ss(line);
				ss >> x;
				ss >> y;
				ss >> z;

				_3dpts.push_back(Eigen::Vector3f(x, y, z));//.normalized());

			}

		}


	}
}
void RBF_3D::ComputeRFB(){
	int N = _3dpts.size();
	DataExtended.resize(2* N, 3);
	/*if (rbf_type == Cubic){
		b_rbf.resize((2 * N)+4, 1);
	}
	else
	{
		b_rbf.resize(2 * N, 1);
	}*/
	b_rbf.resize(2 * N, 1);
	// create extend data set +- normals
	
	   
		for (size_t i = 0; i < N; i++){

			DataExtended(i, 0) = _3dpts.at(i).x();
			DataExtended(i, 1) = _3dpts.at(i).y();
			DataExtended(i, 2) = _3dpts.at(i).z();
			b_rbf(i, 0) = 0;
		}
		size_t j= 0;
		for (size_t i = N; i < 2 * N; i++, j++){

			DataExtended(i, 0) = _3dpts[j][0] + displacement_eps * _3dnormals[j][0];
			DataExtended(i, 1) = _3dpts[j][1] + displacement_eps * _3dnormals[j][1];
			DataExtended(i, 2) = _3dpts[j][2] + displacement_eps *_3dnormals[j][2];
			b_rbf(i, 0) = displacement_eps;
		}
		/*j = 0;
		for (size_t i = 2*N; i < 3 * N; i++, j++){

			DataExtended(i, 0) = _3dpts[j][0] - displacement_eps * _3dnormals[j][0];
			DataExtended(i, 1) = _3dpts[j][1] - displacement_eps * _3dnormals[j][1];
			DataExtended(i, 2) = _3dpts[j][2] - displacement_eps *_3dnormals[j][2];
			b_rbf(i, 0) = -displacement_eps;

		}*/
		switch (rbf_type)
		{
		case Gaussian:
			A_rbf.resize(2 * N, 2 * N);
			break;
		case Linear:
		case Cubic:
			A_rbf.resize((2 * N) , (2 * N) + 4);
			//P_rbf.resize(2*N, 4);
			
			break;
	
		default:
			abort();
			break;
		}

		
		for (size_t i = 0; i < 2* N; i++){

			Eigen::Vector3f Xi(DataExtended(i, 0), DataExtended(i,1), DataExtended(i, 2));

			for (size_t j = 0; j <2 * N; j++){
				//if (i != j) {
					Eigen::Vector3f Xj(DataExtended(j, 0), DataExtended(j, 1), DataExtended(j, 2));
					float r = (Xi - Xj).norm();
					switch (rbf_type)
					{
					case Gaussian:
						r = expf(-0.5 * r*r);
						break;
					case Cubic:
						r = r* r * r;
						A_rbf(i,j+ 1) = 1.0;
						A_rbf(i, j+2) = Xi.x();
						A_rbf(i, j+3) = Xi.y();
						A_rbf(i, j+4) = Xi.z();
						
						break;
					case Linear:
					default:
						break;
					}
					A_rbf(i, j) = r;
					/*if (rbf_type == Cubic) {

						P_rbf(i, 0) = 1.0;
						P_rbf(i, 1) = Xi.x();
						P_rbf(i, 2) = Xi.y();
						P_rbf(i, 3) = Xi.z();
					}*/
				/*}
				else {
					A_rbf(i, j) = 0;
				}*/

			}
			/*if (rbf_type == Cubic){
			A_rbf.topRightCorner(2*N, 4) = P_rbf;
			A_rbf.bottomRightCorner(2*N, 4) = Eigen::MatrixXf::Zero(2*N, 4);
			A_rbf.bottomLeftCorner(4, 2*N) = P_rbf.transpose();
			b_rbf(2*N) = 0.0;
			b_rbf(2* N + 1) = 0.0;
			b_rbf(2*N+ 2) = 0.0;
			b_rbf(2*N+ 3) = 0.0;
			}*/


		}



	


		

}
void RBF_3D::SolveRFB()
{
	//boost::timer::auto_cpu_timer t;
	boost::progress_timer t;  // start timing
	
	
	
	x_rbf = A_rbf.fullPivLu().solve(b_rbf);
	//x_rbf = A_rbf.partialPivLu().solve(b_rbf);
	//x_rbf= A_rbf.lu().solve(b_rbf);
   //x_rbf = A_rbf.llt().solve(b_rbf); //does not work
	//x_rbf = A_rbf.ldlt().solve(b_rbf);
	//Eigen::MatrixXd AInvrbf = A_rbf.inverse();
	//x_rbf = AInvrbf * b_rbf;
	// itrative methods 
	//x_rbf = A_rbf.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_rbf);
	//x_rbf = A_rbf.colPivHouseholderQr().solve(b_rbf);
	//x_rbf = (A_rbf.transpose() * A_rbf).ldlt().solve(A_rbf.transpose() * b_rbf);
	double relative_error = (A_rbf*x_rbf - b_rbf).norm() / b_rbf.norm(); // norm() is L2 norm
	std::cout << "The relative error is:\n" << relative_error << std::endl;

}

void RBF_3D::interpolateRBF()
{
	int N = grid_res;
	int M = _3dpts.size(); // this number of actual data points
    grid4MC.resize(N*N*N);

	int res_x_ = grid_res, res_y_ = grid_res, res_z_ = grid_res;
	//very important step to intialize grid extents 
	getMinMax();
	std::ofstream outGrid("Volume.txt");
	for (int x = 0; x < res_x_; ++x)
		for (int y = 0; y < res_y_; ++y)
			for (int z = 0; z < res_z_; ++z)
			{
				Eigen::Vector3d point;
				point[0] = min_p_[0] + (max_p_[0] - min_p_[0]) * x / res_x_;
				point[1] = min_p_[1] + (max_p_[1] - min_p_[1]) * y / res_y_;
				point[2] = min_p_[2] + (max_p_[2] - min_p_[2]) * z / res_z_;

				float f =  0.0f;
				
				// do function val in a loop x * A
				for (int j = 0; j < M; ++j){
					Eigen::Vector3d Xj(_3dpts.at(j).x(), _3dpts.at(j).y(), _3dpts.at(j).z());
					Eigen::Vector3d XPlusNormals(_3dpts.at(j).x() + displacement_eps * _3dnormals[j][0], _3dpts.at(j).y() + displacement_eps * _3dnormals[j][1], _3dpts.at(j).z() + displacement_eps * _3dnormals[j][1]);
					Eigen::Vector3d XMinusNormals(_3dpts.at(j).x() - displacement_eps * _3dnormals[j][0], _3dpts.at(j).y()- displacement_eps * _3dnormals[j][1], _3dpts.at(j).z() - displacement_eps * _3dnormals[j][1]);
					float r1 = (point - Xj).norm();
					float r2 = (point - XPlusNormals).norm();
					float r3 = (point - XMinusNormals).norm();
					float px = 0.0;
					int k = DataExtended.size();
					switch (rbf_type)
					{
					case Gaussian:
						r1 = expf(-0.5 * r1*r1);
						r2 = expf(-0.5 * r2*r2);
						r3 = expf(-0.5 * r3*r3);
						break;
					case Cubic:
						r1 = r1* r1 * r1;
						r2 = r2* r2 * r2;
						r3 = r3* r3 * r3;
					
						px = Eigen::Vector4f(1, point[0], point[1], point[2]).dot(Eigen::Vector4f(x_rbf(k), x_rbf(k + 1), x_rbf(k + 2), x_rbf(k + 3)));
						break;
					case Linear:
					default:
						break;
					}
					f += r1* x_rbf(j) + px; //+ r2 * x_rbf(j + N) + r3 * x_rbf(j + 2*N);
				}

				Eigen::Vector4f v = Eigen::Vector4f(point[0], point[1], point[2], f);
				//grid4MC[x * res_y_*res_z_ + y * res_z_ + z] = Eigen::Vector4f(point[0], point[1], point[2], f);
				outGrid << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << std::endl;
			}


}

