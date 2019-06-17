#include <vector>
#include <iostream>
#include <fstream>
#include <string>

//#include <cblas.h>
#include <time.h>
#include "normal_search.h"

using namespace std;

void normal_test() {
	int dim = 512;
	std::vector<std::string> random_id = normal_search::get_lfw_id(12000);
	std::vector<float *> feature_data_test;
	std::vector<string> data_all = normal_search::readStringFromFileData("Alignment_LFW_Equalizedlfw_1.feature");
	std::vector<float *> feature_data = normal_search::get_lfw_data(feature_data_test, data_all, dim, 12000, 1000);
	std::cout << feature_data.size() << ":" << feature_data_test.size() << std::endl;

	normal_search::Face_Info face_info;
	face_info.datas = feature_data;
	face_info.dim = dim;
	face_info.face_ids = random_id;

	string id_face;

	clock_t start, finish;
	start = clock();
	int length = feature_data_test.size();
	for (size_t i = 0; i < length; i++)
	{
		id_face = normal_search::get_near_id(face_info, feature_data_test[i], 0);
	}
	finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;

	std::cout << "normal search time is:" << duration << " ms" << std::endl;

	std::cout << "id_face: " << id_face << std::endl;
}

//float GetDistance(float * dp1, float * dp2, int dimNum)
//{
//	//float ip = cblas_sdot(dimNum, dp1, 1, dp2, 1);
//	//float l2p1 = cblas_sdot(dimNum, dp1, 1, dp1, 1);
//	//float l2p2 = cblas_sdot(dimNum, dp2, 1, dp2, 1);
//	//return ip / sqrt(l2p1 * l2p2);
//
//	return cblas_sdot(dimNum, dp1, 1, dp2, 1);
//}
//static std::string get_near_id_blas(normal_search::Face_Info data, float test[], float threhold) {
//	float max_cos = -1;
//	std::string face_id = "-1";
//	unsigned int length = data.face_ids.size();
//
//	for (unsigned int i = 0; i < length; ++i)
//	{
//		float dist_cos = GetDistance(data.datas[i], test, data.dim);
//		if (dist_cos > max_cos && dist_cos > threhold) {
//			max_cos = dist_cos;
//			face_id = data.face_ids[i];
//		}
//	}
//	return face_id;
//}
//
//void openblas_test() {
//	
//	int dim = 512;
//	std::vector<std::string> random_id = normal_search::get_lfw_id(12000);
//	std::vector<float *> feature_data_test;
//	std::vector<string> data_all = normal_search::readStringFromFileData("Alignment_LFW_Equalizedlfw_1.feature");
//	std::vector<float *> feature_data = normal_search::get_lfw_data(feature_data_test, data_all, dim, 12000, 1000);
//	std::cout << feature_data.size() << ":" << feature_data_test.size() << std::endl;
//
//
//	normal_search::Face_Info face_info;
//	face_info.datas = feature_data;
//	face_info.dim = dim;
//	face_info.face_ids = random_id;
//
//	string id_face;
//	clock_t start, finish;
//	start = clock();
//	int length = feature_data_test.size();
//	for (size_t i = 0; i < length; i++)
//	{
//		id_face = get_near_id_blas(face_info, feature_data_test[i], 0);
//	}
//	finish = clock();
//	double duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
//	std::cout << "openblas search time is:" << duration << " ms" << std::endl;
//	std::cout <<"id_face: "<<id_face << std::endl;
//}


#include "sse_neon_search.h"
static std::string get_near_id_sse(normal_search::Face_Info data, float test[], float threhold) {
	float max_cos = -1;
	std::string face_id = "-1";
	unsigned int length = data.face_ids.size();

	for (unsigned int i = 0; i < length; ++i)
	{
		float dist_cos = sse_neon_search::_cal_similarity_avx_neon(data.datas[i], test, data.dim);
		//float dist_cos = GetDistance(, , );
		if (dist_cos > max_cos && dist_cos > threhold) {
			max_cos = dist_cos;
			face_id = data.face_ids[i];
		}
	}
	return face_id;
}
void sse_test() {

	int dim = 512;
	std::vector<std::string> random_id = normal_search::get_lfw_id(12000);
	std::vector<float *> feature_data_test;
	std::vector<string> data_all = normal_search::readStringFromFileData("Alignment_LFW_Equalizedlfw_1.feature");
	std::vector<float *> feature_data = normal_search::get_lfw_data(feature_data_test, data_all, dim, 12000, 1000);
	std::cout << feature_data.size() << ":" << feature_data_test.size() << std::endl;

	normal_search::Face_Info face_info;
	face_info.datas = feature_data;
	face_info.dim = dim;
	face_info.face_ids = random_id;

	string id_face;
	clock_t start, finish;
	start = clock();
	int length = feature_data_test.size();
	for (size_t i = 0; i < length; i++)
	{
		id_face = get_near_id_sse(face_info, feature_data_test[i], 0);
	}
	finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
	std::cout << "sse search time is:" << duration << " ms" << std::endl;
	std::cout << "id_face: " << id_face << std::endl;
}


//#include "none_search.h"
//static std::string get_near_id_none(normal_search::Face_Info data, float test[], float threhold) {
//	float max_cos = -1;
//	std::string face_id = "-1";
//	unsigned int length = data.face_ids.size();
//
//	for (unsigned int i = 0; i < length; ++i)
//	{
//		float dist_cos = none_search::dot(data.datas[i], test, data.dim);
//		//float dist_cos = GetDistance(, , );
//		if (dist_cos > max_cos && dist_cos > threhold) {
//			max_cos = dist_cos;
//			face_id = data.face_ids[i];
//		}
//	}
//	return face_id;
//}
//void none_test() {
//
//	int dim = 512;
//	std::vector<std::string> random_id = normal_search::get_lfw_id(12000);
//	std::vector<float *> feature_data_test;
//	std::vector<string> data_all = normal_search::readStringFromFileData("Alignment_LFW_Equalizedlfw_1.feature");
//	std::vector<float *> feature_data = normal_search::get_lfw_data(feature_data_test, data_all, dim, 12000, 1000);
//	std::cout << feature_data.size() << ":" << feature_data_test.size() << std::endl;
//
//	normal_search::Face_Info face_info;
//	face_info.datas = feature_data;
//	face_info.dim = dim;
//	face_info.face_ids = random_id;
//
//	string id_face;
//	clock_t start, finish;
//	start = clock();
//	int length = feature_data_test.size();
//	for (size_t i = 0; i < length; i++)
//	{
//		id_face = get_near_id_none(face_info, feature_data_test[i], 0);
//	}
//	finish = clock();
//	double duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
//	std::cout << "none search time is:" << duration << " ms" << std::endl;
//	std::cout << "id_face: " << id_face << std::endl;
//}

int main()
{
	normal_test();
	//openblas_test();
	sse_test();
	//none_test();

	system("pause");


	return 0;

}