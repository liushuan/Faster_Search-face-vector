#include <iostream>
#include <vector>
#include <iostream>
#include <sstream>
#include <array>
#include <iomanip>
#include "k_means.hpp"
#include "k_means_util.hpp"
using namespace std;

const int DIM = 512;

//struct Data_Array {
//	unsigned int data_length;
//	unsigned int dim;
//	float * data;
//};

struct Data_KMeans {
	std::vector<std::array<float, DIM>> face_data;
	std::vector<int> face_id;
	std::tuple<std::vector<std::array<float, DIM>>, std::vector<uint32_t>> cluster_data;
};

void normalize(vector<float>& feature1, unsigned int dim) {
	float s = 0;
	for (size_t i = 0; i < dim; i++)
	{
		s += (feature1[i] * feature1[i]);
	}
	s = std::sqrt(s);
	for (size_t i = 0; i < dim; i++)
	{
		feature1[i] = feature1[i] / s;
	}
}
void normalize(float * feature1, unsigned int dim) {
	float s = 0;
	for (size_t i = 0; i < dim; i++)
	{
		s += (feature1[i] * feature1[i]);
	}
	s = std::sqrt(s);
	for (size_t i = 0; i < dim; i++)
	{
		feature1[i] = feature1[i] / s;
	}

}
std::vector<std::array<float, DIM>> get_train_data(const int number) {
	std::vector<std::array<float, DIM>> train_data;

	for (size_t j = 0; j < number; j++)
	{
		std::array<float, DIM> data;
		std::vector<float> vals(DIM);
		for (std::size_t i = 0; i < DIM; ++i) {
			vals[i] = rand() * std::pow(-1, rand() % 2);
		}
		normalize(vals, DIM);

		for (size_t i = 0; i < DIM; i++)
		{
			data[i] = vals[i];
		}
		train_data.push_back(data);
	}
	return train_data;
 }

std::vector<std::array<float, DIM>> get_test_data(const int number) {
	std::vector<std::array<float, DIM>> test_data;

	for (size_t j = 0; j < number; j++)
	{
		std::array<float, DIM> data;
		std::vector<float> vals(DIM);
		for (std::size_t i = 0; i < DIM; ++i) {
			vals[i] = rand() * std::pow(-1, rand() % 2);
		}
		normalize(vals, DIM);

		for (size_t i = 0; i < DIM; i++)
		{
			data[i] = vals[i];
		}
		test_data.push_back(data);
	}
	return test_data;
}

auto build_graph(std::vector<std::array<float, DIM>> train_data, int num_cluser) {
	auto cluster_data = k_means::kmeans_lloyd(train_data, num_cluser);
	return cluster_data;
}

float getAngleDistance(std::array<float, DIM> feature1, std::array<float, DIM> feature2) {
	int length = DIM;
	float distanceCos = 0;
	float nor_feature1 = 1, nor_feature2 = 1;
	for (int j = 0; j < length; j++)
	{
		distanceCos += feature1[j] * feature2[j];
		//nor_feature1 += pow(feature1[j], 2);
		//nor_feature2 += pow(feature2[j], 2);
	}
	distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
	return distanceCos;

	//for (size_t i = 0; i < length; i++)
	//{
	//	distanceCos += std::pow((feature1[i] - feature2[i]), 2);	 
	//}
	//return std::sqrt(distanceCos);
}

int get_near_list(std::vector<std::array<float, DIM>> train_data, std::array<float, DIM> test, float threod) {
	float max_cos = -1;
	int index = -1;
	for (size_t i = 0; i < train_data.size(); i++)
	{
		float dist_cos = getAngleDistance(train_data[i], test);
		if (dist_cos > max_cos && dist_cos > threod) {
			max_cos = dist_cos;
			index = i;
		}
		//std::cout << dist_cos << " ";
	}
	return index;
}

void add_data(Data_KMeans &data_kmeans, std::array<float, DIM> &item, const int id) {
	auto means = std::get<0>(data_kmeans.cluster_data);
	auto labels = std::get<1>(data_kmeans.cluster_data);

	int cluter = -1;
	float max_val = -1;
	int j = 0;
	for (const auto& mean : means) {
		float dist = getAngleDistance(mean, item);
		if (dist > max_val) {
			max_val = dist;
			cluter = j;
		}
		j++;
	}
	//updata means
	for (size_t i = 0; i < DIM; i++)
	{
		means[cluter][i] = (means[cluter][i] * means.size() + item[i]) / (means.size() + 1);
	}
	//add label
	labels.push_back(cluter);
	data_kmeans.cluster_data = std::tuple<std::vector<std::array<float, DIM>>, std::vector<uint32_t>>(means, labels);
	data_kmeans.face_data.push_back(item);
	data_kmeans.face_id.push_back(id);
}

void remove_data(Data_KMeans &data_kmeans, const int id) {
	auto means = std::get<0>(data_kmeans.cluster_data);
	auto labels = std::get<1>(data_kmeans.cluster_data);
	// need get wright id. todo.................
	const int cluster_id = labels[id];
	//updata means
	for (size_t i = 0; i < DIM; i++)
	{
		means[cluster_id][i] = (means[cluster_id][i] * means.size() - data_kmeans.face_data[cluster_id][i]) / (means.size() - 1);
	}
	// need get wright id. todo.................
	labels.erase(labels.begin() + id);
	data_kmeans.face_data.erase(data_kmeans.face_data.begin() + id);
	data_kmeans.cluster_data = std::tuple<std::vector<std::array<float, DIM>>, std::vector<uint32_t>>(means, labels);
	data_kmeans.face_id.erase(data_kmeans.face_id.begin() + id);
}

void k_means_search(Data_KMeans &data_kmeans, const int length,
	std::vector<std::array<float, DIM>> &test_data,
	int *& prim_id2,
	const float threhold) {
	for (size_t i = 0; i < length; i++)
	{
		//int cluter = -1;
		std::vector<int> cluters;
		int current = 2;
		std::cout << "current :" ;
		for (int n = 0; n < 3; n++) {
			float max_val = -1;
			int j = 0;
			for (const auto& mean : std::get<0>(data_kmeans.cluster_data)) {
				float dist = getAngleDistance(mean, test_data[i]);
				if (dist > max_val && dist < current) {
					max_val = dist;
					cluters.push_back(j);
				}
				j++;
			}
			current = max_val;
			std::cout << " " << current;
		}

		int j = 0;
		int max_val = -1;
		int index = -1;
		for (const auto& label : std::get<1>(data_kmeans.cluster_data)) {
			if (label == cluters[0] || label == cluters[1] || label == cluters[2]) {
				float dist_cos = getAngleDistance(data_kmeans.face_data[j], test_data[i]);
				if (dist_cos > max_val && dist_cos > threhold) {
					max_val = dist_cos;
					index = j;
				}
			}
			j++;
		}
		prim_id2[i] = data_kmeans.face_id[index];

		if (i < 10) {
			std::cout << prim_id2[i] << " ";
		}
	}
}

static void  splitString(const string& s, vector<string>& v, const string& c)
{
	string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2 - pos1));
		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length()) {
		v.push_back(s.substr(pos1));
	}
}
std::vector<std::string> get_lfw_id(int number) {
	std::vector<string> random_ids;
	for (size_t i = 0; i < number; i++)
	{
		random_ids.push_back(to_string(i));
	}
	return random_ids;
}
//string transform to int double float...
template <class Type>
Type stringToNum(const string& str)
{
	istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}
std::vector<std::array<float, DIM>> get_lfw_data(std::vector<std::array<float, DIM>> &data_test, std::vector<string>all_datas, int dim, int number_train, int number_test) {
	std::vector<std::array<float, DIM>> train_data;
	for (size_t j = 0; j < (number_train+number_test); j++)
	{
		std::vector<string> str_d;
		splitString(all_datas[j], str_d, ",");

		if (j < number_train) {
			std::array<float, DIM> data;
			std::vector<float> vals(DIM);
			for (std::size_t i = 0; i < DIM; ++i) {
				vals[i] = stringToNum<float>(str_d[i]);
			}
			normalize(vals, DIM);
			for (size_t i = 0; i < DIM; i++)
			{
				data[i] = vals[i];
			}
			train_data.push_back(data);
		}
		else {
			std::array<float, DIM> data;
			std::vector<float> vals(DIM);
			for (std::size_t i = 0; i < DIM; ++i) {
				vals[i] = stringToNum<float>(str_d[i]);
			}
			normalize(vals, DIM);
			for (size_t i = 0; i < DIM; i++)
			{
				data[i] = vals[i];
			}
			data_test.push_back(data);
		}
	}
	return train_data;
}
#include <fstream>
// read file content  to vector.
std::vector<string> readStringFromFileData(string filePath)
{
	std::vector<string> data;
	std::ifstream fileA(filePath);
	if (!fileA)
	{
		return data;
	}
	for (int i = 0; !fileA.eof(); i++)
	{
		string buf;
		getline(fileA, buf, '\n');

		if (buf == "")
		{
			std::cout << "buf is empty." << std::endl;
			continue;
		}
		data.push_back(buf);
	}
	fileA.close();
	return data;
}
int main() {
	srand(time(NULL));
	float threhold = 0.1f;
	int cluter_k = 100;
	int train_length = 5000;
	int dim = 512;

	Data_KMeans data_kmeans;

	//data_kmeans.face_data = get_train_data(train_length);
	std::vector<string> data_all = readStringFromFileData("E:\\work_space\\source_code\\practise\\search\\faster_search_3.0\\Alignment_LFW_Equalizedlfw_1.feature");
	std::vector<std::array<float, DIM>> test_data;
	data_kmeans.face_data = get_lfw_data(test_data, data_all, dim, train_length, 100);
	std::cout << data_kmeans.face_data.size() << ":" << test_data.size() << std::endl;

	for (int i = 0; i < train_length; i++)
	{
		data_kmeans.face_id.push_back(i);
	}
	float sum = 0;
	for (size_t i = 0; i < dim; i++)
	{
		sum += test_data[4][i] * test_data[4][i];
	}
	std::cout << "sum is :" << sum << std::endl;

	clock_t startTime, endTime;
	startTime = clock();
	data_kmeans.cluster_data = build_graph(data_kmeans.face_data, cluter_k);

	//auto means = std::get<0>(data_kmeans.cluster_data);
	//auto labels = std::get<1>(data_kmeans.cluster_data);
	//for (auto e : labels)
	//{
	//	std::cout << e << " ";
	//}

	endTime = clock();
	std::cout << "build time is:" << (endTime - startTime) << "ms" << std::endl;

	//std::vector<std::array<float, DIM>> test_data = get_test_data(1000);
	const int  length = test_data.size();
	std::cout << "length :" << data_kmeans.face_data.size() << "  " << test_data.size() << std::endl;
	int * prim_id1 = new int[length];
	int * prim_id2 = new int[length];

	// 1. method.
	startTime = clock();
	for (size_t i = 0; i < length; i++)
	{
		prim_id1[i] = data_kmeans.face_id[get_near_list(data_kmeans.face_data, test_data[i], threhold)];
		if (i < 10) {
			std::cout << prim_id1[i] << " ";
		}
	}
	endTime = clock();
	std::cout << "1 all run time is: " << (double)(endTime - startTime) << "ms" << endl;

	/*std::cout << "1 index is:" << std::endl;
	for (size_t i = 0; i < length; i++)
	{
		std::cout << prim_id1[i] << " ";
	}
	std::cout << std::endl;*/

	//2. k_means
	

	/*std::cout << "2 index is:" << std::endl;
	for (size_t i = 0; i < length; i++)
	{
		std::cout << prim_id2[i] << " ";
	}
	std::cout << std::endl;*/

	startTime = clock();
	k_means_search(data_kmeans, length, test_data, prim_id2, threhold);
	endTime = clock();
	std::cout << "2 all run time is: " << (double)(endTime - startTime) << "ms" << endl;
	
	float accuracy = 0;
	for (size_t i = 0; i < length; i++)
	{
		if (prim_id1[i] == prim_id2[i]) {
			accuracy += 1;
		}
	}
	std::cout <<"accuracy is:"<<accuracy / length<< std::endl;


	//////////////////////////////////
	int add_num = 3000;
	std::vector<std::array<float, DIM>> need_add = get_test_data(add_num);
	startTime = clock();
	for (size_t i = 0; i < add_num; i++)
	{
		add_data(data_kmeans, need_add[i], train_length + i);
	}
	endTime = clock();
	std::cout << "add one item time is: " << (double)(endTime - startTime)/add_num << "ms" << endl;
	
	for (size_t i = 0; i < length; i++)
	{
		prim_id1[i] = data_kmeans.face_id[get_near_list(data_kmeans.face_data, test_data[i], threhold)];
	}
	k_means_search(data_kmeans, length, test_data, prim_id2, threhold);
	accuracy = 0;
	for (size_t i = 0; i < length; i++)
	{
		if (prim_id1[i] == prim_id2[i]) {
			accuracy += 1;
		}
	}
	std::cout << "add accuracy is:" << accuracy / length << std::endl;
	/////////////////////////////////

	/////////////////////////////////
	int remove_num = 3000;
	startTime = clock();
	for (size_t i = 0; i < add_num; i++)
	{
		remove_data(data_kmeans, rand() % (train_length - i));
	}
	endTime = clock();
	std::cout << "remove one item time is: " << (double)(endTime - startTime) / remove_num << "ms" << endl;

	for (size_t i = 0; i < length; i++)
	{
		prim_id1[i] = data_kmeans.face_id[get_near_list(data_kmeans.face_data, test_data[i], threhold)];
	}
	k_means_search(data_kmeans, length, test_data, prim_id2, threhold);
	accuracy = 0;
	for (size_t i = 0; i < length; i++)
	{
		if (prim_id1[i] == prim_id2[i]) {
			accuracy += 1;
		}
	}
	std::cout << "remove accuracy is:" << accuracy / length << std::endl;



	/////////////////////////////////




	system("pause");
}