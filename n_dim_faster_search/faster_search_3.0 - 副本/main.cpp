#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include <map>
#include <fstream>
#include <sstream>
#include "faster_search.h"

using namespace std;

// read file content  to vector.
static std::vector<string> readStringFromFileData(string filePath)
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

//string transform to int double float...
template <class Type>
static Type stringToNum(const string& str)
{
	istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}
// string split 
static vector<string> split_string(const string& s, const string& c)
{
	vector<string> v;
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
	return v;
}

std::vector<float *> get_feats_data(string file_name) {
	std::vector<float *> feature_data;
	std::vector<string> features = readStringFromFileData(file_name);
	int length = features.size();
	for (size_t i = 0; i < length; i++)
	{
		std::vector<string> fe= split_string(features[i], ",");
		int dim = fe.size();
		float * data = new float[dim];
		for (size_t j = 0; j < dim; j++)
		{
			data[j] = stringToNum<float>(fe[j]);
		}
		feature_data.push_back(data);
	}
	return feature_data;
}

static void normalize(float * feature1, unsigned int dim) {
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
std::vector<float *> get_random_data(const unsigned int number, const unsigned int dim) {
	std::vector<float *> random_data;
	for (size_t j = 0; j < number; j++)
	{
		float * data = new float[dim];
		for (std::size_t i = 0; i < dim; ++i) {
			data[i] = rand() * std::pow(-1, rand() % 2);
		}
		normalize(data, dim);
		random_data.push_back(data);
	}
	return random_data;
}
std::vector<string> get_random_id(const unsigned int number) {
	std::vector<string> random_ids;
	for (size_t i = 0; i < number; i++)
	{
		random_ids.push_back(to_string(i));
	}
	return random_ids;
}
static float getAngleDistance(float feature1[], float feature2[], unsigned int size) {
	float distanceCos = 0;
	float nor_feature1 = 1, nor_feature2 = 1;
	for (int j = 0; j < size; j++)
	{
		distanceCos += feature1[j] * feature2[j];
	}
	distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
	return distanceCos;
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
std::vector<float *> get_lfw_data(std::vector<float *> &data_test,std::vector<string>all_datas ,int dim, int number_train, int number_test) {
	std::vector<float *> random_data;
	for (size_t j = 0; j < (number_train+number_test); j++)
	{
		std::vector<string> str_d;
		splitString(all_datas[j], str_d, ",");
		if (j < number_train) {
			float * data = new float[dim];
			for (std::size_t i = 0; i < dim; ++i) {
				data[i] = stringToNum<float>(str_d[i]);
			}
			normalize(data, dim);
			random_data.push_back(data);
		}
		else {
			float * data = new float[dim];
			for (std::size_t i = 0; i < dim; ++i) {
				data[i] = stringToNum<float>(str_d[i]);
			}
			normalize(data, dim);
			data_test.push_back(data);
		}
	}
	return random_data;
}


int main()
{
	srand(time(NULL));
	int maxLength = 14000;
	int person_size = 12000;
	int dim = 512;
	
	//generate data.
	//std::vector<std::string> random_id = get_random_id(person_size);
	//std::vector<float *>feature_data = get_random_data(person_size, dim);
	//std::vector<float *>feature_data_test = get_random_data(1000, dim);
	
	std::vector<std::string> random_id = get_lfw_id(12000);
	std::vector<float *> feature_data_test;
	std::vector<string> data_all = readStringFromFileData("E:\\work_space\\source_code\\practise\\search\\faster_search_3.0\\Alignment_LFW_Equalizedlfw_1.feature");
	std::vector<float *> feature_data = get_lfw_data(feature_data_test, data_all, dim, 12000, 100);
	std::cout << feature_data.size() << ":" << feature_data_test.size() << std::endl;
	
	for (size_t i = 0; i < 10; i++)
	{
		std::cout<<"dis :"<<getAngleDistance(feature_data[i], feature_data[i+1], dim) << std::endl;
	}
	
	//use faster search.
	float threhold = 0.38f; //Whether it's the same person
	search::Faster_Search faster_search(dim, threhold, maxLength);
	faster_search.load_data(feature_data, random_id);
	
	clock_t startTime, endTime;
	startTime = clock();
	faster_search.build_graph();
	endTime = clock();
	std::cout << "build graph time is: " << (double)(endTime - startTime) << "ms" << std::endl;

	string result_face_id = faster_search.search(feature_data_test[1]);
	std::cout << "find person id is:" << result_face_id << std::endl;
	faster_search.test_graph_accuracy(feature_data_test);

//	faster_search.save_face_data("face_data.bin", feature_data_test);

//	faster_search.storage_graph("graph.bin");

	/*std::vector<float *>feature_data_test;
	//std::vector<float *>feature_data_test = get_random_data(1000, dim);
	//use faster search.
	float threhold = 0.2f; //Whether it's the same person
	Faster_Search faster_search(dim, threhold, maxLength);
	faster_search.load_face_data("face_data.bin", feature_data_test);
	
	//feature_data_test = get_random_data(1000, dim);

	faster_search.restorage_graph("graph.bin");
	string result_face_id = faster_search.search(feature_data_test[1]);
	std::cout << "find person id is:" << result_face_id << std::endl;
	faster_search.test_graph_accuracy(feature_data_test);
	*/


	/*
	//动态添加人
	std::vector<float *>random_add_data = get_random_data(1000, dim);
	for (size_t i = 0; i < 1000; i++)
	{
		string new_data_id = to_string(i);
		bool is_success = faster_search.add_one_person(random_add_data[i], new_data_id);
	}

	//动态删除人
	int size = 0;
	for (size_t i = 0; i < 500; i++)
	{
		string data_id = to_string(i);
		int remove_num = faster_search.remove_one_person(data_id);
		size += remove_num;
	}
	std::cout << "remove size is:" << size << std::endl;
	faster_search.test_graph_accuracy(random_test_data);
	*/
	system("pause");

	return 0;
}