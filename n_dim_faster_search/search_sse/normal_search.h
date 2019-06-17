#ifndef _NORMAL_SEARCH__
#define _NORMAL_SEARCH__

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
//#include <io.h>

using namespace std;

namespace normal_search {

	struct Face_Info {
		std::vector <float *> datas;
		std::vector <std::string> face_ids;
		int dim;
	};

	//string transform to int double float...
	template <class Type>
	static Type stringToNum(string& str)
	{
		istringstream iss(str);
		Type num;
		iss >> num;
		return num;
	}

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
			std::vector<string> fe = split_string(features[i], ",");
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

	static float getAngleDistance(float feature1[], float feature2[], unsigned int size) {
		float distanceCos = 0;
		//float nor_feature1 = 1, nor_feature2 = 1;
		for (int j = 0; j < size; j++)
		{
			distanceCos += feature1[j] * feature2[j];
		}
		//distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
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
	std::vector<float *> get_lfw_data(std::vector<float *> &data_test, std::vector<string>all_datas, int dim, int number_train, int number_test) {
		std::vector<float *> random_data;
		for (size_t j = 0; j < (number_train + number_test); j++)
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

	static std::string get_near_id(Face_Info data, float test[], float threhold) {
		float max_cos = -1;
		std::string face_id = "-1";
		unsigned int length = data.face_ids.size();

		for (unsigned int i = 0; i < length; ++i)
		{
			float dist_cos = getAngleDistance(data.datas[i], test, data.dim);
			if (dist_cos > max_cos && dist_cos > threhold) {
				max_cos = dist_cos;
				face_id = data.face_ids[i];
			}
		}
		return face_id;
	}

}

#endif