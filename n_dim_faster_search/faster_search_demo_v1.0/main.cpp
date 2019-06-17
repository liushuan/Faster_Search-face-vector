#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include <map>

#include "faster_search.h"

using namespace std;

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
std::vector<float *> get_random_data(const unsigned int number, const unsigned int dim) {
	std::vector<float *> random_data;
	for (size_t j = 0; j < number; j++)
	{
		float * data = new float[dim];
		for (std::size_t i = 0; i < dim; ++i) {
			data[i] = rand() * std::pow(-1, rand() % 2) / 100;
		}
		//normalize(data, dim);
		random_data.push_back(data);
	}
	return random_data;
}
std::vector<unsigned int> get_random_id(const unsigned int number) {
	std::vector<unsigned int> random_ids;
	for (size_t i = 0; i < number; i++)
	{
		random_ids.push_back(rand() % number);
	}
	return random_ids;
}

int main()
{
	srand(time(NULL));
	int maxLength = 20000;
	int person_size = 10000;
	int dim = 512;

	//generate data.
	std::vector<unsigned int> random_id = get_random_id(person_size);
	std::vector<float *>random_data = get_random_data(person_size, dim);
	std::vector<float *>random_test_data = get_random_data(1000, dim);

	//use faster search.
	float threhold = 0.3f; //Whether it's the same person
	Faster_Search faster_search(dim, threhold, maxLength);
	faster_search.load_data(random_data, random_id);
	faster_search.build_graph();
	//search person_id .
	unsigned int result_face_id = faster_search.search(random_test_data[0]);
	std::cout << "find person id is:" << result_face_id << std::endl;

	//faster_search.test_graph_accuracy(random_test_data);

	//动态添加人
	std::vector<float *>random_add_data = get_random_data(1, dim);
	unsigned int new_data_id = 10;
	bool is_success = faster_search.add_one_person(random_add_data[0], new_data_id);

	//动态删除人
	int remove_num = faster_search.remove_one_person(new_data_id);

	faster_search.test_graph_accuracy(random_test_data);

	system("pause");

	return 0;
}