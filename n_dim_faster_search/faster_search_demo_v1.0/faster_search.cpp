/*************************************************
Copyright:tievd.com
Author:shuan.liu
Date:2019-01-10
Description:face vector faster search.
**************************************************/
#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include "faster_search.h"

int my_lower_bound(float *array, unsigned int size, float key)
{
	int first = 0, last = size - 1;
	int middle, pos = 0;       //需要用pos记录第一个大于等于key的元素位置

	while (first < last)
	{
		middle = (first + last) / 2;
		if (array[middle] < key) {      //若中位数的值小于key的值，我们要在右边子序列中查找，这时候pos可能是右边子序列的第一个
			first = middle + 1;
			pos = first;
		}
		else {
			last = middle;           //若中位数的值大于等于key，我们要在左边子序列查找，但有可能middle处就是最终位置，所以我们不移动last,
			pos = last;              //而是让first不断逼近last。
		}
	}
	return pos;
}

int my_upper_bound(float * array, unsigned int size, float key)
{
	int first = 0, last = size - 1;
	int middle, pos = 0;

	while (first < last)
	{
		middle = (first + last) / 2;
		if (array[middle] > key) {     //当中位数大于key时，last不动，让first不断逼近last
			last = middle;
			pos = last;
		}
		else {
			first = middle + 1;     //当中位数小于等于key时，将first递增，并记录新的位置
			pos = first;
		}
	}
	return pos;
}

float getAngleDistance(float feature1[], float feature2[], unsigned int size) {
	float distanceCos = 0;
	float nor_feature1 = 0, nor_feature2 = 0;
	for (int j = 0; j < size; j++)
	{
		distanceCos += feature1[j] * feature2[j];
		nor_feature1 += pow(feature1[j], 2);
		nor_feature2 += pow(feature2[j], 2);
	}
	distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
	return distanceCos;
}

long long get_near_id(std::shared_ptr<Data_FasterSearch> data_graph, float test[], float threhold) {
	float max_cos = -1;
	long long index = -1;
	unsigned int length = data_graph->face_info.datas.size();
	for (unsigned int i = 0; i < length; ++i)
	{
		float dist_cos = getAngleDistance(data_graph->face_info.datas[i], test, data_graph->dim);
		if (dist_cos > max_cos && dist_cos > threhold) {
			max_cos = dist_cos;
			index = data_graph->face_info.face_ids[i];
		}
		//std::cout << dist_cos << " ";
	}
	//std::cout << "dist_cos:" << max_cos<< " "<< index << std::endl;
	return index;
}

void Faster_Search::load_data(std::vector<float *> datas, std::vector<unsigned int> face_ids) {
	data_graph = std::make_shared<Data_FasterSearch>(dim);
	int length = datas.size();
	for (unsigned int j = 0; j < length; j++)
	{
		data_graph->face_info.datas.push_back(datas[j]);
		data_graph->face_info.face_ids.push_back(face_ids[j]);
	}
}

void Faster_Search::build_graph() {
	std::vector<std::multimap<float, int>> all_index;
	int length = data_graph->face_info.datas.size();

	for (unsigned int i = 0; i < data_graph->dim; i++)
	{
		base[i] = 1;
		std::multimap<float, int> cos_index;
		for (unsigned int j = 0; j < length; j++)
		{
			float dis = getAngleDistance(data_graph->face_info.datas[j], base, data_graph->dim);
			cos_index.insert(std::make_pair(dis, j));
		}

		all_index.push_back(cos_index);
		base[i] = 0;
	}
	for (unsigned int i = 0; i < data_graph->dim; i++)
	{
		float * cos = new float[MaxLength];
		unsigned int * index = new unsigned int[MaxLength];
		int j = 0;
		for (auto element : all_index[i])
		{
			cos[j] = element.first;
			index[j] = element.second;
			j++;
		}
		data_graph->base_coss.push_back(cos);
		data_graph->base_indexs.push_back(index);
	}

}

//search some one .
long long Faster_Search::search(float *test) {
	unsigned int length = data_graph->face_info.datas.size();

	int *results = new int[data_graph->face_info.datas.size()];
	memset(results, 0, data_graph->face_info.datas.size() * sizeof(int));

	unsigned int start = 0, end = length;
	for (unsigned int i = 0; i < data_graph->dim; ++i)
	{
		base[i] = 1;
		float diss = getAngleDistance(base, test, data_graph->dim);
		if ((diss - threhold2) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(data_graph->base_coss[i], length, diss - threhold2);
		}

		if (diss + threhold2 > 1) {
			end = length;
		}
		else {
			end = my_upper_bound(data_graph->base_coss[i], length, diss + threhold2) + 1;
		}
		// 15ms
		unsigned int * cursor = data_graph->base_indexs[i] + start;
		for (unsigned int j = start; j < end; ++j)
		{
			results[*(cursor++)]++;
		}
		base[i] = 0;
		if ((i + 1) == search_use_dim) {
			break;
		}
	}
	float max_cos = -1;
	long long resu_dex = -1;
	//int num_recall = 0;
	//10ms
	for (unsigned int i = 0; i < length; ++i)
	{
		if (results[i] == search_use_dim) {
			float dist = getAngleDistance(data_graph->face_info.datas[i], test, data_graph->dim);
			if ((dist > max_cos) && (dist > threhold1)) {
				max_cos = dist;
				resu_dex = data_graph->face_info.face_ids[i];
				//resu_dex = i;
			}
			//num_recall++;
		}
	}
	delete[] results;
	//std::cout << "  " << num_recall;
	return resu_dex;
}

unsigned int get_insert_pos(float * arr, unsigned int size, float data) {
	int first = 0, last = size - 1;
	int middle, pos = 0;
	while (first < last)
	{
		middle = (first + last) / 2;
		if (arr[middle] < data) {
			first = middle + 1;
			pos = first;
		}
		else {
			last = middle;
			pos = last;
		}
	}
	return pos;
}
template <typename T>
void insert_array(T*arr, unsigned int size, T data, unsigned int pos) {
	for (unsigned int i = size; i > pos; i--)
	{
		arr[i] = arr[i - 1];
	}
	arr[pos] = data;
}
bool Faster_Search::add_one_person(float * add_data, unsigned int face_id) {
	unsigned int length = data_graph->face_info.datas.size();
	if (length == MaxLength) {
		return false;
	}
	//add data.
	data_graph->face_info.datas.push_back(add_data);
	data_graph->face_info.face_ids.push_back(face_id);

	for (unsigned int i = 0; i < data_graph->dim; i++)
	{
		base[i] = 1;
		//std::multimap<float, int> cos_index;

		float dis = getAngleDistance(add_data, base, data_graph->dim);
		//cos_index.insert(std::make_pair(dis, data_faster_search.face_info.face_ids[j]));
		unsigned int pos = get_insert_pos(data_graph->base_coss[i], length, dis);
		//int pos = insert_array_data<float>(data_faster_search.base_coss[i], length, dis);
		insert_array<float>(data_graph->base_coss[i], length, dis, pos);
		insert_array<unsigned int>(data_graph->base_indexs[i], length, length, pos);
		base[i] = 0;
	}
	return true;
}

unsigned int found_array_data(unsigned int *arr, unsigned int size, unsigned int position) {
	for (unsigned int i = 0; i < size; ++i)
	{
		if (arr[i] == position) {
			return i;
		}
	}
	std::cout << "count't found face_id . error " << std::endl;
}
void delete_array_data_cos(float * arr, unsigned int size, unsigned int index) {
	for (unsigned int i = index; i < size - 1; i++)
	{
		arr[i] = arr[i + 1];
	}
}

void delete_array_data_index(unsigned int * arr, unsigned int size, unsigned int index, unsigned int position) {
	unsigned int last = size - 1;
	for (unsigned int i = 0; i < size; i++)
	{
		if (i >= index && i < last) {
			arr[i] = arr[i + 1];
		}
		if (arr[i] > position) {
			arr[i]--;
		}

	}
}

int Faster_Search::remove_one_person(unsigned int face_id) {
	unsigned int length = data_graph->face_info.datas.size();
	int remove_num = 0;
	//find data and erase data.
	while (true) {
		auto iter = std::find(data_graph->face_info.face_ids.begin(), data_graph->face_info.face_ids.end(), face_id);
		if (iter != data_graph->face_info.face_ids.end()) {
			unsigned int nPosition = distance(data_graph->face_info.face_ids.begin(), iter);
			data_graph->face_info.face_ids.erase(iter);

			auto begin = data_graph->face_info.datas.begin();
			//std::cout << "-----delete nPosition:" << nPosition << std::endl;
			data_graph->face_info.datas.erase(begin + nPosition);

			for (unsigned int i = 0; i < data_graph->dim; i++)
			{
				unsigned int index = found_array_data(data_graph->base_indexs[i], length, nPosition);
				delete_array_data_cos(data_graph->base_coss[i], length, index);
				delete_array_data_index(data_graph->base_indexs[i], length, index, nPosition);
			}
			remove_num++;
		}
		else {
			//std::cout << "can't found face_id:" << face_id << std::endl;
			break;
		}
	}
	return remove_num;
}

float Faster_Search::test_graph_accuracy(std::vector<float *>test_data) {
	const int test_size = test_data.size();
	//Face_Info test_data;
	long long * prim_id1 = new long long[test_size];
	long long * prim_id2 = new long long[test_size];

	std::cout << "all people length is:" << data_graph->face_info.datas.size() << std::endl;

	clock_t startTime, endTime;
	startTime = clock();//计时开始
	for (unsigned int i = 0; i < test_size; i++)
	{
		prim_id1[i] = get_near_id(data_graph, test_data[i], threhold1);
		//std::cout << " " << prim_id1[i];
	}
	//std::cout << std::endl;
	endTime = clock();//计时结束
	std::cout << "1 data run time is: " << (double)(endTime - startTime) << "ms" << std::endl;
	startTime = clock();//计时开始
	for (unsigned int i = 0; i < test_size; i++)
	{
		prim_id2[i] = search(test_data[i]);
		//std::cout << " " << prim_id2[i];
	}
	//std::cout << std::endl;
	endTime = clock();//计时结束
	std::cout << "2 data run time is: " << (double)(endTime - startTime) << "ms" << std::endl;
	float accuracy = 0;
	for (unsigned int i = 0; i < test_size; i++)
	{
		if (prim_id1[i] == prim_id2[i]) {
			accuracy += 1;
		}
	}
	std::cout << "accuracy is:" << accuracy / test_size << std::endl;
	delete[] prim_id1;
	delete[] prim_id2;
	return accuracy / test_size;
}

Faster_Search::~Faster_Search() {
	int dim = data_graph->base_coss.size();
	for (unsigned int i = 0; i < dim; i++)
	{
		delete[] data_graph->base_coss[i];
		delete[] data_graph->base_indexs[i];
	}
	int length = data_graph->face_info.datas.size();
	for (unsigned int i = 0; i < length; i++)
	{
		delete[] data_graph->face_info.datas[i];
	}
	delete[] base;
}