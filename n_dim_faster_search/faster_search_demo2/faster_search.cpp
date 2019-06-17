/*************************************************
Copyright:tievd.com
Author:shuan.liu
Date:2019-01-10
Description:face vector faster search.
**************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <time.h>
#include <string>
#include "faster_search.h"

static int my_lower_bound(float *array, unsigned int size, float key)
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

static int my_upper_bound(float * array,unsigned int size, float key)
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

static std::string get_near_id(std::shared_ptr<Data_FasterSearch> data_graph, float test[], float threhold) {
	float max_cos = -1;
	std::string face_id = "-1";
	unsigned int length = data_graph->face_info.datas.size();
	for (unsigned int i = 0; i < length; ++i)
	{
		float dist_cos = getAngleDistance(data_graph->face_info.datas[i], test, data_graph->dim);
		if (dist_cos > max_cos && dist_cos > threhold) {
			max_cos = dist_cos;
			face_id = data_graph->face_info.face_ids[i];
		}
		//std::cout << dist_cos << " ";
	}
	//std::cout << "dist_cos:" << max_cos<< " "<< face_id << std::endl;
	return face_id;
}

void Faster_Search::load_data(std::vector<float *> datas, std::vector<std::string> face_ids) {
	int length = datas.size();
	for (unsigned int j = 0; j < length; j++)
	{
		data_graph->face_info.datas.push_back(datas[j]);
		data_graph->face_info.face_ids.push_back(face_ids[j]);
	}
}

template <typename T>
static float ave(T a[], const int n)
{
	float sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i];
	}
	return sum / n;
}
template <typename T>
static float variance(T a[], const int n)
{
	float sum = 0;
	float average = ave<T>(a, n);
	for (int i = 0; i < n; i++) {
		sum += (a[i] - average)*(a[i] - average);
	}
	return sum/n;
}
static void bubble_sort1(float a[], int b[], int n)
{
	int i, j;
	for (i = n - 1; i>0; i--)
	{
		// 将a[0...i]中最大的数据放在...
		for (j = 0; j<i; j++)
		{
			if (a[j] > a[j + 1]) {
				std::swap(a[j], a[j + 1]);
				std::swap(b[j], b[j + 1]);
			}
		}
	}

}

static void get_max_var_index(std::vector<float * >base_coss, int *& sort_var_index, const int length) {
	int split_size = 200; //[-1 ，1 ] split 100.
	int half_size = split_size / 2;
	//statistic analysis.
	int *static_num = new int[split_size];
	//memset(static_num, 0, split_size * sizeof(int));
	
	int dim = base_coss.size();
	float *index_var = new float[dim];

	sort_var_index = new int[dim];
	for (size_t i = 0; i < dim; i++)
	{
		sort_var_index[i] = i;
	}
	for (size_t i = 0; i < dim; i++)
	{
		memset(static_num, 0, split_size * sizeof(int));
		for (size_t j = 0; j < length; j++)
		{
			static_num[(int)std::floor(base_coss[i][j] * half_size) + half_size]++;
		}
		index_var[i] = variance<int>(static_num, split_size);
		//index_var[i] = base_coss[i][length - 1] - base_coss[i][0];
	}
	bubble_sort1(index_var, sort_var_index, dim);
}

void Faster_Search::build_graph() {
	std::vector<std::multimap<float, int>> all_index;
	int length = data_graph->face_info.datas.size();

	for (unsigned int i = 0; i < search_use_dim; i++)
	{
		base[i] = 1;
		std::multimap<float, int> cos_index;
		for (unsigned int j = 0; j < length; j++)
		{
			float dis = getAngleDistance(data_graph->face_info.datas[j], base, data_graph->dim);
			cos_index.insert(std::make_pair(dis, j));
			if (j < 10) {
				std::cout << dis << " ";
			}
		}

		all_index.push_back(cos_index);
		base[i] = 0;
	}
	for (unsigned int i = 0; i < search_use_dim; i++)
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
	get_max_var_index(data_graph->base_coss, data_graph->sort_var_index, length);

	//for (size_t i = 0; i < dim; i++)
	//{
	//	std::cout << " " << data_graph->sort_var_index[i];
	//}
}

//search some one .
std::string Faster_Search::search(float *test) {
	unsigned int length = data_graph->face_info.datas.size();
	//std::cout << "length is:" << length << std::endl;
	bool *results = new bool[length];
	memset(results, 1, length * sizeof(bool));

	unsigned int start = 0, end = length;
	for (unsigned int i = 0; i < search_use_dim; ++i)
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
		//std::cout << "start :end " << start << " : " << end << std::endl;
		// 15ms
		
		/*unsigned int * cursor = data_graph->base_indexs[i] + start;
		for (unsigned int j = start; j < end; ++j)
		{
			results[*(cursor++)]++;
		}*/

		unsigned int * cursor = data_graph->base_indexs[i];
		for (unsigned int j = 0; j < start; ++j)
		{
			results[*(cursor++)] = false;
		}
		unsigned int * cursor2 = data_graph->base_indexs[i] + end;
		for (unsigned int j = end; j < length; ++j)
		{
			results[*(cursor2++)] = false;
		}

		base[i] = 0;
		
		/*if (current_start > start) {
			start = current_start;
		}
		if (current_end < end) {
			end = current_end;
		}*/
		
		/*if ((i + 1) == search_use_dim) {
			break;
		}*/
	}
	float max_cos = -1;
	std::string face_id = "-1";
	int num_recall = 0;
	//10ms
	//std::cout << "start :end " << start << " : " << end<<std::endl;
	for (unsigned int i = 0; i < length; ++i)
	{
		if (results[i]) {
			//std::cout << i << std::endl;
			float dist = getAngleDistance(data_graph->face_info.datas[i], test, data_graph->dim);
			if ((dist > max_cos) && (dist > threhold1)) {
				max_cos = dist;
				face_id = data_graph->face_info.face_ids[i];
				//resu_dex = i;
			}
			num_recall++;
		}
	}
	delete[] results;
	std::cout << "  " << num_recall;
	return face_id;
}

std::string Faster_Search::search_sse_neon(float *test) {
	unsigned int length = data_graph->face_info.datas.size();

	float max_cos = -1;
	std::string face_id = "-1";

	for (unsigned int i = 0; i < length; ++i)
	{
		float dist_cos = sse_neon_search::_cal_similarity_avx_neon(data_graph->face_info.datas[i], test, data_graph->dim);
		//float dist_cos = GetDistance(, , );
		if (dist_cos > max_cos && dist_cos > threhold1) {
			max_cos = dist_cos;
			face_id = data_graph->face_info.face_ids[i];
		}
	}
	return face_id;
}



std::string Faster_Search::search2(float *test) {
	unsigned int length = data_graph->face_info.datas.size();

	int *results = new int[data_graph->face_info.datas.size()];
	memset(results, 0, data_graph->face_info.datas.size() * sizeof(int));

	unsigned int start = 0, end = length;
	for (unsigned int i = 0; i < search_use_dim; ++i)
	{
		base[data_graph->sort_var_index[i]] = 1;
		float diss = getAngleDistance(base, test, data_graph->dim);
		if ((diss - threhold2) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(data_graph->base_coss[data_graph->sort_var_index[i]], length, diss - threhold2);
		}

		if (diss + threhold2 > 1) {
			end = length;
		}
		else {
			end = my_upper_bound(data_graph->base_coss[data_graph->sort_var_index[i]], length, diss + threhold2) + 1;
		}
		// 15ms
		unsigned int * cursor = data_graph->base_indexs[data_graph->sort_var_index[i]] + start;
		for (unsigned int j = start; j < end; ++j)
		{
			results[*(cursor++)]++;
		}
		base[data_graph->sort_var_index[i]] = 0;
		
		/*if ((i + 1) == search_use_dim) {
			break;
		}*/
	}
	float max_cos = -1;
	std::string face_id = "-1";
	//int num_recall = 0;
	//10ms
	for (unsigned int i = 0; i < length; ++i)
	{
		if (results[i] == search_use_dim) {
			float dist = getAngleDistance(data_graph->face_info.datas[i], test, data_graph->dim);
			if ((dist > max_cos) && (dist > threhold1)) {
				max_cos = dist;
				face_id = data_graph->face_info.face_ids[i];
				//resu_dex = i;
			}

			//num_recall++;
		}
	}
	delete[] results;
	//std::cout << "  " << num_recall;
	return face_id;
}

static unsigned int get_insert_pos(float * arr, unsigned int size, float data) {
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
static void insert_array(T*arr, unsigned int size, T data, unsigned int pos) {
	for (unsigned int i = size; i > pos; i--)
	{
		arr[i] = arr[i - 1];
	}
	arr[pos] = data;
}
bool Faster_Search::add_one_person(float * add_data, std::string face_id) {
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

static unsigned int found_array_data(unsigned int *arr, unsigned int size, unsigned int position) {
	for (unsigned int i = 0; i < size; ++i)
	{
		if (arr[i] == position) {
			return i;
		}
	}
	std::cout << "count't found face_id . error " << std::endl;
}
static void delete_array_data_cos(float * arr, unsigned int size, unsigned int index) {
	for (unsigned int i = index; i < size - 1; i++)
	{
		arr[i] = arr[i + 1];
	}
}

static void delete_array_data_index(unsigned int * arr, unsigned int size, unsigned int index, unsigned int position) {
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

int Faster_Search::remove_one_person(std::string face_id) {
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
	std::string * prim_id1 = new std::string[test_size];
	std::string * prim_id2 = new std::string[test_size];
	std::string * prim_id3 = new std::string[test_size];
	std::cout << "all people length is:" << data_graph->face_info.datas.size() << std::endl;

	clock_t startTime, endTime;
	startTime = clock();//计时开始
	for (unsigned int i = 0; i < test_size; i++)
	{
		prim_id1[i] = get_near_id(data_graph, test_data[i], threhold1);
		
		/*if (i < 10)
			std::cout << " " << prim_id1[i];*/
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

	startTime = clock();//计时开始
	for (unsigned int i = 0; i < test_size; i++)
	{
		prim_id3[i] = search2(test_data[i]);
	}
	endTime = clock();//计时结束
	std::cout << "3 data run time is: " << (double)(endTime - startTime) << "ms" << std::endl;

	float accuracy = 0;
	for (unsigned int i = 0; i < test_size; i++)
	{
		if (prim_id1[i] == prim_id2[i]) {
			accuracy += 1;
		}
	}
	std::cout << "accuracy is:" << accuracy / test_size << std::endl;
	
	accuracy = 0;
	for (unsigned int i = 0; i < test_size; i++)
	{
		if (prim_id1[i] == prim_id3[i]) {
			accuracy += 1;
		}
	}
	std::cout << "2 accuracy is:" << accuracy / test_size << std::endl;
	delete[] prim_id1;
	delete[] prim_id2;
	delete[] prim_id3;
	return accuracy / test_size;
}


bool Faster_Search::storage_graph(std::string graph_file) {
	std::ofstream os(graph_file, std::ios::binary);
	if (!os) {
		return false;
	}
	uint32_t length = data_graph->face_info.datas.size();

	os.write(reinterpret_cast<char const *>(&data_graph->dim), sizeof(data_graph->dim));
	os.write(reinterpret_cast<char const *>(&length), sizeof(length));

	for (size_t i = 0; i < search_use_dim; i++)
	{
		os.write(reinterpret_cast<char const *>(&data_graph->base_coss[i][0]), length * sizeof(data_graph->base_coss[i][0]));
		os.write(reinterpret_cast<char const *>(&data_graph->base_indexs[i][0]), length * sizeof(data_graph->base_indexs[i][0]));
	}
	os.write(reinterpret_cast<char const *>(&data_graph->sort_var_index[0]), search_use_dim * sizeof(data_graph->sort_var_index[0]));
	//check with read.
	//for (size_t i = 0; i < 50; i++)
	//{
	//	std::cout << data_graph->base_coss[i][1] << " " << data_graph->base_indexs[i][1] << " " << data_graph->sort_var_index[i] << " ";
	//}
	//std::cout << std::endl;
	//for (size_t i = 0; i < 50; i++)
	//{
	//	std::cout << data_graph->base_coss[i][length-1] << " " << data_graph->base_indexs[i][length-1] << " ";
	//}
	//std::cout << std::endl;
	//os.flush();
	os.close();

	return true;
}

bool Faster_Search::restorage_graph(std::string graph_file, const unsigned int max_length) {
	std::ifstream is(graph_file.c_str(), std::ios::binary);
	if (!is) {
		return false;
	}
	//is.seekg(0, std::ios::beg);
	//is.seekg(sizeof(dim), std::ios::cur);
	//is.seekg(0, std::ios::end);
	uint32_t dim=0, length=0;
	is.read((char*)&dim, sizeof(dim));
	is.read((char*)&length, sizeof(length));
	//std::cout << "dim is:" << dim << " length is:" << length << std::endl;

	for (size_t i = 0; i < search_use_dim; i++)
	{
		float * cos = new float[max_length];
		unsigned int * index = new unsigned int[max_length];
		is.read((char*)&cos[0], length * sizeof(cos[0]));
		is.read((char*)&index[0], length * sizeof(index[0]));
		data_graph->base_coss.push_back(cos);
		data_graph->base_indexs.push_back(index);
	}
	data_graph->sort_var_index = new int[dim];
	is.read((char*)&data_graph->sort_var_index[0], search_use_dim * sizeof(data_graph->sort_var_index[0]));
	is.close();
	//check with write.
	//for (size_t i = 0; i < 50; i++)
	//{
	//	std::cout << data_graph->base_coss[i][1] << " " << data_graph->base_indexs[i][1]<<" "<< data_graph->sort_var_index[i]<<" ";
	//}
	//std::cout << std::endl;

	//for (size_t i = 0; i < 50; i++)
	//{
	//	std::cout << data_graph->base_coss[i][length - 1] << " " << data_graph->base_indexs[i][length - 1] << " ";
	//}
	//std::cout << std::endl;
	return true;
}

bool Faster_Search::save_face_data(std::string face_data, std::vector<float *>feature_data_test) {
	std::ofstream os(face_data, std::ios::binary);
	if (!os) {
		return false;
	}
	uint32_t length = data_graph->face_info.datas.size();
	os.write(reinterpret_cast<char const *>(&data_graph->dim), sizeof(data_graph->dim));
	os.write(reinterpret_cast<char const *>(&length), sizeof(length));
	for (size_t i = 0; i < length; i++)
	{
		os.write(reinterpret_cast<char const *>(&data_graph->face_info.datas[i][0]),
			data_graph->dim * sizeof(data_graph->face_info.datas[i][0]));

		os.write(reinterpret_cast<char const *>(&data_graph->face_info.face_ids[i]), sizeof(data_graph->face_info.face_ids[i]));
	}
	for (size_t i = 0; i < 1000; i++)
	{
		os.write(reinterpret_cast<char const *>(&feature_data_test[i][0]),
			data_graph->dim * sizeof(feature_data_test[i][0]));
	}
	//os.flush();
	os.close();
	//check data.
	for (size_t i = 0; i < 10; i++)
	{
		std::cout << data_graph->face_info.datas[i][0] << " " << data_graph->face_info.face_ids[i] << " ";
	}
	std::cout << std::endl;
	return true;
}
bool Faster_Search::load_face_data(std::string face_data, std::vector<float *> &feature_data_test) {
	std::ifstream is(face_data.c_str(), std::ios::binary);
	if (!is) {
		return false;
	}
	uint32_t dim = 0, length = 0;
	is.read((char*)&dim, sizeof(dim));
	is.read((char*)&length, sizeof(length));
	std::cout << "dim is:" << dim << " length is:" << length << std::endl;

	for (size_t i = 0; i < length; i++)
	{
		float * data = new float[dim];
		is.read((char*)&data[0], dim * sizeof(data[0]));
		
		std::string face_id;
		is.read((char*)&face_id, sizeof(face_id));
		
		data_graph->face_info.datas.push_back(data);
		data_graph->face_info.face_ids.push_back(face_id);
	}
	for (size_t i = 0; i < 1000; i++)
	{
		float * data = new float[dim];
		is.read((char*)&data[0], dim * sizeof(data[0]));
		feature_data_test.push_back(data);
	}

	is.close();
	//check data.
	//for (size_t i = 0; i < 10; i++)
	//{
	//	std::cout << data_graph->face_info.datas[i][0] << " " << data_graph->face_info.face_ids[i] << " ";
	//}
	//std::cout << std::endl;
	return true;
}

Faster_Search::~Faster_Search() {
	int dim = data_graph->base_coss.size();
	for (unsigned int i = 0; i < dim; i++)
	{
		delete [] data_graph->base_coss[i];
		delete[] data_graph->base_indexs[i];
	}
	int length = data_graph->face_info.datas.size();
	for (unsigned int i = 0; i < length; i++)
	{
		delete[] data_graph->face_info.datas[i];
	}
	delete[] data_graph->sort_var_index;
	delete[] base;
	
}