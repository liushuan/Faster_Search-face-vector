#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <assert.h>
#include <time.h>
#include <map>

using namespace std;

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

struct Face_Info {
	std::vector <float *> datas;
	std::vector <unsigned int> data_ids;
};
struct Data_FasterSearch {
	const unsigned int dim;

	Face_Info face_info;
	//std::vector < std::array<float, MaxLength>> base_coss;
	//std::vector <std::forward_list<float>> base_coss;
	//std::vector <std::forward_list<unsigned int>> base_indexs;
	//std::vector < std::array<unsigned int, MaxLength>> base_indexs;
	std::vector<unsigned int * >base_indexs;
	std::vector<float * >base_coss;
	Data_FasterSearch(const unsigned int dim) :dim(dim) {}
};

float getAngleDistance(float feature1[], float feature2[], int size) {
	float distanceCos = 0;
	float nor_feature1 = 1, nor_feature2 = 1;
	for (int j = 0; j < size; j++)
	{
		distanceCos += feature1[j] * feature2[j];
	}
	distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
	return distanceCos;
}

int get_near_list(Face_Info face_infos, float test[], float threod, int dim) {
	float max_cos = -1;
	unsigned int index = -1;
	int length = face_infos.datas.size();
	for (size_t i = 0; i < length; ++i)
	{
		float dist_cos = getAngleDistance(face_infos.datas[i], test, dim);
		if (dist_cos > max_cos && dist_cos > threod) {
			max_cos = dist_cos;
			index = face_infos.data_ids[i];
		}
		//std::cout << dist_cos << " ";
	}
	return index;
}

void build_graph(Data_FasterSearch & data_faster_search, unsigned int max_length) {

	std::vector<std::multimap<float, int>> all_index;
	float *base = new float[data_faster_search.dim];
	memset(base, 0.0f, data_faster_search.dim * sizeof(float));
	int length = data_faster_search.face_info.datas.size();
	for (size_t i = 0; i < data_faster_search.dim; i++)
	{
		base[i] = 1;
		std::multimap<float, int> cos_index;
		for (size_t j = 0; j < length; j++)
		{
			float dis = getAngleDistance(data_faster_search.face_info.datas[j], base, data_faster_search.dim);
			cos_index.insert(std::make_pair(dis, j));
		}

		all_index.push_back(cos_index);
		base[i] = 0;
	}
	for (size_t i = 0; i < data_faster_search.dim; i++)
	{
		float * cos = new float[max_length];
		unsigned int * index = new unsigned int[max_length];
		int j = 0;
		for (auto element : all_index[i])
		{
			cos[j] = element.first;
			index[j] = element.second;
			//array_cos[j] = element.first;
			//array_index[j] = element.second;
			j++;
		}
		data_faster_search.base_coss.push_back(cos);
		data_faster_search.base_indexs.push_back(index);
	}
}

int my_lower_bound(float *array, int size, float key)
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

int my_upper_bound(float * array, int size, float key)
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


void get_data(Face_Info & face_info, const int number, const int dim) {

	for (size_t j = 0; j < number; j++)
	{
		float * data = new float[dim];

		std::vector<float> vals(dim);
		for (std::size_t i = 0; i < dim; ++i) {
			//unsigned v = std::lround(n(e)); //取整-最近的整数
			vals[i] = rand() * std::pow(-1, rand() % 2);
		}
		normalize(vals, dim);
		for (size_t i = 0; i < dim; i++)
		{
			data[i] = vals[i];
		}
		face_info.datas.push_back(data);
		face_info.data_ids.push_back(j);
	}
}

unsigned int found_array_data(unsigned int *arr, unsigned int size, unsigned int position) {
	for (unsigned int i = 0; i < size; ++i)
	{
		if (arr[i] == position) {
			return i;
		}
	}
	std::cout << "count't found face_id . error " << std::endl;
	assert(false);
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

void remove_data(Data_FasterSearch & data_faster_search, unsigned int face_id) {
	unsigned int length = data_faster_search.face_info.datas.size();
	//find data and erase data.
	auto iter = std::find(data_faster_search.face_info.data_ids.begin(), data_faster_search.face_info.data_ids.end(), face_id);
	if (iter != data_faster_search.face_info.data_ids.end()) {
		unsigned int nPosition = distance(data_faster_search.face_info.data_ids.begin(), iter);
		data_faster_search.face_info.data_ids.erase(iter);

		auto begin = data_faster_search.face_info.datas.begin();
		//std::cout << "-----delete nPosition:" << nPosition << std::endl;
		data_faster_search.face_info.datas.erase(begin + nPosition);

		for (size_t i = 0; i < data_faster_search.dim; i++)
		{
			unsigned int index = found_array_data(data_faster_search.base_indexs[i], length, nPosition);
			delete_array_data_cos(data_faster_search.base_coss[i], length, index);
			delete_array_data_index(data_faster_search.base_indexs[i], length, index, nPosition);
		}
	}
	else {
		std::cout << "can't found face_id:" << face_id << std::endl;
	}
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
	for (size_t i = size; i > pos; i--)
	{
		arr[i] = arr[i - 1];
	}
	arr[pos] = data;
}

void add_data(Data_FasterSearch & data_faster_search, float * add_data, unsigned int data_id) {
	unsigned int length = data_faster_search.face_info.datas.size();
	//add data.
	data_faster_search.face_info.datas.push_back(add_data);
	data_faster_search.face_info.data_ids.push_back(data_id);

	float *base = new float[data_faster_search.dim];
	memset(base, 0.0f, data_faster_search.dim * sizeof(float));
	for (size_t i = 0; i < data_faster_search.dim; i++)
	{
		base[i] = 1;
		//std::multimap<float, int> cos_index;

		float dis = getAngleDistance(add_data, base, data_faster_search.dim);
		//cos_index.insert(std::make_pair(dis, data_faster_search.face_info.data_ids[j]));
		unsigned int pos = get_insert_pos(data_faster_search.base_coss[i], length, dis);
		//int pos = insert_array_data<float>(data_faster_search.base_coss[i], length, dis);
		insert_array<float>(data_faster_search.base_coss[i], length, dis, pos);
		insert_array<unsigned int>(data_faster_search.base_indexs[i], length, length, pos);
		base[i] = 0;
	}


}

int search(Data_FasterSearch & data_faster_search, float *test, float threod) {
	float th = 0.08f;
	unsigned int length = data_faster_search.face_info.datas.size();
	int *results = new int[length];
	memset(results, 0, length * sizeof(int));

	float *base = new float[data_faster_search.dim];
	memset(base, 0, data_faster_search.dim * sizeof(float));

	unsigned int start = 0, end = length;
	int use_dim = 0;
	for (size_t i = 0; i < data_faster_search.dim; ++i)
	{
		base[i] = 1;
		float diss = getAngleDistance(base, test, data_faster_search.dim);
		if ((diss - th) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(data_faster_search.base_coss[i], length, diss - th);
		}

		if (diss + th > 1) {
			end = length;
		}
		else {
			end = my_upper_bound(data_faster_search.base_coss[i], length, diss + th) + 1;
		}
		// 13ms
		unsigned int * cursor = data_faster_search.base_indexs[i] + start;
		for (unsigned int j = start; j < end; ++j)
		{
			results[*(cursor++)]++;
		}
		base[i] = 0;
		use_dim++;
		if (use_dim == 20) {
			break;
		}
	}

	float max_cos = -1;
	int resu_dex = -1;
	int num_recall = 0;
	//7ms
	for (int i = 0; i < length; ++i)
	{
		//if (i % 20 == 0)
		//std::cout << " " << results[i];
		if (results[i] == use_dim) {
			float dist = getAngleDistance(data_faster_search.face_info.datas[i], test, data_faster_search.dim);
			if ((dist > max_cos) && (dist > threod)) {
				max_cos = dist;
				resu_dex = data_faster_search.face_info.data_ids[i];
			}
			num_recall++;
		}
	}
	std::cout << "  " << num_recall;

	return resu_dex;
}

void check_data(Data_FasterSearch & data_faster_search) {
	std::cout << "face data info size:" << data_faster_search.face_info.datas.size() << " " << data_faster_search.face_info.data_ids.size() << std::endl;
	int length = data_faster_search.face_info.datas.size();

	for (size_t i = 0; i < length; i++)
	{
		std::cout << "face_id:" << data_faster_search.face_info.data_ids[i] << " " << data_faster_search.face_info.datas[i][0] << std::endl;
	}

	for (size_t i = 0; i < data_faster_search.dim; i++)
	{
		for (size_t j = 0; j < length; j++) {
			std::cout << " " << j << ":" << data_faster_search.base_coss[i][j];
			std::cout << " " << data_faster_search.base_indexs[i][j];
		}
		std::cout << std::endl;
	}

}

int main()
{
	srand(time(NULL));
	float threhold = 0.1f;
	const int MaxLength = 100000;
	const int data_num = 10000;
	const int data_dim = 512;

	//typedef float data_type;
	Data_FasterSearch data_faster_search(data_dim);
	get_data(data_faster_search.face_info, data_num, data_dim);

	const int test_size = 100;
	Face_Info test_data;
	get_data(test_data, test_size, data_dim);
	int * prim_id1 = new int[test_size];
	int * prim_id2 = new int[test_size];

	clock_t startTime, endTime;
	startTime = clock();//计时开始
	for (size_t i = 0; i < test_size; i++)
	{
		prim_id1[i] = get_near_list(data_faster_search.face_info, test_data.datas[i], threhold, data_dim);
	}
	endTime = clock();//计时结束
	cout << "all data run time is: " << (double)(endTime - startTime) << "ms" << endl;

	startTime = clock();//计时开始
	float ** array_index;
	//int * vars_index;
	build_graph(data_faster_search, MaxLength);
	endTime = clock();//计时结束
	cout << "Build time is: " << (double)(endTime - startTime) << "ms" << endl;

	startTime = clock();//计时开始
	for (size_t j = 0; j < test_size; j++)
	{
		prim_id2[j] = search(data_faster_search, test_data.datas[j], threhold);
	}
	endTime = clock();//计时结束
	cout << "all search time is: " << (double)(endTime - startTime) << "ms" << endl;

	float accuracy = 0;
	for (size_t i = 0; i < test_size; i++)
	{
		if (prim_id1[i] == prim_id2[i]) {
			accuracy += 1;
		}
	}
	std::cout << "add accuracy is:" << accuracy / test_size << std::endl;

	//check_data(data_faster_search);

	//startTime = clock();//计时开始
	for (size_t i = 0; i < 500; i++)
	{
		remove_data(data_faster_search, i);
	}
	endTime = clock();//计时结束
	cout << "remove one time is: " << (double)(endTime - startTime) / 500 << "ms" << endl;

	//check_data(data_faster_search);

	startTime = clock();//计时开始
	Face_Info data_;
	get_data(data_, 1000, data_dim);
	for (size_t i = 0; i < 1000; i++)
	{
		add_data(data_faster_search, data_.datas[i], data_.data_ids[i]);
	}
	endTime = clock();//计时结束
	cout << "add one time is: " << (double)(endTime - startTime) / 1000 << "ms" << endl;


	startTime = clock();//计时开始
	for (size_t i = 0; i < test_size; i++)
	{
		prim_id1[i] = get_near_list(data_faster_search.face_info, test_data.datas[i], threhold, data_dim);
		//std::cout << " " << prim_id1[i];
	}
	//std::cout << std::endl;
	endTime = clock();//计时结束
	cout << "2 all data run time is: " << (double)(endTime - startTime) << "ms" << endl;

	startTime = clock();
	for (size_t j = 0; j < test_size; j++)
	{
		prim_id2[j] = search(data_faster_search, test_data.datas[j], threhold);
		//std::cout << " " << prim_id2[j];
	}
	//std::cout << std::endl;
	endTime = clock();//计时结束
	cout << "2 all search time is: " << (double)(endTime - startTime) << "ms" << endl;

	accuracy = 0;
	for (size_t i = 0; i < test_size; i++)
	{
		if (prim_id1[i] == prim_id2[i]) {
			accuracy += 1;
		}
	}
	std::cout << "2 add accuracy is:" << accuracy / test_size << std::endl;

	system("pause");

	return 0;
}