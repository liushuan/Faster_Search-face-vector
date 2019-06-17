#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <map>
#include <array>
using namespace std;

struct Data_Array {
	unsigned int data_length;
	unsigned int dim;
	float * data;
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

//float getAngleDistance(vector<float> feature1, vector<float> feature2) {
//	int length = feature1.size();
//	float distanceCos = 0;
//	float nor_feature1 = 1, nor_feature2 = 1;
//	for (int j = 0; j < length; j++)
//	{
//		distanceCos += feature1[j] * feature2[j];
//		//nor_feature1 += pow(feature1[j], 2);
//		//nor_feature2 += pow(feature2[j], 2);
//	}
//	distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
//	return distanceCos;
//}
float getAngleDistance(float *feature1, float *feature2, const int dim) {
	float distanceCos = 0;
	float nor_feature1 = 1, nor_feature2 = 1;
	for (int j = 0; j < dim; j++)
	{
		distanceCos += (*feature1) * (*feature2);
		feature1++;
		feature2++;
		//distanceCos += feature1[j] * feature2[j];
		//nor_feature1 += pow(feature1[j], 2);
		//nor_feature2 += pow(feature2[j], 2);
	}
	distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
	return distanceCos;
}

int get_near_list(Data_Array data_array, float test[], float threod) {
	float max_cos = -1;
	int index = -1;
	for (size_t i = 0; i < data_array.data_length; i++)
	{
		float dist_cos = getAngleDistance(data_array.data + i*data_array.dim, test, data_array.dim);
		if (dist_cos > max_cos && dist_cos > threod) {
			max_cos = dist_cos;
			index = i;
		}
		//std::cout << dist_cos << " ";
	}
	return index;
}
float ave(float a[], int n)
{
	float sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i];
	}
	return sum / n;
}
float variance(float a[], int n)
{
	float sum = 0;
	float average = ave(a, n);
	for (int i = 0; i < n; i++) {
		sum += (a[i] - average)*(a[i] - average);
	}
	return sum;  //没有 用 sum/n;
}

void bubble_sort1(float a[], int b[], int n)
{
	int i, j;
	for (int i = 0; i < n; i++)
	{
		b[i] = i;
	}

	for (i = n - 1; i>0; i--)
	{
		// 将a[0...i]中最大的数据放在...
		for (j = 0; j<i; j++)
		{
			if (a[j] < a[j + 1]) {
				swap(a[j], a[j + 1]);
				swap(b[j], b[j + 1]);
			}
		}
	}

}
void build_graph2(Data_Array data_array, float **& array_index, int *& vars_index) {

	std::vector<std::multimap<float, int>> all_index;
	float *base = new float[data_array.dim];
	memset(base, 0.0f, data_array.dim * sizeof(float));

	for (size_t i = 0; i < data_array.dim; i++)
	{
		base[i] = 1;
		std::multimap<float, int> cos_index;
		for (size_t j = 0; j < data_array.data_length; j++)
		{
			float dis = getAngleDistance(data_array.data + j * data_array.dim, base, data_array.dim);
			cos_index.insert(std::make_pair(dis, j));
		}

		all_index.push_back(cos_index);

		base[i] = 0;
	}
	array_index = new float*[2];
	array_index[0] = new float[data_array.dim * data_array.data_length];
	array_index[1] = new float[data_array.dim * data_array.data_length];
	//array_index[2] = new float[data_mat.dim];
	vars_index = new int[data_array.dim];

	for (size_t i = 0; i < data_array.dim; i++)
	{
		int j = 0;
		for (auto element : all_index[i])
		{
			array_index[0][i*data_array.data_length + j] = element.first;
			array_index[1][i*data_array.data_length + j] = element.second;
			j++;
		}
	}

	//get var_dims;
	float *var = new float[data_array.dim];
	for (size_t i = 0; i < data_array.dim; i++)
	{
		var[i] = variance(array_index[0] + i * data_array.data_length, data_array.data_length);
	}
	bubble_sort1(var, vars_index, data_array.dim);
}

template<size_t _Size>
void build_graph3(Data_Array data_array, std::vector<std::array<float, _Size>>& array_coss,  std::vector<std::array<int, _Size>>& array_indexs) {

	std::vector<std::multimap<float, int>> all_index;
	float *base = new float[data_array.dim];
	memset(base, 0.0f, data_array.dim * sizeof(float));

	for (size_t i = 0; i < data_array.dim; i++)
	{
		base[i] = 1;
		std::multimap<float, int> cos_index;
		for (size_t j = 0; j < data_array.data_length; j++)
		{
			float dis = getAngleDistance(data_array.data + j * data_array.dim, base, data_array.dim);
			cos_index.insert(std::make_pair(dis, j));
		}

		all_index.push_back(cos_index);
		base[i] = 0;
	}

	for (size_t i = 0; i < data_array.dim; i++)
	{
		int j = 0;
		std::array<float, _Size> array_cos;
		std::array<int, _Size> array_index;
		for (auto element : all_index[i])
		{
			array_cos[j] = element.first;
			array_index[j] = element.second;
			//array_index[0][i*data_array.data_length + j] = element.first;
			//array_index[1][i*data_array.data_length + j] = element.second;
			j++;
		}
		array_coss.push_back(array_cos);
		array_indexs.push_back(array_index);
	}

}


std::vector<int> merge(std::vector<int> t1, std::vector<int> t2) {
	int length1 = t1.size();
	int length2 = t2.size();
	std::vector<int> merg;
	for (size_t i = 0; i < length1; i++)
	{
		for (size_t j = 0; j < length2; j++)
		{
			if (t1[i] == t2[j]) {
				merg.push_back(t1[i]);
				break;
			}
		}
	}
	return merg;
}
template <size_t _Size>
int my_lower_bound(std::array<float, _Size> array, int size, float key)
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
template <size_t _Size>
int my_upper_bound(std::array<float, _Size> array, int size, float key)
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


int my_lower_bound(float *arr, int size, float key)
{
	int first = 0, last = size - 1;
	int middle, pos = 0;       //需要用pos记录第一个大于等于key的元素位置

	while (first < last)
	{
		middle = (first + last) / 2;
		if (arr[middle] < key) {      //若中位数的值小于key的值，我们要在右边子序列中查找，这时候pos可能是右边子序列的第一个
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
int my_upper_bound(float *arr, int size, float key)
{
	int first = 0, last = size - 1;
	int middle, pos = 0;

	while (first < last)
	{
		middle = (first + last) / 2;
		if (arr[middle] > key) {     //当中位数大于key时，last不动，让first不断逼近last
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

int search2(Data_Array &data_array, float **& array_index, int *& vars_index, float test2[], float threod) {
	float th = 0.08f;
	int *results = new int[data_array.data_length];
	memset(results, 0, data_array.data_length * sizeof(int));

	
	float *base = new float[data_array.dim];
	memset(base, 0, data_array.dim * sizeof(float));

	int start = 0, end = data_array.data_length;
	int use_dim = 0;
	for (size_t i = 0; i < data_array.dim; ++i)
	{
		base[i] = 1;
		float diss = getAngleDistance(base, test2, data_array.dim);

		if ((diss - th) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(array_index[0] + i * data_array.data_length, data_array.data_length, diss - th);
		}

		if (diss + th > 1) {
			end = data_array.data_length;
		}
		else {
			end = my_upper_bound(array_index[0] + i * data_array.data_length, data_array.data_length, diss + th) + 1;
		}
		//std::cout << end - start << " " << vars_index[i] << std::endl;
		for (int j = start; j < end; ++j)
		{
			results[(int)(array_index[1][i * data_array.data_length + j])]++;
		}

		base[i] = 0;
		use_dim++;
		if (use_dim == 20) {
			break;
		}
	}

	float max_cos = -1;
	int resu_dex = -1;
	//int num_recall = 0;
	for (int i = 0; i < data_array.data_length; ++i)
	{
		if (results[i] == use_dim) {
			float dist = getAngleDistance(data_array.data + i * data_array.dim, test2, data_array.dim);
			if ((dist > max_cos) && (dist > threod)) {
				max_cos = dist;
				resu_dex = i;
			}
			//num_recall++;
		}
	}
	//std::cout << "----------search2 num_recall is:" << num_recall << std::endl;
	return resu_dex;
}

int search3(Data_Array data_array, float **& array_index, int *& vars_index, float test2[], float threod) {
	float th = 0.08f;
	float max_cos = -1.0f;
	int resu_dex = -1;
	int start = 0, end = data_array.data_length;
	int use_dim = 0;
	int *results = new int[data_array.data_length];
	memset(results, 0, data_array.data_length * sizeof(int));
	float *base = new float[data_array.dim];
	memset(base, 0, data_array.dim * sizeof(float));

	for (size_t i = 0; i < data_array.dim; i++)
	{
		base[vars_index[i]] = 1;
		float diss = getAngleDistance(base, test2, data_array.dim);

		if ((diss - th) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(array_index[0] + vars_index[i] * data_array.data_length, data_array.data_length, diss - th);
		}

		if (diss + th > 1) {
			end = data_array.data_length;
		}
		else {
			end = my_upper_bound(array_index[0] + vars_index[i] * data_array.data_length, data_array.data_length, diss + th) + 1;
		}
		//std::cout << end - start << " " << vars_index[i] << std::endl;
		for (int j = start; j < end; ++j)
		{
			results[(int)(array_index[1][vars_index[i] * data_array.data_length + j])]++;
		}
		base[vars_index[i]] = 0;
		use_dim++;
		if (use_dim == 20) {
			break;
		}
	}

	//int num_recall = 0;
	for (int i = 0; i < data_array.data_length; ++i)
	{
		if (results[i] == use_dim) {
			float dist = getAngleDistance(data_array.data + i * data_array.dim, test2, data_array.dim);
			if ((dist > max_cos) && (dist > threod)) {
				max_cos = dist;
				resu_dex = i;
			}
			//num_recall++;
		}
	}
	//std::cout << "----------search3 num_recall is:" << num_recall << std::endl;
	return resu_dex;
}

template <size_t _Size>
int search4(Data_Array &data_array,
	std::vector<std::array<float, _Size>> &array_coss,
	std::vector<std::array<int, _Size>> & array_indexs,
	float test2[], float threod) {
	float th = 0.08f;
	int *results = new int[data_array.data_length];
	memset(results, 0, data_array.data_length * sizeof(int));


	float *base = new float[data_array.dim];
	memset(base, 0, data_array.dim * sizeof(float));

	int start = 0, end = data_array.data_length;
	int use_dim = 0;
	for (size_t i = 0; i < data_array.dim; ++i)
	{
		base[i] = 1;
		float diss = getAngleDistance(base, test2, data_array.dim);

		if ((diss - th) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(array_coss[i], data_array.data_length, diss - th);
		}

		if (diss + th > 1) {
			end = data_array.data_length;
		}
		else {
			end = my_upper_bound(array_coss[i], data_array.data_length, diss + th) + 1;
		}
		//std::cout << end - start << " " << vars_index[i] << std::endl;
		for (int j = start; j < end; ++j)
		{
			results[array_indexs[i][j]]++;
		}

		base[i] = 0;
		use_dim++;
		if (use_dim == 20) {
			break;
		}
	}

	float max_cos = -1;
	int resu_dex = -1;
	//int num_recall = 0;
	for (int i = 0; i < data_array.data_length; ++i)
	{
		if (results[i] == use_dim) {
			float dist = getAngleDistance(data_array.data + i * data_array.dim, test2, data_array.dim);
			if ((dist > max_cos) && (dist > threod)) {
				max_cos = dist;
				resu_dex = i;
			}
			//num_recall++;
		}
	}
	//std::cout << "----------search2 num_recall is:" << num_recall << std::endl;
	return resu_dex;
}


int search5(Data_Array data_array, float **& array_index, int *& vars_index, float test2[], float threod) {
	float th = 0.08f;
	int *results = new int[data_array.data_length];
	memset(results, 0, data_array.data_length * sizeof(int));

	float *base = new float[data_array.dim];
	memset(base, 0, data_array.dim * sizeof(float));

	int start = 0, end = data_array.data_length;
	int use_dim = 0;
	for (size_t i = 0; i < data_array.dim; ++i)
	{
		base[i] = 1;
		float diss = getAngleDistance(base, test2, data_array.dim);

		if ((diss - th) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(array_index[0] + i * data_array.data_length, data_array.data_length, diss - th);
		}

		if (diss + th > 1) {
			end = data_array.data_length;
		}
		else {
			end = my_upper_bound(array_index[0] + i * data_array.data_length, data_array.data_length, diss + th) + 1;
		}

		// 13ms
		float * cursor = array_index[1] + (i * data_array.data_length) + start;
		for (int j = start; j < end; ++j)
		{
			results[int(*(cursor++))]++;
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
	for (int i = 0; i < data_array.data_length; ++i)
	{
		if (results[i] == use_dim) {
			float dist = getAngleDistance(data_array.data + i * data_array.dim, test2, data_array.dim);
			if ((dist > max_cos) && (dist > threod)) {
				max_cos = dist;
				resu_dex = i;
			}
		//	//num_recall++;
		}
	}
	//std::cout << "----------search2 num_recall is:" << num_recall << std::endl;
	return resu_dex;
}


Data_Array get_data(const int number, const int dim) {
	Data_Array data_array;
	data_array.data_length = number;
	data_array.dim = dim;
	data_array.data = new float[(long)data_array.data_length * data_array.dim];
	for (size_t j = 0; j < data_array.data_length; j++)
	{
		std::vector<float> vals(data_array.dim);
		for (std::size_t i = 0; i < data_array.dim; ++i) {
			//unsigned v = std::lround(n(e)); //取整-最近的整数
			vals[i] = rand() * std::pow(-1, rand() % 2);
		}
		normalize(vals, data_array.dim);

		for (size_t i = 0; i < data_array.dim; i++)
		{
			data_array.data[j*data_array.dim + i] = vals[i];
		}
	}
	return data_array;
}

int main()
{
	srand(time(NULL));
	float threhold = 0.1f;
	const int data_num = 10000;
	const int data_dim = 512;
	Data_Array train_data = get_data(data_num, data_dim);
	const int test_size = 100;
	Data_Array test_data = get_data(test_size, data_dim);
	int * prim_id1 = new int[test_size];
	int * prim_id2 = new int[test_size];

	clock_t startTime, endTime;
	startTime = clock();//计时开始
	for (size_t i = 0; i < test_size; i++)
	{
		prim_id1[i] = get_near_list(train_data, test_data.data + i * test_data.dim, threhold);
	}
	endTime = clock();//计时结束
	cout << "all data run time is: " << (double)(endTime - startTime) << "ms" << endl;
	
	startTime = clock();//计时开始
	float ** array_index;
	int * vars_index;
	build_graph2(train_data, array_index, vars_index);
	endTime = clock();//计时结束
	cout << "Build time is: " << (double)(endTime - startTime) << "ms" << endl;

	std::vector<std::array<float, data_num>> array_coss;
	std::vector<std::array<int, data_num>> array_indexs;
	build_graph3<data_num>(train_data, array_coss, array_indexs);

	int *results = new int[train_data.data_length];
	memset(results, 0, train_data.data_length * sizeof(int));
	float *base = new float[train_data.dim];
	memset(base, 0, train_data.dim * sizeof(float));

	startTime = clock();//计时开始
	for (size_t j = 0; j < test_size; j++)
	{
		prim_id2[j] = search2(train_data, array_index, vars_index, test_data.data + j * test_data.dim, threhold);
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

	startTime = clock();//计时开始
	for (size_t j = 0; j < test_size; j++)
	{
		prim_id2[j] = search5(train_data, array_index, vars_index, test_data.data + j * test_data.dim, threhold);
	}
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