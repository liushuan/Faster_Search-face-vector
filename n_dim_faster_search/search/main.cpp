#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <map>

using namespace std;

struct Data_Matrix {
	unsigned int data_length;
	unsigned int dim;
	std::vector<std::vector<float>> data;
	std::vector<int> var_dims;
};

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

float getAngleDistance(vector<float> feature1, vector<float> feature2) {
	int length = feature1.size();
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
}
float getAngleDistance(float feature1[], float feature2[], int dim) {
	float distanceCos = 0;
	float nor_feature1 = 1, nor_feature2 = 1;
	for (int j = 0; j < dim; j++)
	{
		distanceCos += feature1[j] * feature2[j];
		//nor_feature1 += pow(feature1[j], 2);
		//nor_feature2 += pow(feature2[j], 2);
	}
	distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
	return distanceCos;
}
int get_near_list(Data_Matrix data_mat, std::vector<float> test, float threod) {
	float max_cos = -1;
	int index = -1;
	for (size_t i = 0; i < data_mat.data_length; i++)
	{
		float dist_cos = getAngleDistance(data_mat.data[i], test);
		if (dist_cos > max_cos && dist_cos > threod) {
			max_cos = dist_cos;
			index = i;
		}
		//std::cout << dist_cos << " ";
	}
	return index;
}

std::vector<std::multimap<float, int>> build_graph(Data_Matrix data_mat, float **& array_index) {

	std::vector<std::multimap<float, int>> all_index;
	for (size_t i = 0; i < data_mat.dim; i++)
	{
		std::vector<float> base(data_mat.dim, 0);
		base[i] = 1;
		std::multimap<float, int> cos_index;
		for (size_t j = 0; j < data_mat.data_length; j++)
		{
			float dis = getAngleDistance(data_mat.data[j], base);
			cos_index.insert(std::make_pair(dis, j));
		}

		all_index.push_back(cos_index);

	}
	array_index = new float*[2];
	array_index[0] = new float[data_mat.dim * data_mat.data_length];
	array_index[1] = new float[data_mat.dim * data_mat.data_length];


	for (size_t i = 0; i < data_mat.dim; i++)
	{
		int j = 0;
		for (auto element : all_index[i])
		{
			array_index[0][i*data_mat.data_length + j] = element.first;
			array_index[1][i*data_mat.data_length + j] = element.second;
			j++;
		}
	}

	return all_index;
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

int search(Data_Matrix data_mat, std::vector<std::multimap<float, int>> all_index, std::vector<float> test, float threod) {
	std::vector<int> results;

	clock_t startTime, endTime;
	startTime = clock();

	float th = 0.1f;
	std::vector<float> base(data_mat.dim, 0);

	std::multimap<float, int>::iterator iter_lower;
	std::multimap<float, int>::iterator iter_upper;
	std::vector<int> res;
	for (size_t i = 0; i < data_mat.dim; i++)
	{
		if (i == 0) {
			base[0] = 1;
		}
		else {
			base[i - 1] = 0;
			base[i] = 1;
		}
		float diss = getAngleDistance(base, test);

		if ((diss - th) < -1) {
			iter_lower = all_index[i].begin();
		}
		else {
			iter_lower = all_index[i].lower_bound(diss - th);
		}

		if (diss + th > 1) {
			iter_upper = all_index[i].end();
		}
		else {
			iter_upper = all_index[i].upper_bound(diss + th);
		}

		res.resize(0);
		for (std::multimap<float, int>::iterator iter = iter_lower; iter != iter_upper; iter++)
		{
			int index = iter->second;
			res.push_back(index);
		}
		if (i == 0) {
			results = res;
		}
		else {
			results = merge(results, res);
		}

		if (results.size() < 10) {
			std::cout << "results size < 10 , and break ." << std::endl;
			break;
		}

	}
	endTime = clock();

	std::cout << " results size is:" << results.size() << " time" << (endTime - startTime) << std::endl;

	startTime = clock();
	float max_cos = -1;
	int resu_dex = -1;
	for (int i = 0; i < results.size(); i++)
	{
		float dist = getAngleDistance(data_mat.data[results[i]], test);
		if (dist > max_cos && dist > threod) {
			max_cos = dist;
			resu_dex = results[i];
		}
	}
	endTime = clock();
	std::cout << " time is:" << endTime - startTime << std::endl;
	return resu_dex;
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
int my_upper_bound(float *array, int size, float key)
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

void merge(int arr1[], int arr2[], int &L, int K) {

	int index = 0;
	for (size_t i = 0; i < L; i++)
	{
		for (size_t j = 0; j < K; j++)
		{
			if (arr1[i] == arr2[j]) {
				arr1[index++] = arr1[j];
				break;
			}
		}
	}
	L = index;
	//std::cout << "length is:" << L << " " << " " << K << " " << index << std::endl;
}

int search(Data_Matrix data_mat, float **& array_index, std::vector<float> test, float threod) {

	clock_t startTime, endTime;
	startTime = clock();
	int * results = new int[data_mat.data_length];
	int *res = new int[data_mat.data_length];
	float th = 0.1f;
	std::vector<float> base(data_mat.dim);
	int start = 0, end = data_mat.data_length;
	int K = 0, L = 0;
	for (size_t i = 0; i < data_mat.dim; i++)
	{
		if (i == 0) {
			base[0] = 1;
		}
		else {
			base[i - 1] = 0;
			base[i] = 1;
		}
		float diss = getAngleDistance(base, test);

		if ((diss - th) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(array_index[0] + i * data_mat.data_length, data_mat.data_length, diss - th);
		}

		if (diss + th > 1) {
			end = data_mat.data_length;
		}
		else {
			end = my_upper_bound(array_index[0] + i * data_mat.data_length, data_mat.data_length, diss + th) + 1;
		}
		K = 0;
		for (int j = start; j < end; j++)
		{
			res[K] = array_index[1][i * data_mat.data_length + j];
			if (i == 0) {
				results[L] = res[K];
				L++;
			}
			K++;
		}

		if (i > 0) {
			int index = 0;
			for (size_t ii = 0; ii < L; ii++)
			{
				for (size_t jj = 0; jj < K; jj++)
				{
					if (results[ii] == res[jj]) {
						results[index++] = res[jj];
						break;
					}
				}
			}
			L = index;
		}

		if (L < 10) {
			//std::cout << "results size < 10 , and break ."<< std::endl;
			break;
		}

	}
	endTime = clock();

	std::cout << " results size is:" << L << " time" << (endTime - startTime) << std::endl;

	startTime = clock();
	float max_cos = -1;
	int resu_dex = -1;
	for (int i = 0; i < L; i++)
	{
		float dist = getAngleDistance(data_mat.data[results[i]], test);
		if (dist > max_cos && dist > threod) {
			max_cos = dist;
			resu_dex = results[i];
		}
	}
	endTime = clock();
	std::cout << " time is:" << endTime - startTime << std::endl;
	return resu_dex;
}

int search1(Data_Matrix data_mat, float **& array_index, std::vector<float> test, float threod) {

	clock_t startTime, endTime;
	int *results = new int[data_mat.data_length];
	memset(results, 0, data_mat.data_length*sizeof(int));

	float th = 0.1f;
	std::vector<float> base(data_mat.dim, 0);
	int start = 0, end = data_mat.data_length;
	for (size_t i = 0; i < data_mat.dim; i++)
	{
		if (i == 0) {
			base[0] = 1;
		}
		else {
			base[i - 1] = 0;
			base[i] = 1;
		}
		float diss = getAngleDistance(base, test);

		if ((diss - th) < -1) {
			start = 0;
		}
		else {
			start = my_lower_bound(array_index[0] + i * data_mat.data_length, data_mat.data_length, diss - th);
		}

		if (diss + th > 1) {
			end = data_mat.data_length;
		}
		else {
			end = my_upper_bound(array_index[0] + i * data_mat.data_length, data_mat.data_length, diss + th) + 1;
		}

		for (int j = start; j < end; j++)
		{
			results[(int)array_index[1][i * data_mat.data_length + j]]++;
		}
	}
	float max_cos = -1;
	int resu_dex = -1;
	for (int i = 0; i < data_mat.data_length; i++)
	{
		if (results[i] == data_mat.dim) {
			float dist = getAngleDistance(data_mat.data[i], test);
			if (dist > max_cos && dist > threod) {
				max_cos = dist;
				resu_dex = i;
			}
		}
	}
	return resu_dex;
}

int search2(Data_Array data_array, float **& array_index, int *& vars_index,  float test2[], float threod) {
	clock_t startTime, endTime;
	int *results = new int[data_array.data_length];
	memset(results, 0, data_array.data_length*sizeof(int));

	float th = 0.1f;
	float *base = new float[data_array.dim];
	memset(base, 0, data_array.dim*sizeof(float));

	int start = 0, end = data_array.data_length;
	int use_dim = 0;
	for (size_t i = 0; i < data_array.dim; i++)
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
		for (int j = start; j < end; j++)
		{
			results[(int)(array_index[1][i * data_array.data_length + j])]++;
		}

		base[i] = 0;
		use_dim++;
		if (use_dim == 50) {
			break;
		}
	}

	float max_cos = -1;
	int resu_dex = -1;
	int num_recall = 0;
	for (int i = 0; i < data_array.data_length; i++)
	{
		if (results[i] == use_dim) {
			float dist = getAngleDistance(data_array.data + i * data_array.dim, test2, data_array.dim);
			if (dist > max_cos && dist > threod) {
				max_cos = dist;
				resu_dex = i;
			}
			num_recall++;
		}
	}
	//std::cout << "----------search2 num_recall is:" << num_recall << std::endl;
	return resu_dex;
}

int search3(Data_Array data_array, float **& array_index, int *& vars_index, float test2[], float threod) {
	float th = 0.1f;

	int *results = new int[data_array.data_length];
	memset(results, 0, data_array.data_length * sizeof(int));
	float *base = new float[data_array.dim];
	memset(base, 0, data_array.dim * sizeof(float));

	int start = 0, end = data_array.data_length;
	int use_dim = 0;
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
		if (use_dim == 40) {
			break;
		}
	}

	float max_cos = -1;
	int resu_dex = -1;
	int num_recall = 0;
	for (int i = 0; i < data_array.data_length; i++)
	{
		if (results[i] == use_dim) {
			float dist = getAngleDistance(data_array.data + i * data_array.dim, test2, data_array.dim);
			if (dist > max_cos && dist > threod) {
				max_cos = dist;
				resu_dex = i;
			}
			num_recall++;
		}
	}
	//std::cout << "----------search3 num_recall is:" << num_recall << std::endl;
	return resu_dex;
}

int main()
{
	//std::default_random_engine e; //引擎
	//std::normal_distribution<double> n(4, 1.5); //均值, 方差
	srand(time(NULL));
	Data_Matrix data_mat;
	Data_Array data_array;
	data_mat.data_length = 10000;
	data_mat.dim = 512;

	data_array.data_length = data_mat.data_length;
	data_array.dim = data_mat.dim;
	data_array.data = new float[(long)data_array.data_length * data_array.dim];

	for (size_t j = 0; j < data_mat.data_length; j++)
	{
		std::vector<float> vals(data_mat.dim);
		for (std::size_t i = 0; i < data_mat.dim; ++i) {
			//unsigned v = std::lround(n(e)); //取整-最近的整数
			vals[i] = rand() * std::pow(-1, rand() % 2);
		}
		normalize(vals, data_mat.dim);

		for (size_t i = 0; i < data_array.dim; i++)
		{
			data_array.data[j*data_array.dim + i] = vals[i];
		}

		data_mat.data.push_back(vals);
	}

	std::vector<float> test(data_mat.dim);
	float *test2 = new float[data_array.dim];
	for (std::size_t i = 0; i < data_mat.dim; ++i) {
		//unsigned v = std::lround(n(e)); //取整-最近的整数
		test[i] = rand() * std::pow(-1, rand() % 2);
	}
	normalize(test, data_mat.dim);
	for (size_t i = 0; i < data_array.dim; i++)
	{
		test2[i] = test[i];
	}

	clock_t startTime, endTime;
	startTime = clock();//计时开始
	int index = get_near_list(data_mat, test, 0.4);
	endTime = clock();//计时结束
	cout << "The run time is: " << (double)(endTime - startTime) << "ms" << endl;
	std::cout << "run index is:" << index << " dis is:" << getAngleDistance(data_mat.data[index], test) << std::endl;

	startTime = clock();//计时开始
	float ** array_index;
	std::vector<std::multimap<float, int>> all_index = build_graph(data_mat, array_index);
	endTime = clock();//计时结束
	cout << "Build time is: " << (double)(endTime - startTime) << "ms" << endl;

	startTime = clock();//计时开始
	int resu_dex = search(data_mat, all_index, test, 0.4);
	endTime = clock();//计时结束
	cout << "1 search time is: " << (double)(endTime - startTime) << "ms" << endl;
	std::cout << "1 search index is:" << resu_dex << " dis is:" << getAngleDistance(data_mat.data[resu_dex], test) << std::endl;


	startTime = clock();//计时开始
	int resu_dex2 = search(data_mat, array_index, test, 0.4f);
	endTime = clock();//计时结束
	cout << "2 search time is: " << (double)(endTime - startTime) << "ms" << endl;
	std::cout << "2 search index is:" << resu_dex2 << " dis is:" << getAngleDistance(data_mat.data[resu_dex2], test) << std::endl;

	startTime = clock();//计时开始
	int resu_dex3 = search1(data_mat, array_index, test, 0.4f);
	endTime = clock();//计时结束
	cout << "3 search time is: " << (double)(endTime - startTime) << "ms" << endl;
	std::cout << "3 search index is:" << resu_dex3 << " dis is:" << getAngleDistance(data_mat.data[resu_dex3], test) << std::endl;

	startTime = clock();//计时开始
	float ** array_index2;
	int * vars_index;
	build_graph2(data_array, array_index2, vars_index);
	endTime = clock();//计时结束
	cout << "4 Build time is: " << (double)(endTime - startTime) << "ms" << endl;

	startTime = clock();//计时开始
	int resu_dex4 = search2(data_array, array_index2, vars_index, test2, 0.4f);
	endTime = clock();//计时结束
	cout << "4 search time is: " << (double)(endTime - startTime) << "ms" << endl;
	std::cout << " " << resu_dex4 << std::endl;
	std::cout << "4 search index is:" << resu_dex4 << " dis is:" << getAngleDistance(data_array.data + resu_dex4*data_array.dim, test2, data_array.dim) << std::endl;


	int resu_dex5 = 0;
	startTime = clock();//计时开始
	for (size_t i = 0; i < 10; i++)
	{
		resu_dex5 = search3(data_array, array_index2, vars_index, test2, 0.4f);
	}
	endTime = clock();//计时结束
	cout << "5 search time is: " << (double)(endTime - startTime) << "ms" << endl;
	std::cout << " " << resu_dex5 << std::endl;
	std::cout << "5 search index is:" << resu_dex5 << " dis is:" << getAngleDistance(data_array.data + resu_dex5*data_array.dim, test2, data_array.dim) << std::endl;



	system("pause");

	return 0;
}