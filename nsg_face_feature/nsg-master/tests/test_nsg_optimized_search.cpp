//
// Created by 付聪 on 2017/6/21.
//

#include <efanna2e/index_nsg.h>
#include <efanna2e/util.h>
#include <chrono>
#include <string>
#include<opencv2\opencv.hpp>

int ReverseInt(int i)
{
	unsigned char ch1, ch2, ch3, ch4;
	ch1 = i & 255;
	ch2 = (i >> 8) & 255;
	ch3 = (i >> 16) & 255;
	ch4 = (i >> 24) & 255;
	return((int)ch1 << 24) + ((int)ch2 << 16) + ((int)ch3 << 8) + ch4;
}
using std::cout;
using std::endl;
void load_mnist_data(std::string filename, float*& data, unsigned& num, unsigned& dim) {
	std::ifstream file(filename, std::ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		int n_rows = 0;
		int n_cols = 0;
		unsigned char label;
		file.read((char*)&magic_number, sizeof(magic_number));
		file.read((char*)&number_of_images, sizeof(number_of_images));
		file.read((char*)&n_rows, sizeof(n_rows));
		file.read((char*)&n_cols, sizeof(n_cols));
		magic_number = ReverseInt(magic_number);
		number_of_images = ReverseInt(number_of_images);
		n_rows = ReverseInt(n_rows);
		n_cols = ReverseInt(n_cols);

		cout << "magic number = " << magic_number << endl;
		cout << "number of images = " << number_of_images << endl;
		cout << "rows = " << n_rows << endl;
		cout << "cols = " << n_cols << endl;

		/////////////////////////////////////////
		dim = n_rows*n_cols;
		num = number_of_images;
		data = new float[(size_t)num * (size_t)dim];

		for (int i = 0; i < number_of_images; i++)
		{
			//cv::Mat dst(28, 28, cv::CV_8UC1);
			//cv::Mat dst(28, 28, CV_8UC1);
			for (int r = 0; r < n_rows; r++)
			{
				for (int c = 0; c < n_cols; c++)
				{
					unsigned char image = 0;
					file.read((char*)&image, sizeof(image));
					//tp.push_back(image);
					data[dim * i + r*n_cols + c] = static_cast<float>(image);
					//dst.at<uchar>(r, c) = image;
				}
			}
			//cv::imshow("img", dst);
			//cv::waitKey(0);
		}
	}
}

void read_Mnist_Label(std::string filename, std::vector<double>&labels)
{
	std::ifstream file(filename, std::ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		file.read((char*)&magic_number, sizeof(magic_number));
		file.read((char*)&number_of_images, sizeof(number_of_images));
		magic_number = ReverseInt(magic_number);
		number_of_images = ReverseInt(number_of_images);
		cout << "magic number = " << magic_number << endl;
		cout << "number of images = " << number_of_images << endl;
		for (int i = 0; i < number_of_images; i++)
		{
			unsigned char label = 0;
			file.read((char*)&label, sizeof(label));
			labels.push_back((double)label);
		}
	}
}


using std::string;
template <class Type>
Type stringToNum(const std::string& str)
{
	std::istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}
static void  splitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
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
static std::vector<std::string> readStringFromFileData(std::string filePath)
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
			cout << "buf is empty." << endl;
			continue;
		}
		data.push_back(buf);
	}
	fileA.close();
	return data;
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

void load_face_feature_data_train(std::string filename, int dim, float *&data, unsigned &num) {
	std::vector<std::string> data_all = readStringFromFileData(filename);
	int train_size = data_all.size() - 1000;
	num = train_size;
	data = new float[(size_t)num * (size_t)dim];

	for (size_t i = 0; i < train_size; i++)
	{
		std::vector<string> str_d;
		splitString(data_all[i], str_d, ",");
		/////////////////////////////////////////
		for (size_t j = 0; j < dim; j++)
		{
			data[i*dim + j] = stringToNum<float>(str_d[j]);
		}
		normalize(data + i*dim, dim);
	}
}

void load_face_feature_data_test(std::string filename, int dim, float *&data, unsigned &num) {
	std::vector<std::string> data_all = readStringFromFileData(filename);
	int test_start = data_all.size() - 1000;
	num = 1000;
	data = new float[(size_t)num * (size_t)dim];

	for (size_t i = test_start; i < data_all.size(); i++)
	{
		std::vector<string> str_d;
		splitString(data_all[i], str_d, ",");
		/////////////////////////////////////////
		for (size_t j = 0; j < dim; j++)
		{
			data[(i - test_start)*dim + j] = stringToNum<float>(str_d[j]);
		}
		normalize(data + (i - test_start)*dim, dim);
	}
}



void load_data(char* filename, float*& data, unsigned& num,
	unsigned& dim) {  // load data with sift10K pattern
	std::ifstream in(filename, std::ios::binary);
	if (!in.is_open()) {
		std::cout << "open file error" << std::endl;
		exit(-1);
	}
	in.read((char*)&dim, 4);

	in.seekg(0, std::ios::end);
	std::ios::pos_type ss = in.tellg();
	size_t fsize = (size_t)ss;
	num = (unsigned)(fsize / (dim + 1) / 4);

	std::cout << "data dimension: " << dim << " num:" << num << std::endl;

	data = new float[(size_t)num * (size_t)dim];

	in.seekg(0, std::ios::beg);
	for (size_t i = 0; i < num; i++) {
		in.seekg(4, std::ios::cur);
		in.read((char*)(data + i * dim), dim * 4);
	}
	in.close();
}

void save_result(const char* filename,
	std::vector<std::vector<unsigned> >& results) {
	std::ofstream out(filename, std::ios::binary | std::ios::out);

	for (unsigned i = 0; i < results.size(); i++) {
		unsigned GK = (unsigned)results[i].size();
		out.write((char*)&GK, sizeof(unsigned));
		out.write((char*)results[i].data(), GK * sizeof(unsigned));
	}
	out.close();
}

void save_result_text(const char* filename,
	std::vector<std::vector<unsigned> >& results, std::vector<double>label_train, std::vector<double>label_test) {
	std::ofstream out;
	out.open(filename, std::ios::trunc);

	for (unsigned i = 0; i < results.size(); i++) {
		unsigned GK = (unsigned)results[i].size();

		out << "Kernel:" << label_test[i] << std::endl;

		if (i < 5) {
			std::cout << "test :" << label_test[i] << std::endl;
		}

		std::string line = "";
		int length = results[i].size();
		for (size_t j = 0; j < length; j++)
		{
			line.append(std::to_string(label_train[results[i][j]]) + " ");

			if (i < 5) {
				std::cout << " " << label_train[results[i][j]];
			}
		}
		out << " Nei:" << line << std::endl;
		if (i < 5) {
			std::cout << std::endl;
		}
	}
	out.close();
}

void Show(std::vector<std::vector<unsigned>> final_graph_, std::vector<double>labels) {

	unsigned nd_ = final_graph_.size();
	for (unsigned i = 0; i < nd_; i++) {
		if (i < 5)
			std::cout << "Kernel:" << labels[i] << std::endl;
		unsigned GK = (unsigned)final_graph_[i].size();

		if (i < 5)
			for (size_t j = 0; j < GK; j++)
			{
				std::cout << " " << labels[final_graph_[i][j]];
			}
		if (i < 5) {
			std::cout << std::endl;
		}
	}
}

float compare_L2(const float a[],  const float  b[], unsigned size){
	float sum = 0;
	for (size_t i = 0; i < size; i++)
	{
		sum += std::pow((a[i] - b[i]), 2);
	}
	return std::sqrt(sum);
}
int get_near_list(const float * data_load, const float* test, int num_set, int dim) {
	float max_cos = 10000;
	int index = -1;
	for (size_t i = 0; i < num_set; i++)
	{
		float dist_cos = compare_L2(data_load + i * dim, test, dim);
		if (dist_cos < max_cos) {
			max_cos = dist_cos;
			index = i;
		}
	}
	return index;
}

int main(int argc, char** argv) {
	if (argc != 7) {
		std::cout << argv[0]
			<< " data_file query_file nsg_path search_L search_K result_path"
			<< std::endl;
		exit(-1);
	}

	//std::string filenames1 = "D:\\work_space\\work_code\\nsg\\source_code\\kgraph\\mnist\\train-images.idx3-ubyte";
	//std::string label1 = "D:\\work_space\\work_code\\nsg\\source_code\\kgraph\\mnist\\train-labels.idx1-ubyte";
	//std::string filenames2 = "D:\\work_space\\work_code\\nsg\\source_code\\kgraph\\mnist\\t10k-images.idx3-ubyte";
	//std::string label2 = "D:\\work_space\\work_code\\nsg\\source_code\\kgraph\\mnist\\t10k-labels.idx1-ubyte";
	std::string filenames = "E:\\work_space\\source_code\\face_search\\nsg_windows\\face.feature";
	
	float* data_load = NULL;
	unsigned points_num, dim=512;
	//load_data(argv[1], data_load, points_num, dim);
	//load_mnist_data(filenames1, data_load, points_num, dim);
	load_face_feature_data_train(filenames, dim, data_load, points_num);
	//std::vector<double> labels_train;
	//read_Mnist_Label(label1, labels_train);

	float* query_load = NULL;
	unsigned query_num, query_dim=512;
	//load_data(argv[2], query_load, query_num, query_dim);
	load_face_feature_data_test(filenames, query_dim, query_load, query_num);
	//std::vector<double> labels_test;
	//read_Mnist_Label(label2, labels_test);
	std::cout << "load finish :" << points_num << " : " << query_num << "  data: "<<data_load[0]<< " "<< query_load[0] << std::endl;
	assert(dim == query_dim);

	unsigned L = (unsigned)atoi(argv[4]);
	unsigned K = (unsigned)atoi(argv[5]);

	if (L < K) {
		std::cout << "search_L cannot be smaller than search_K!" << std::endl;
		exit(-1);
	}
	// data_load = efanna2e::data_align(data_load, points_num, dim);//one must
	// align the data before build query_load = efanna2e::data_align(query_load,
	// query_num, query_dim);
	efanna2e::IndexNSG index(dim, points_num, efanna2e::FAST_L2, nullptr);
	index.Load(argv[3]);
	index.OptimizeGraph(data_load);

	efanna2e::Parameters paras;
	paras.Set<unsigned>("L_search", L);
	paras.Set<unsigned>("P_search", L);

	std::vector<std::vector<unsigned> > res(query_num);
	for (unsigned i = 0; i < query_num; i++) res[i].resize(K);

	std::vector<int> res_id;


	load_face_feature_data_train(filenames, dim, data_load, points_num);

	auto s = std::chrono::high_resolution_clock::now();
	for (unsigned i = 0; i < query_num; i++) {
		index.SearchWithOptGraph(query_load + i * dim, K, paras, res[i].data());

		
	}
	auto e = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = e - s;
	std::cout << "search time: " << diff.count() << "\n";

	auto s1 = std::chrono::high_resolution_clock::now();
	for (unsigned i = 0; i < query_num; i++) {
		res_id.push_back(get_near_list(data_load, query_load + (i*dim), points_num, dim));
	}
	auto e1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff1 = e1 - s1;
	std::cout << "search2 time: " << diff1.count() << "\n";
	

	//save_result(argv[6], res);
	//print search result .
	//for (size_t i = 0; i < 10; i++)
	//{
	//	int length = res[i].size();
	//	std::cout << "near  is:" << res_id[i] << std::endl;
	//	for (size_t j = 0; j < length; j++)
	//	{
	//		std::cout << " " << res[i][j];
	//	}
	//	std::cout << std::endl;
	//}

	float acc = 0;
	for (size_t i = 0; i < query_num; i++)
	{
		int length = res[i].size();
		for (size_t j = 0; j < length; j++)
		{
			//说明最近的 K 个 包含 最near的元素
			if (res_id[i] == res[i][j]) {
				acc++;
				break;
			}
		}
	}
	std::cout << "acc rate is:" << acc / query_num << std::endl;


	//save_result_text(argv[6], res, labels_train, labels_test);

	system("pause");
	return 0;
}
