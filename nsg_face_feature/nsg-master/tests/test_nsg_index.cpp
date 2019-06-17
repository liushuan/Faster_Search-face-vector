//
// Created by 付聪 on 2017/6/21.
//

#include <efanna2e/index_nsg.h>
#include <efanna2e/util.h>

#include <fstream>

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

void load_face_feature_data(std::string filename, int dim, float *&data, unsigned &num) {
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

		normalize(data+ i*dim, dim);

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
	data = new float[(size_t)num * (size_t)dim];

	in.seekg(0, std::ios::beg);
	for (size_t i = 0; i < num; i++) {
		in.seekg(4, std::ios::cur);
		in.read((char*)(data + i * dim), dim * 4);
	}
	in.close();
}

int main(int argc, char** argv) {
	srand((unsigned)time(NULL));
	if (argc != 7) {
		std::cout << argv[0] << " data_file nn_graph_path L R C save_graph_file"
			<< std::endl;
		exit(-1);
	}
	float* data_load = NULL;
	unsigned points_num, dim = 512;
	//load_data(argv[1], data_load, points_num, dim);
	//std::string filenames = "D:\\work_space\\work_code\\nsg\\source_code\\kgraph\\mnist\\train-images.idx3-ubyte";
	//load_mnist_data(filenames, data_load, points_num, dim);
	
	std::string filenames = "E:\\work_space\\source_code\\face_search\\nsg_windows\\face.feature";
	load_face_feature_data(filenames, dim, data_load, points_num);

	std::cout << "dim:" << dim << " num:" << points_num << std::endl;
	std::string nn_graph_path(argv[2]);
	std::cout << "file is:" << argv[2] << std::endl;
	unsigned L = (unsigned)atoi(argv[3]);
	unsigned R = (unsigned)atoi(argv[4]);
	unsigned C = (unsigned)atoi(argv[5]);
	std::cout << "L:" << L << " R:" << R << " C:" << C << std::endl;

	// data_load = efanna2e::data_align(data_load, points_num, dim);//one must
	// align the data before build
	efanna2e::IndexNSG index(dim, points_num, efanna2e::L2, nullptr);

	auto s = std::chrono::high_resolution_clock::now();
	efanna2e::Parameters paras;
	paras.Set<unsigned>("L", L);
	paras.Set<unsigned>("R", R);
	paras.Set<unsigned>("C", C);
	paras.Set<std::string>("nn_graph_path", nn_graph_path);

	index.Build(points_num, data_load, paras);
	auto e = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = e - s;

	std::cout << "indexing time: " << diff.count() << "\n";

	index.Save(argv[6]);

	/*std::string label1 = "D:\\work_space\\work_code\\nsg\\source_code\\kgraph\\mnist\\train-labels.idx1-ubyte";
	std::vector<double>labels;
	read_Mnist_Label(label1, labels);
	index.Show(labels);*/
	//std::vector<unsigned> tmp(5);
	//index.Search(data_load, data_load, 5, paras, tmp.data());
	//std::cout <<"label is:"<< tmp[0] << std::endl;


	system("pause");
	return 0;
}
