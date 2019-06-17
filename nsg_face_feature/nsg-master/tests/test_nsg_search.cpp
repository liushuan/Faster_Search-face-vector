//
// Created by 付聪 on 2017/6/21.
//

#include <efanna2e/index_nsg.h>
#include <efanna2e/util.h>

using namespace std;

void load_data(char* filename, float*& data, unsigned& num,
	unsigned& dim) {  // load data with sift10K pattern
	std::ifstream in(filename, std::ios::binary);
	if (!in.is_open()) {
		std::cout << "open file error" << std::endl;
		exit(-1);
	}
	in.read((char*)&dim, 4);
	std::cout << "data dimension: " << dim << std::endl;
	in.seekg(0, std::ios::end);
	std::ios::pos_type ss = in.tellg();
	size_t fsize = (size_t)ss;
	num = (unsigned)(fsize / (dim + 1) / 4);
	data = new float[num * dim * sizeof(float)];

	in.seekg(0, std::ios::beg);
	for (size_t i = 0; i < num; i++) {
		in.seekg(4, std::ios::cur);
		in.read((char*)(data + i * dim), dim * 4);
	}
	in.close();
}

void save_result(char* filename, std::vector<std::vector<unsigned> >& results) {
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
			data[(i- test_start)*dim + j] = stringToNum<float>(str_d[j]);
		}
	}
}


int main(int argc, char** argv) {
	if (argc != 7) {
		std::cout << argv[0]
			<< " data_file query_file nsg_path search_L search_K result_path"
			<< std::endl;
		exit(-1);
	}
	
	float* data_load = NULL;
	unsigned points_num, dim=512;
	std::string filenames = "E:\\work_space\\source_code\\face_search\\nsg_windows\\face.feature";
	load_face_feature_data_train(filenames, dim, data_load, points_num);
	std::cout << "num:" << points_num << " : " << dim << std::endl;
	
	float* query_load = NULL;
	unsigned query_num, query_dim=512;
	load_face_feature_data_test(filenames, query_dim, query_load, query_num);
	std::cout << "num:" << query_num << " : " << query_dim << std::endl;
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
	efanna2e::IndexNSG index(dim, points_num, efanna2e::L2, nullptr);
	index.Load(argv[3]);

	efanna2e::Parameters paras;
	paras.Set<unsigned>("L_search", L);
	paras.Set<unsigned>("P_search", L);

	auto s = std::chrono::high_resolution_clock::now();
	std::vector<std::vector<unsigned> > res;
	for (unsigned i = 0; i < query_num; i++) {
		std::vector<unsigned> tmp(K);
		index.Search(query_load + i * dim, data_load, K, paras, tmp.data());
		res.push_back(tmp);
		for (size_t j = 0; j < K; j++)
		{
			std::cout << tmp[j] << " ";
		}
	}
	auto e = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = e - s;
	std::cout << "search time: " << diff.count() << "\n";

	save_result(argv[6], res);

	//save_result_text(argv[6], res, labels_train, labels_test);



	system("pause");
	return 0;
}
