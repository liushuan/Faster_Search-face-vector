/* 
    Copyright (C) 2013,2014 Wei Dong <wdong@wdong.org>. All Rights Reserved.
*/
#ifndef KGRAPH_VALUE_TYPE
#define KGRAPH_VALUE_TYPE float
#endif

#include <time.h>
#include <cctype>
#include <random>
#include <iomanip>
#include <type_traits>
#include <boost/timer/timer.hpp>
#include <boost/tr1/random.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "kgraph.h"
#include "kgraph-data.h"
#include <opencv2\opencv.hpp>
#include <iostream>
using namespace std;
using namespace boost;
using namespace boost::timer;
using namespace kgraph;

namespace po = boost::program_options; 

typedef KGRAPH_VALUE_TYPE value_type;

static float getAngleDistance(float* feature1, float* feature2, int dim) {
	int length = dim;
	float distanceCos = 0;
	float nor_feature1 = 0, nor_feature2 = 0;
	for (int j = 0; j < length; j++)
	{
		distanceCos += feature1[j] * feature2[j];
		nor_feature1 += pow(feature1[j], 2);
		nor_feature2 += pow(feature2[j], 2);
	}
	distanceCos = distanceCos / (sqrt(nor_feature1)*sqrt(nor_feature2));
	return distanceCos;
}

int main (int argc, char *argv[]) {
    string data_path;
    string output_path;
    KGraph::IndexParams params;
    unsigned D;
    unsigned skip;
    unsigned gap;
    unsigned synthetic;
    float noise;

    bool lshkit = true;

    po::options_description desc_visible("General options");
    desc_visible.add_options()
    ("help,h", "produce help message.")
    ("version,v", "print version information.")
    ("data", po::value(&data_path), "input path")
    ("output", po::value(&output_path), "output path")
    (",K", po::value(&params.K)->default_value(default_K), "number of nearest neighbor")
    ("controls,C", po::value(&params.controls)->default_value(default_controls), "number of control pounsigneds");

    po::options_description desc_hidden("Expert options");
    desc_hidden.add_options()
    ("iterations,I", po::value(&params.iterations)->default_value(default_iterations), "")
    (",S", po::value(&params.S)->default_value(default_S), "")
    (",R", po::value(&params.R)->default_value(default_R), "")
    (",L", po::value(&params.L)->default_value(default_L), "")
    ("delta", po::value(&params.delta)->default_value(default_delta), "")
    ("recall", po::value(&params.recall)->default_value(default_recall), "")
    ("prune", po::value(&params.prune)->default_value(default_prune), "")
    ("reverse", po::value(&params.reverse)->default_value(default_reverse), "")
    ("noise", po::value(&noise)->default_value(0), "noise")
    ("seed", po::value(&params.seed)->default_value(default_seed), "")
    ("dim,D", po::value(&D), "dimension, see format")
    ("skip", po::value(&skip)->default_value(0), "see format")
    ("gap", po::value(&gap)->default_value(0), "see format")
    ("raw", "read raw binary file, need to specify D.")
    ("synthetic", po::value(&synthetic)->default_value(0), "generate synthetic data, for performance evaluation only, specify number of points")
    ("l2norm", "l2-normalize data, so as to mimic cosine similarity");

    po::options_description desc("Allowed options");
    desc.add(desc_visible).add(desc_hidden);

    po::positional_options_description p;
    p.add("data", 1);
    p.add("output", 1);

    po::variables_map vm; 
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm); 

    if (vm.count("raw") == 1) {
        lshkit = false;
    }

    if (vm.count("version")) {
        std::cout << "KGraph version " << KGraph::version() << endl;
        return 0;
    }

    if (vm.count("help")
            || (synthetic && (vm.count("dim") == 0 || vm.count("data")))
            || (!synthetic && (vm.count("data") == 0 || (vm.count("dim") == 0 && !lshkit)))) {
        std::cout << "Usage: index [OTHER OPTIONS]... INPUT [OUTPUT]" << endl;
        std::cout << desc_visible << endl;
		std::cout << desc_hidden << endl;
        return 0;
    }

    if (params.S == 0) {
        params.S = params.K;
    }

    if (lshkit && (synthetic == 0)) {   // read dimension information from the data file
        static const unsigned LSHKIT_HEADER = 3;
        ifstream is(data_path.c_str(), ios::binary);
        unsigned header[LSHKIT_HEADER]; /* entry size, row, col */
        is.read((char *)header, sizeof header);
        BOOST_VERIFY(is);
        BOOST_VERIFY(header[0] == sizeof(value_type));
        is.close();
        D = header[2];
        skip = LSHKIT_HEADER * sizeof(unsigned);
        gap = 0;
    }

	std::vector<double> labels;
    Matrix<value_type> data;

    if (synthetic) {
        if (!std::is_floating_point<value_type>::value) {
            throw std::runtime_error("synthetic data not implemented for non-floating-point values.");
        }
        data.resize(synthetic, D);
        std::cout << "Generating synthetic data..." << endl;
        default_random_engine rng(params.seed);
        uniform_real_distribution<double> distribution(-1.0, 1.0);
        data.zero(); // important to do that
        for (unsigned i = 0; i < synthetic; ++i) {
            value_type *row = data[i];
            for (unsigned j = 0; j < D; ++j) {
                row[j] = distribution(rng);
            }
        }
    }
    else {
		//D = 128;
		//gap = 4;
		//skip = 0;
		//std::cout << data_path <<" D:"<<D <<" skip:"<<skip<<" gap:"<<gap<< std::endl;
        //data.load(data_path, D, skip, gap);

		//std::string filenames = "D:\\work_space\\work_code\\nsg\\source_code\\kgraph\\mnist\\train-images.idx3-ubyte";
		//data.load_mnist_data(filenames, 28*28);

		std::string filenames = "E:\\work_space\\source_code\\face_search\\nsg_windows\\face.feature";
		data.load_face_feature_data_train(filenames, 512, labels);
		data.normalize2();
	}
    if (noise != 0) {
        if (!std::is_floating_point<value_type>::value) {
            throw std::runtime_error("noise injection not implemented for non-floating-point value.");
        }
        tr1::ranlux64_base_01 rng;
        double sum = 0, sum2 = 0;
        for (unsigned i = 0; i < data.size(); ++i) {
            for (unsigned j = 0; j < data.dim(); ++j) {
                value_type v = data[i][j];
                sum += v;
                sum2 += v * v;
            }
        }
        double total = double(data.size()) * data.dim();
        double avg2 = sum2 / total, avg = sum / total;
        double dev = sqrt(avg2 - avg * avg);
        std::cout << "Adding Gaussian noise w/ " << noise << "x sigma(" << dev << ")..." << endl;
        std::normal_distribution<double> gaussian(0, noise * dev);
        for (unsigned i = 0; i < data.size(); ++i) {
            for (unsigned j = 0; j < data.dim(); ++j) {
                data[i][j] += gaussian(rng);
            }
        }
    }
    if (vm.count("l2norm")) {
        std::cout << "L2-normalizing data..." << endl;
        data.normalize2();
    }

    MatrixOracle<value_type, metric::l2sqr> oracle(data);
    KGraph::IndexInfo info;
    KGraph *kgraph = KGraph::create(); //(oracle, params, &info);
    {
        auto_cpu_timer timer;
        kgraph->build(oracle, params, &info);
        std::cout << info.stop_condition << endl;
    }
	std::cout << "out_path:" << output_path << std::endl;
    if (output_path.size()) {
        kgraph->save(output_path.c_str(), labels, KGraph::FORMAT_NO_DIST);
    }

	float acc = 0;
	for (size_t i = 0; i < data.size(); i++)
	{
		if (kgraph->found_id(i, (labels[data.getnearDistance(i)]))) {
			acc++;
		}
	}
	std::cout << " acc is:" << acc / data.size() << std::endl;

	delete kgraph;

	system("pause");
    return 0;
}

