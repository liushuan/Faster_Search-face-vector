#ifndef WDONG_KGRAPH_DATA
#define WDONG_KGRAPH_DATA

#include <opencv2/opencv.hpp>

#include <cmath>
#include <cstring>
#include <malloc.h>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <boost/assert.hpp>
//#include <io.h>

#ifdef __GNUC__
#ifdef __AVX__
#define KGRAPH_MATRIX_ALIGN 32
#else
#ifdef __SSE2__
#define KGRAPH_MATRIX_ALIGN 16
#else
#define KGRAPH_MATRIX_ALIGN 4
#endif
#endif
#endif

#define KGRAPH_MATRIX_ALIGN 4

namespace kgraph {

    /// L2 square distance with AVX instructions.
    /** AVX instructions have strong alignment requirement for t1 and t2.
     */
    extern float float_l2sqr_avx (float const *t1, float const *t2, unsigned dim);
    /// L2 square distance with SSE2 instructions.
    extern float float_l2sqr_sse2 (float const *t1, float const *t2, unsigned dim);
    extern float float_l2sqr_sse2 (float const *, unsigned dim);
    extern float float_dot_sse2 (float const *, float const *, unsigned dim);
    /// L2 square distance for uint8_t with SSE2 instructions (for SIFT).
    extern float uint8_l2sqr_sse2 (uint8_t const *t1, uint8_t const *t2, unsigned dim);

    extern float float_l2sqr (float const *, float const *, unsigned dim);
    extern float float_l2sqr (float const *, unsigned dim);
    extern float float_dot (float const *, float const *, unsigned dim);


    using std::vector;

    /// namespace for various distance metrics.
    namespace metric {
        /// L2 square distance.
        struct l2sqr {
            template <typename T>
            /// L2 square distance.
            static float apply (T const *t1, T const *t2, unsigned dim) {
                float r = 0;
                for (unsigned i = 0; i < dim; ++i) {
                    float v = float(t1[i]) - float(t2[i]);
                    v *= v;
                    r += v;
                }
                return r;
            }

            /// inner product.
            template <typename T>
            static float dot (T const *t1, T const *t2, unsigned dim) {
                float r = 0;
                for (unsigned i = 0; i < dim; ++i) {
                    r += float(t1[i]) *float(t2[i]);
                }
                return r;
            }

            /// L2 norm.
            template <typename T>
            static float norm2 (T const *t1, unsigned dim) {
                float r = 0;
                for (unsigned i = 0; i < dim; ++i) {
                    float v = float(t1[i]);
                    v *= v;
                    r += v;
                }
                return r;
            }
        };

        struct l2 {
            template <typename T>
            static float apply (T const *t1, T const *t2, unsigned dim) {
                return sqrt(l2sqr::apply<T>(t1, t2, dim));
            }
        };
    }

    /// Matrix data.
    template <typename T, unsigned A = KGRAPH_MATRIX_ALIGN>
    class Matrix {
        unsigned col;
        unsigned row;
        size_t stride;
        char *data;

        void reset (unsigned r, unsigned c) {
            row = r;
            col = c;
			std::cout << " T size is:" << sizeof(T) << std::endl;
			//stride = (sizeof(T) * c + A - 1) / A * A;
			stride = sizeof(T) * c;
            /*
            data.resize(row * stride);
            */
            if (data) free(data);
			std::cout << " A:" << A << " row:" << row << " stride:" << stride << std::endl;
            //data = (char *)memalign(A, row * stride); // SSE instruction needs data to be aligned
			//data = (char *)_aligned_malloc(A, row * stride);
			data = (char *)malloc(row * stride);
			//data = new char[r * stride];
            //if (!data) throw runtime_error("memalign");
        }
    public:
        Matrix (): col(0), row(0), stride(0), data(0) {}
        Matrix (unsigned r, unsigned c): data(0) {
            reset(r, c);
        }
        ~Matrix () {
            if (data) free(data);
        }
        unsigned size () const {
            return row;
        }
        unsigned dim () const {
            return col;
        }
        size_t step () const {
            return stride;
        }
        void resize (unsigned r, unsigned c) {
            reset(r, c);
        }

		T* getData() {
			return (T*)data;
		}

        T const *operator [] (unsigned i) const {
            return reinterpret_cast<T const *>(&data[stride * i]);
        }
        T *operator [] (unsigned i) {
            return reinterpret_cast<T *>(&data[stride * i]);
        }
        void zero () {
            memset(data, 0, row * stride);
        }

        void normalize2 () {
#pragma omp parallel for
            for (int i = 0; i < row; ++i) {
                T *p = operator[](i);
                double sum = metric::l2sqr::norm2(p, col);
                sum = std::sqrt(sum);
                for (unsigned j = 0; j < col; ++j) {
                    p[j] /= sum;
                }
            }
        }
        
        void load (const std::string &path, unsigned dim, unsigned skip = 0, unsigned gap = 0) {
            std::ifstream is(path.c_str(), std::ios::binary);
            //if (!is) throw io_error(path);
            is.seekg(0, std::ios::end);
            size_t size = is.tellg();
            size -= skip;
			
            unsigned line = sizeof(T) * dim + gap;

			std::cout <<"size:"<< size << " line:" << line <<" gap:"<<gap<< std::endl;

            unsigned N =  size / line;
			
            reset(N, dim);
            zero();
            is.seekg(0, std::ios::beg);
            for (unsigned i = 0; i < N; ++i) {
				is.seekg(gap, std::ios::cur);
                is.read(&data[stride * i], sizeof(T) * dim);
            }
            //if (!is) throw io_error(path);
        }

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
			ifstream fileA(filePath);
			if (!fileA)
			{
				cout << "没有找到需要读取的  " << filePath << " 请将文件放到指定位置再次运行本程序。" << endl << "  按任意键以退出";
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
		float getAngleDistance(float* feature1, float* feature2, int dim) {
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

		int getnearDistance(int j ) {

			float * f_data = (float*)data;
			float max_cos = -1;
			int face_id = -1;
			unsigned int length = row;
			for (unsigned int i = 0; i < length; ++i)
			{
				if (i == j) {
					continue;
				}
				float dist_cos = getAngleDistance(f_data + i * col, f_data + j*col, col);
				if (dist_cos > max_cos) {
					max_cos = dist_cos;
					face_id = i;
				}
			}
			return face_id;
		}

		void load_face_feature_data_train(std::string filename, int dim, std::vector<double>&labels) {
			std::vector<std::string> data_all = readStringFromFileData(filename);
			int train_size = data_all.size() - 1000;
			int number_of_images = train_size;
			//int n_rows = 1;
			//int n_cols = 512;
			reset(number_of_images, dim);
			zero();
			float *float_data = (float*)data;
			for (size_t i = 0; i < train_size; i++)
			{
				std::vector<string> str_d;
				splitString(data_all[i], str_d, ",");
				/////////////////////////////////////////
				for (size_t j = 0; j < dim; j++)
				{
					float_data[i*dim + j] = stringToNum<float>(str_d[j]);
				}
				labels.push_back(i);
			}

			std::cout <<"1."<< getAngleDistance(float_data, float_data + 9927 * dim, dim) << std::endl;
			std::cout << "2." << getAngleDistance(float_data, float_data + 6419 * dim, dim) << std::endl;
			std::cout << "3." << getAngleDistance(float_data, float_data + 6866 * dim, dim) << std::endl;
			std::cout << "4." << getAngleDistance(float_data, float_data + 7488 * dim, dim) << std::endl;
		}
		void load_face_feature_data_test(std::string filename, int dim) {
			std::vector<std::string> data_all = readStringFromFileData(filename);
			int test_start = data_all.size() - 1000;
			int number_of_images = 1000;
			reset(number_of_images, dim);
			zero();
			float *float_data = (float*)data;
			for (size_t i = test_start; i < data_all.size(); i++)
			{
				std::vector<string> str_d;
				splitString(data_all[i], str_d, ",");
				/////////////////////////////////////////
				for (size_t j = 0; j < dim; j++)
				{
					float_data[(i- test_start)*dim + j] = stringToNum<float>(str_d[j]);
				}
			}
		}

        void load_lshkit (std::string const &path) {
            static const unsigned LSHKIT_HEADER = 3;
            std::ifstream is(path.c_str(), std::ios::binary);
            unsigned header[LSHKIT_HEADER]; /* entry size, row, col */
            is.read((char *)header, sizeof header);
            //if (!is) throw io_error(path);
            //if (header[0] != sizeof(T)) throw io_error(path);
            is.close();
            unsigned D = header[2];
            unsigned skip = LSHKIT_HEADER * sizeof(unsigned);
            unsigned gap = 0;
            load(path, D, skip, gap);
        }

        void save_lshkit (std::string const &path) {
            std::ofstream os(path.c_str(), std::ios::binary);
            unsigned header[3];
            assert(sizeof header == 3*4);
            header[0] = sizeof(T);
            header[1] = row;
            header[2] = col;
            os.write((const char *)header, sizeof(header));
            for (unsigned i = 0; i < row; ++i) {
                os.write(&data[stride * i], sizeof(T) * col);
            }
        }
    };

    /// Matrix proxy to interface with 3rd party libraries (FLANN, OpenCV, NumPy).
    template <typename DATA_TYPE, unsigned A = KGRAPH_MATRIX_ALIGN>
    class MatrixProxy {
        unsigned rows;
        unsigned cols;      // # elements, not bytes, in a row, 
        size_t stride;    // # bytes in a row, >= cols * sizeof(element)
        uint8_t const *data;
    public:
        MatrixProxy (Matrix<DATA_TYPE> const &m)
            : rows(m.size()), cols(m.dim()), stride(m.step()), data(reinterpret_cast<uint8_t const *>(m[0])) {
        }

#ifndef __AVX__
#ifdef FLANN_DATASET_H_
        /// Construct from FLANN matrix.
        MatrixProxy (flann::Matrix<DATA_TYPE> const &m)
            : rows(m.rows), cols(m.cols), stride(m.stride), data(m.data) {
            if (stride % A) throw invalid_argument("bad alignment");
        }
#endif
#ifdef __OPENCV_CORE_HPP__
        /// Construct from OpenCV matrix.
        MatrixProxy (cv::Mat const &m)
            : rows(m.rows), cols(m.cols), stride(m.step), data(m.data) {
            if (stride % A) throw invalid_argument("bad alignment");
        }
#endif
#ifdef NPY_NDARRAYOBJECT_H
        /// Construct from NumPy matrix.
        MatrixProxy (PyArrayObject *obj) {
            if (!obj || (obj->nd != 2)) throw invalid_argument("bad array shape");
            rows = obj->dimensions[0];
            cols = obj->dimensions[1];
            stride = obj->strides[0];
            data = reinterpret_cast<uint8_t const *>(obj->data);
            if (obj->descr->elsize != sizeof(DATA_TYPE)) throw invalid_argument("bad data type size");
            if (stride % A) throw invalid_argument("bad alignment");
            if (!(stride >= cols * sizeof(DATA_TYPE))) throw invalid_argument("bad stride");
        }
#endif
#endif
        unsigned size () const {
            return rows;
        }
        unsigned dim () const {
            return cols;
        }
        DATA_TYPE const *operator [] (unsigned i) const {
            return reinterpret_cast<DATA_TYPE const *>(data + stride * i);
        }
        DATA_TYPE *operator [] (unsigned i) {
            return const_cast<DATA_TYPE *>(reinterpret_cast<DATA_TYPE const *>(data + stride * i));
        }
    };

    /// Oracle for Matrix or MatrixProxy.
    /** DATA_TYPE can be Matrix or MatrixProxy,
    * DIST_TYPE should be one class within the namespace kgraph.metric.
    */
    template <typename DATA_TYPE, typename DIST_TYPE>
    class MatrixOracle: public kgraph::IndexOracle {
        MatrixProxy<DATA_TYPE> proxy;
    public:
        class SearchOracle: public kgraph::SearchOracle {
            MatrixProxy<DATA_TYPE> proxy;
            DATA_TYPE const *query;
        public:
            SearchOracle (MatrixProxy<DATA_TYPE> const &p, DATA_TYPE const *q): proxy(p), query(q) {
            }
            virtual unsigned size () const {
                return proxy.size();
            }
            virtual float operator () (unsigned i) const {
                return DIST_TYPE::apply(proxy[i], query, proxy.dim());
            }
        };
        template <typename MATRIX_TYPE>
        MatrixOracle (MATRIX_TYPE const &m): proxy(m) {
        }
        virtual unsigned size () const {
            return proxy.size();
        }
        virtual float operator () (unsigned i, unsigned j) const {
            return DIST_TYPE::apply(proxy[i], proxy[j], proxy.dim());
        }
        SearchOracle query (DATA_TYPE const *query) const {
            return SearchOracle(proxy, query);
        }
    };

    inline float AverageRecall (Matrix<float> const &gs, Matrix<float> const &result, unsigned K = 0) {
        if (K == 0) {
            K = result.dim();
        }
       /* if (!(gs.dim() >= K)) throw invalid_argument("gs.dim() >= K");
        if (!(result.dim() >= K)) throw invalid_argument("result.dim() >= K");
        if (!(gs.size() >= result.size())) throw invalid_argument("gs.size() > result.size()");*/
        float sum = 0;
        for (unsigned i = 0; i < result.size(); ++i) {
            float const *gs_row = gs[i];
            float const *re_row = result[i];
            // compare
            unsigned found = 0;
            unsigned gs_n = 0;
            unsigned re_n = 0;
            while ((gs_n < K) && (re_n < K)) {
                if (gs_row[gs_n] < re_row[re_n]) {
                    ++gs_n;
                }
                else if (gs_row[gs_n] == re_row[re_n]) {
                    ++found;
                    ++gs_n;
                    ++re_n;
                }
                else {
                    //throw runtime_error("distance is unstable");
                }
            }
            sum += float(found) / K;
        }
        return sum / result.size();
    }


}

#ifndef KGRAPH_NO_VECTORIZE
#ifdef __GNUC__
#ifdef __AVX__
#if 0
namespace kgraph { namespace metric {
        template <>
        inline float l2sqr::apply<float> (float const *t1, float const *t2, unsigned dim) {
            return float_l2sqr_avx(t1, t2, dim);
        }
}}
#endif
#else
#ifdef __SSE2__
namespace kgraph { namespace metric {
        template <>
        inline float l2sqr::apply<float> (float const *t1, float const *t2, unsigned dim) {
            return float_l2sqr_sse2(t1, t2, dim);
        }
        template <>
        inline float l2sqr::dot<float> (float const *t1, float const *t2, unsigned dim) {
            return float_dot_sse2(t1, t2, dim);
        }
        template <>
        inline float l2sqr::norm2<float> (float const *t1, unsigned dim) {
            return float_l2sqr_sse2(t1, dim);
        }
        template <>
        inline float l2sqr::apply<uint8_t> (uint8_t const *t1, uint8_t const *t2, unsigned dim) {
            return uint8_l2sqr_sse2(t1, t2, dim);
        }
}}
#endif
#endif
#endif
#endif



#endif

