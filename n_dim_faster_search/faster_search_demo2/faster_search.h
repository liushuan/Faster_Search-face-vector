/*************************************************
Copyright:tievd.com
Author:shuan.liu
Date:2019-01-10
Description:face vector faster search.
**************************************************/
#ifndef  _FASTER_SEARCH_H
#define _FASTER_SEARCH_H

#include<vector>
#include <memory>

/** storage face eigenvector and face_id.*/
struct Face_Info {
	std::vector <float *> datas;
	std::vector <std::string> face_ids;
};
struct Data_FasterSearch {
	const unsigned int dim;
	Face_Info face_info;
	std::vector<unsigned int * >base_indexs;
	std::vector<float * >base_coss;
	int * sort_var_index;
	Data_FasterSearch(const unsigned int dim) :dim(dim) {}
};

class Faster_Search {
public :
	/**
	param@dim : face vector µÄÎ¬¶È.
	param@threhold1: the min cos(theta) value, to judge the same person.
	param@max_length : support max length people.
	*/
	Faster_Search(const int dim, const float threhold1=0.0f, const unsigned int max_length=10000)
		:dim(dim), threhold1(threhold1), MaxLength(max_length){
		data_graph = std::make_shared<Data_FasterSearch>(dim);
		base = new float[dim];
		memset(base, 0.0f, dim * sizeof(float));
	};
	/** step 1 load data */
	void load_data(std::vector<float *> datas, std::vector<std::string> face_ids);
	
	/** step 2 build graph*/
	void build_graph();
	
	/** step 3 search person id 
	param@test:person face eigenvector.
	*/
	std::string search(float *test);
	std::string search2(float *test);

	std::string search_sse_neon(float * test);

	/* add one person 
	param@add_data:person face eigenvector.
	param@face_id: person face id.
	*/
	bool add_one_person(float * add_data, std::string face_id);
	
	/**delete one person 
	param@face_id: person face id.
	*/
	int remove_one_person(std::string face_id);
	
	/*test graph accuracy 
	param@test_data: test face eigenvector.
	*/
	float test_graph_accuracy(std::vector<float *>test_data);
	
	const std::shared_ptr<Data_FasterSearch> get_data_graph() {
		return data_graph;
	}

	bool storage_graph(std::string graph_file);
	bool restorage_graph(std::string graph_file, const unsigned int max_length = 10000);

	bool save_face_data(std::string face_data, std::vector<float *> feature_data_test);
	bool load_face_data(std::string face_data, std::vector<float *> &feature_data_test);

	~Faster_Search();
private:
	std::shared_ptr<Data_FasterSearch> data_graph;
	float * base;
	//Data_FasterSearch *data_faster_search;
	const unsigned int MaxLength;
	const int dim;
	const float threhold1;
	const float threhold2 = 0.1f;
	const int search_use_dim = 20;
};

#endif // ! _FASTER_SEARCH_H


