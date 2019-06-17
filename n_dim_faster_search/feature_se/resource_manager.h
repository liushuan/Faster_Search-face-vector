/*
 * resource_manager.h
 *
 *  Created on: 2016年7月25日
 *      Author: guyadong
 */

#ifndef COMMON_SOURCE_CPP_RESOURCE_MANAGER_H_
#define COMMON_SOURCE_CPP_RESOURCE_MANAGER_H_
#include <thread>
#include <type_traits>
#include <memory>
#include <vector>
#include <queue>
#include <utility>
#include <initializer_list>
#include "assert_macros.h"
#include "threadsafe_queue.h"
#include "threadsafe_unordered_map.h"
#include "raii.h"
namespace gdface {
/*
 * 多线程类空间
 */
inline namespace mt{
/*
* 无资源异常
* 当资源数为0时抛出此异常时
* */
class no_resource_except :public std::logic_error {
public:
	// 继承基类构造函数
	using std::logic_error::logic_error;
	no_resource_except(const std::string &at, const std::string &msg) :logic_error(std::string(at).append(":").append(msg)) {}
	no_resource_except(const std::string &at, const std::exception&e) :logic_error(std::string(at).append(":").append(e.what())) {}
	no_resource_except(const std::string &at, const std::exception&e, const std::string &msg) :logic_error(std::string(at).append(":").append(e.what()).append(":").append(msg)) {}
};
/*
 * 多线程环境共享资源管理类
 * 禁止移动/复制构造函数
 * 禁止移动/复制赋值操作符
 * 所有被管理的资源(R)存放在数组中
 * acquire申请资源,当无资源可用时阻塞
 * release释放资源
 * 同一线程内允许嵌套执行acquire/release
 * acquire/release必须配对使用,否则会造成资源泄漏
 * */
template<typename R>
class resource_manager {
private:
	// 资源索引类型
	using resource_index_type=size_t;
	// 资源队列类型
	using resource_queue_type=threadsafe_queue<resource_index_type>;
public:
	// 返回类型,R为标量类型时直接返回R的值,否则返回引用
	using return_type=typename std::conditional<std::is_scalar<R>::value,R,R&>::type;
	resource_manager()=default;
	/* 禁止复制构造 */
	resource_manager(const resource_manager&)=delete;
	/* 禁止移动构造 */
	resource_manager(resource_manager&&)=delete;
	/* 禁止复制赋值 */
	resource_manager&operator=(const resource_manager&)=delete;
	/* 禁止移动赋值 */
	resource_manager&operator=(resource_manager&&)=delete;
	/*
	* 基本构造函数,
	* 使用迭代器为参数的构造函数,适用所有容器对象,
	* lock_count初始化为0
	* occupy_thread初始化为空
	* free_queue中包含所有资源索引
	* */
	template<typename _InputIterator>
	resource_manager(_InputIterator first, _InputIterator last) :
		resource(first, last)
		,lock_count(resource.size(), typename decltype(lock_count)::value_type(0))
		,occupy_thread()
		,free_queue(make_free_queue(resource.size())){
	}
	/*
	 * 对于类型为整数的资源,提供一个简便的构造函数
	 * count 资源数目
	 * start 整数起始值
	 * 根据这两个参数构建一个start开始count个整数作为资源数组
	 * */
	template<typename Enable=typename std::enable_if<std::is_integral<R>::value>::type>
	resource_manager(size_t count,R start=0):resource_manager(make_vector<R>(count,start)){
		/*std::cout<<"resource:"<<resource.size()<<"{";
		for(size_t i=0;i<resource.size();++i){std::cout<<resource[i]<<",";}
		std::cout<<"}"<<std::endl;
		std::cout<<"lock_count:"<<lock_count.size()<<"{";
		for(size_t i=0;i<lock_count.size();++i){std::cout<<lock_count[i]<<",";}
		std::cout<<"}"<<std::endl;
		std::cout<<"occupy_thread size="<<occupy_thread.size()<<std::endl;
		std::cout<<"free_queue:"<<free_queue->size()<<std::endl;*/
	}
	/*
	 * std::vector类型的资源数组为参数的构造函数
	 * */
	resource_manager(const std::vector<R>& res):resource_manager(res.begin(),res.end()){}
	/*
	 * 使用初始化列表为参数的构造函数
	 * */
	resource_manager(std::initializer_list<R> list):resource_manager(list.begin(),list.end()){}

	virtual ~resource_manager(){
		// 将资源数组清空,如果还有线程请求资源会导致抛出no_resource_except异常
		resource.clear();
	}

	/*
	 * 返回一个自动化的资源管理对象(不可跨线程使用)
	 * raii_var对象构造时会自动申请资源
	 * raii_var对象析构时会自动释放资源
	 * raii_var对象的生命周期必须在当前对象生命周期内,否则在执行资源释放时this指针无效
	 * */
	raii_var<return_type> resource_guard(){
		return raii_var<return_type>(
				[this]{return this->acquire();},
				[this](return_type){this->release();}
			);
	}
private:
	std::vector<R> resource;
	// 占用资源的线程中的加锁计数
	std::vector<size_t> lock_count;
	// 保存每个占用资源的线程id和所占用资源索引的映射,初始为空
	threadsafe_unordered_map<std::thread::id,resource_index_type> occupy_thread;
	// 空闲资源(索引)队列,队列中保存的是资源在resource中的索引,初始为resource全部索引
	std::shared_ptr<resource_queue_type> free_queue;
	template<typename Enable=typename std::enable_if<std::is_integral<R>::value>::type>
	static std::vector<R>
	make_vector(size_t count,R start=0){
		std::vector<R> v(count);
		for(size_t i=0;i<count;++i){
			v[i]=R(i+start);
		}
		return v;
	}
	/*
	 * 创建并初始化资源索引队列,将所有资源索引加入队列
	 * */
	static std::shared_ptr<resource_queue_type>
	make_free_queue(size_t size) {
		using v_type = typename resource_queue_type::value_type;
		std::vector<v_type>v(size);
		//创建索引数组
		for (size_t i = 0; i<size; ++i) { v[i] = v_type(i); }
		return std::make_shared<resource_queue_type>(v.begin(), v.end());;
	}
	/*
	 * 阻塞方式从队列中获取可用的资源
	 * 资源数为0时抛出no_resource_except异常
	 */
	return_type acquire(){
		auto this_thread_id=std::this_thread::get_id();
		resource_index_type resource_index;
		// 当前线程重复加锁时不需要再申请资源,将加lock_cout+1,然后返指定的对象
		if(!occupy_thread.find(this_thread_id,resource_index)){
			// 向空闲队列申请资源
			resource_index=free_queue->wait_and_pop();
			throw_except_if_msg(std::logic_error,0>resource_index||0!=lock_count[resource_index],"invalid resource status");
			// 将申请到的资源索引加入线程占用表,代表当前线程已经使用了这个资源
			occupy_thread.insert({this_thread_id,resource_index});
		}
		++lock_count[resource_index];
		// 资源数目为0时抛出异常
		throw_except_if(no_resource_except,resource.empty());
		return resource[resource_index];
	}
	/*
	 * 释放资源
	 */
	void release(){
		resource_index_type resource_index;
		auto thread_id=std::this_thread::get_id();
		throw_except_if_msg(std::logic_error,
				!occupy_thread.find(thread_id,resource_index)||resource_index>=lock_count.size()||lock_count[resource_index]<=0
				,"invalid acquire/release nest");
		--lock_count[resource_index];
		// 当前线程所有嵌套解锁后才将资源重新加入free_queue
		if(0==lock_count[resource_index]){
			// 从线程map中删除当前线程
			occupy_thread.erase(thread_id);
			// 将释放出来的索引号加入free队列
			free_queue->push(resource_index);
		}
	}
};
}/* namespace mt */
} /* namespace gdface */

#endif /* COMMON_SOURCE_CPP_RESOURCE_MANAGER_H_ */
