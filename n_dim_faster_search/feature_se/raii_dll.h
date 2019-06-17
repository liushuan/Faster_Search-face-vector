#ifndef COMMON_SOURCE_CPP_RAII_DLL_H_
#define COMMON_SOURCE_CPP_RAII_DLL_H_
namespace gdface {
	/* 用于dll分配的资源T的raii管理类，析构时自动正确释放资源
	* T为资源类型，外部不可修改
	*/
	template<typename T>
	class raii_dll {
	public:
		typedef raii_dll<T> _Self;
		typedef T resource_type;
		/* 默认构造函数 */
		raii_dll() :_resource() {}
		/**
		* res 资源对象
		*/
		explicit raii_dll(const T& res) :
			_resource(res) {
		}
		virtual ~raii_dll() {}
		/* 获取资源引用 */
		const T& get() const noexcept { return _resource; }
		const T& operator*() const noexcept { return get(); }
		const T& operator()() const noexcept { return get(); }
		/* 成员指针引用运算符 */
		const T* operator->()const noexcept { return &get(); }
	private:
		/* 封装的资源对象，外部不可修改 */
		T _resource;
	}; /* raii_dll */
} /* namespace gdface */
#endif // !COMMON_SOURCE_CPP_RAII_DLL_H_

