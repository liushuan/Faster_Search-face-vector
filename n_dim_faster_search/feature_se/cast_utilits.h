#ifndef COMMON_SOURCE_CPP_CAST_UTILITS_H_
#define COMMON_SOURCE_CPP_CAST_UTILITS_H_
#include <type_traits>
#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>
#include <ctime>
#include <string>
#include <map>
#include <list>
#include <set>
#include <sstream>
/** 类型转换模板函数族，实现通用类型L->R的转换 */
namespace net {
	namespace gdface {
		namespace utils {
			// 相同类型直接转发
			template<typename L>
			L
			cast(L left, typename std::decay<L>::type *) {
				return std::forward<L>(left);
			}
			// 数字类型和enum都使用强制类型转换
			template<typename L, typename R>
			typename std::enable_if < !std::is_same<L, R>::value && std::is_arithmetic<L>::value && std::is_arithmetic<R>::value, R > ::type
			cast(L left, R*right) {
				return (R)left;
			}
			template<typename L, typename R>
			typename std::enable_if< !std::is_same<L, R>::value && std::is_enum<L>::value && std::is_enum<R>::value, R>::type
			cast(L left, R*right) {
				return (R)left;
			}
			// number -> string
			template<typename L>
			typename std::enable_if<std::is_arithmetic<L>::value, std::string>::type
			cast(L left, std::string *right) {
				std::ostringstream ss;
				ss << left;
				return ss.str();
			}
			// string -> number
			template<typename R>
			typename std::enable_if<std::is_arithmetic<R>::value, R>::type
			cast(const std::string &left, R *) {
				std::istringstream ss(left);
				R right;
				ss >> right;
				return right;
			}
			std::tm
			inline cast(uint64_t left, std::tm*right) {
				auto t = cast(left, (std::time_t*)nullptr);
				return *std::gmtime(&t);
			}
			uint64_t
			inline cast(const std::tm &left, uint64_t*right) {
				auto l = left;
				auto time = std::mktime(&l);
				return cast(time, (uint64_t*)nullptr);
			}
			std::string
			inline cast(const std::vector<uint8_t> &left, std::string*right) {
				return std::string(left.data(), left.data() + left.size());
			}
			std::vector<uint8_t>
			inline cast(const std::string &left, std::vector<uint8_t>*right) {
				return std::vector<uint8_t>(left.data(), left.data() + left.size());
			}
			template<typename L, typename R>
			typename std::enable_if<!std::is_same<L, R>::value, std::list<R>>::type
			cast(std::list<L>&&left, std::list<R>*) {
				std::list<R> right;
				std::transform(left.begin(), left.end(), std::back_inserter(right), [](L l)->R {return cast(std::move(l), (R*)nullptr); });
				return std::move(right);
			}
			template<typename L, typename R>
			typename std::enable_if<!std::is_same<L, R>::value, std::list<R>>::type
			cast(const std::list<L>&left, std::list<R>*) {
				std::list<R> right;
				std::transform(left.begin(), left.end(), std::back_inserter(right), [](L l)->R {return cast(l, (R*)nullptr); });
				return std::move(right);
			}
			template<typename L, typename R>
			typename std::enable_if<!std::is_same<L, R>::value, std::set<R>>::type
			cast(std::set<L>&&left, std::set<R>*) {
				std::set<R> right;
				std::transform(left.begin(), left.end(), std::insert_iterator<std::set<R>>(right, right.begin()), [](L l)->R {return cast(std::move(l), (R*)nullptr); });
				return std::move(right);
			}
			template<typename L, typename R>
			typename std::enable_if<!std::is_same<L, R>::value, std::set<R>>::type
			cast(const std::set<L>&left, std::set<R>*) {
				std::set<R> right;
				std::transform(left.begin(), left.end(), std::insert_iterator<std::set<R>>(right, right.begin()), [](L l)->R {return cast(l, (R*)nullptr); });
				return std::move(right);
			}
			template<typename L, typename R>
			typename std::enable_if<!std::is_same<L, R>::value, std::vector<R>>::type
			cast(std::vector<L>&&left, std::vector<R>*) {
				std::vector<R> right;
				std::transform(left.begin(), left.end(), std::insert_iterator<std::vector<R>>(right, right.begin()), [](L l)->R {return cast(std::move(l), (R*)nullptr); });
				return std::move(right);
			}
			template<typename L, typename R>
			typename std::enable_if<!std::is_same<L, R>::value, std::vector<R>>::type
			cast(const std::vector<L>&left, std::vector<R>*) {
				std::vector<R> right;
				std::transform(left.begin(), left.end(), std::insert_iterator<std::vector<R>>(right, right.begin()), [](L l)->R {return cast(l, (R*)nullptr); });
				return std::move(right);
			}
			template<typename KL, typename VL, typename KR, typename VR>
			typename std::enable_if<!std::is_same<KL, KR>::value || !std::is_same<VL, VR>::value, std::map<KR, VR>>::type
			cast(std::map<KL, VL>&&left, std::map<KR, VR>*) {
				std::map<KR, VR> right;
				std::transform(left.begin(), left.end(), std::insert_iterator<std::map<KR, VR>>(right, right.begin()),
					[](std::pair<KL, VL> l)->std::pair<KR, VR> {return std::pair<KR, VR>(cast(std::move(l.first), (KR*)nullptr), cast(std::move(l.second), (VR*)nullptr)); });
				return std::move(right);
			}
			template<typename KL, typename VL, typename KR, typename VR>
			typename std::enable_if<!std::is_same<KL, KR>::value || !std::is_same<VL, VR>::value, std::map<KR, VR>>::type
			cast(const std::map<KL, VL>&left, std::map<KR, VR>*) {
				std::map<KR, VR> right;
				std::transform(left.begin(), left.end(), std::insert_iterator<std::map<KR, VR>>(right, right.begin()),
					[](std::pair<KL, VL> l)->std::pair<KR, VR> {return std::pair<KR, VR>(cast(l.first, (KR*)nullptr), cast(l.second, (VR*)nullptr)); });
				return std::move(right);
			}
		} /* namespace utils */
	} /* namespace gdface */
} /* namespace net */
#endif // !COMMON_SOURCE_CPP_CAST_UTILITS_H_

