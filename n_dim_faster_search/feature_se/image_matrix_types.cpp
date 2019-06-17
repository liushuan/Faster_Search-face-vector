#include  <cstring>
#include <stdexcept>
#include "image_matrix_types.h"
#include "assert_macros.h"
uint32_t fs_get_row_stride(const fs_image_matrix_ptr matrix)
{
	if (nullptr == matrix) {
		return -1;
	}
	// align只取低4位
	uint32_t m = (1U << (matrix->align & 0x0f)) - 1;
	return uint32_t((matrix->width + m)&(~m));
}

uint32_t fs_get_matrix_size(const fs_image_matrix_ptr matrix)
{
	if (nullptr == matrix) {
		return -1;
	}
	auto stride = fs_get_row_stride(matrix);
	return stride * matrix->height * matrix->channels;
}

void fs_fill_channels(const fs_image_matrix_ptr matrix, uint8_t * dst_ptr)
{
	if (nullptr == matrix) {
		return ;
	}
	auto row_stride = fs_get_row_stride(matrix);
	if (nullptr == dst_ptr || nullptr == matrix->pixels) {
		return;
	}
	for (uint8_t ch = 0; ch < matrix->channels; ++ch, dst_ptr += matrix->width*matrix->height) {
		auto dst_offset = dst_ptr;
		auto src_ptr = matrix->pixels;
		for (unsigned int y = 0; y < matrix->height; ++y, src_ptr += row_stride*matrix->channels, dst_offset += matrix->width) {
			for (unsigned int x = 0; x<matrix->width; ++x) {
				dst_offset[x] = src_ptr[x*matrix->channels + ch];
			}
		}
	}
}

int fs_matrix_is_NULL(const fs_image_matrix_ptr matrix)
{
	if (nullptr == matrix) {
		return 1;
	}
	return 0 == matrix->width
		&& 0 == matrix->height
		&& 0 == matrix->channels;
}

int fs_make_matrix(fs_image_matrix_ptr matrix, uint32_t with, uint32_t height, uint8_t channels, FS_COLOR_SPACE color_space, uint8_t align, void * pixels)
{
	if (nullptr != matrix) {
		fs_destruct_matrix(matrix);
		matrix->width = with;
		matrix->height = height;
		matrix->channels = channels;
		matrix->color_space = color_space;
		matrix->align = align;
		matrix->shared = nullptr != pixels;
		matrix->pixels = nullptr;
		uint32_t size = fs_get_matrix_size(matrix);
		if (size > 0) {			
			matrix->pixels = (uint8_t*)(matrix->shared ? pixels : malloc(size));
			// 内存分配失败返回
			return nullptr != matrix->pixels;
		}
	}
	return 0;
}

fs_image_matrix_ptr fs_new_matrix(uint32_t width, uint32_t height, uint8_t channels, FS_COLOR_SPACE color_space, uint8_t align, void * pixels) {
	fs_image_matrix_ptr matrix = (fs_image_matrix_ptr)calloc(1, sizeof(fs_image_matrix));
	fs_make_matrix(matrix,width, height, channels, color_space, align, pixels);
	return matrix;
}

fs_image_matrix_ptr fs_new_matrix_s(uint32_t with, uint32_t height, FS_COLOR_SPACE color_space) {
	auto channels = fs_color_depth(color_space);
	if (channels > 0) {
		return fs_new_matrix(with, height, channels, color_space, 0, nullptr);
	}
	return nullptr;
}

void fs_free_matrix(fs_image_matrix_ptr matrix) {
	fs_destruct_matrix(matrix);
	free(matrix);
}

void fs_destruct_matrix(fs_image_matrix_ptr matrix)
{
	if (nullptr != matrix) {
		if (!matrix->shared) {
			free(matrix->pixels);
		}
	}
}

uint8_t fs_color_depth(FS_COLOR_SPACE color_space) {
	switch (color_space) {
	case FSC_GRAYSCALE:
		return 1;
	case FSC_RGB565:
		return 2;
	case FSC_EXT_RGB:
	case FSC_YCbCr:
	case FSC_RGB:
		return 3;
		break;
	case FSC_EXT_RGBA:
	case FSC_EXT_RGBX:
	case FSC_EXT_ARGB:
	case FSC_EXT_XRGB:
	case FSC_EXT_BGRA:
	case FSC_EXT_BGRX:
	case FSC_EXT_ABGR:
	case FSC_EXT_XBGR:
	case FSC_YCCK:
	case FSC_CMYK:
		return 4;
	default:
		return -1;
	}
}

_fs_image_matrix::_fs_image_matrix(uint32_t width, uint32_t height, uint8_t channels, FS_COLOR_SPACE color_space, uint8_t align, void * pixels) {
	fs_make_matrix(this, width, height, channels, color_space, align, pixels);
}

_fs_image_matrix::_fs_image_matrix(uint32_t width, uint32_t height, FS_COLOR_SPACE color_space) {
	fs_make_matrix(this, width, height, fs_color_depth(color_space), color_space, 0, nullptr);
}

_fs_image_matrix::_fs_image_matrix() :width(0), height(0), channels(0), color_space(FSC_UNKNOWN), align(0), shared(0), pixels(nullptr) {
}

_fs_image_matrix::_fs_image_matrix(const _fs_image_matrix & rv){
	auto b = fs_make_matrix(this, rv.width, rv.height, rv.channels, rv.color_space, rv.align, nullptr);
	throw_except_if_msg(std::logic_error, !b, "copy constructor function error")
	if (nullptr != pixels && nullptr != rv.pixels) {
		std::memcpy(pixels, rv.pixels, get_matrix_size());
	}
}

_fs_image_matrix & _fs_image_matrix::operator=(const _fs_image_matrix & other)
{
	fs_destruct_matrix(this);
	auto b = fs_make_matrix(this, other.width, other.height, other.channels, other.color_space, other.align, nullptr);
	throw_except_if_msg(std::logic_error, !b, "copy constructor function error")
	if (nullptr != pixels && nullptr != other.pixels) {
		std::memcpy(pixels, other.pixels, get_matrix_size());
	}
	return *this;
}

_fs_image_matrix::~_fs_image_matrix() {
	fs_destruct_matrix(this);
}

uint32_t _fs_image_matrix::get_matrix_size() const{
	return fs_get_matrix_size((fs_image_matrix_ptr)this);
}

uint32_t _fs_image_matrix::get_row_stride()const {
	return fs_get_row_stride((fs_image_matrix_ptr)this);
}

bool _fs_image_matrix::is_NULL()const{
	auto b = fs_matrix_is_NULL((fs_image_matrix_ptr)this);
	return b ? true : false;
}

void _fs_image_matrix::fill_channels(uint8_t * dst_ptr) const{
	fs_fill_channels((fs_image_matrix_ptr)this, dst_ptr);
}

_fs_image_matrix _fs_image_matrix::to_gray() const
{
	_fs_image_matrix out;
	auto b = fs_to_gray_image_matrix((fs_image_matrix_ptr)this, std::addressof(out));
	throw_except_if_msg(std::logic_error,!b,"fail fs_to_gray_image_matrix")
	return std::move(out);
}


template<FS_COLOR_SPACE COLOR_SPACE>
typename std::enable_if<COLOR_SPACE == FSC_RGB || COLOR_SPACE == FSC_EXT_RGBA || COLOR_SPACE == FSC_EXT_RGB || COLOR_SPACE == FSC_EXT_RGBX, uint8_t>::type
gray_value(const uint8_t *ptr) {
	// red=ptr[0] green=ptr[1] blud=ptr[2]
	return uint8_t((ptr[0] * 299 + ptr[1] * 587 + ptr[2] * 114 + 500) / 1000);
}
template<FS_COLOR_SPACE COLOR_SPACE>
typename std::enable_if<COLOR_SPACE == FSC_EXT_ARGB || COLOR_SPACE == FSC_EXT_XRGB, uint8_t>::type
gray_value(const uint8_t *ptr) {
	// red=ptr[1] green=ptr[2] blud=ptr[3]
	return uint8_t((ptr[1] * 299 + ptr[2] * 587 + ptr[3] * 114 + 500) / 1000);
}

template<FS_COLOR_SPACE COLOR_SPACE>
typename std::enable_if<COLOR_SPACE == FSC_EXT_BGRA || COLOR_SPACE == FSC_EXT_BGR || COLOR_SPACE == FSC_EXT_BGRX, uint8_t>::type
gray_value(const uint8_t *ptr) {
	// red=ptr[2] green=ptr[1] blud=ptr[0]
	return uint8_t((ptr[2] * 299 + ptr[1] * 587 + ptr[0] * 114 + 500) / 1000);
}
template<FS_COLOR_SPACE COLOR_SPACE>
typename std::enable_if<COLOR_SPACE == FSC_EXT_ABGR || COLOR_SPACE == FSC_EXT_XBGR, uint8_t>::type
gray_value(const uint8_t *ptr) {
	// red=ptr[3] green=ptr[2] blud=ptr[1]
	return uint8_t((ptr[3] * 299 + ptr[2] * 587 + ptr[1] * 114 + 500) / 1000);
}

#define RGB565_MASK_RED 0xF800
#define RGB565_MASK_GREEN 0x07E0
#define RGB565_MASK_BLUE 0x001F
template<FS_COLOR_SPACE COLOR_SPACE>
typename std::enable_if<COLOR_SPACE == FSC_RGB565, uint8_t>::type
gray_value(const uint8_t *ptr) {
	// 这个转换在大端系统存在问题
	auto wPixel = uint16_t(*ptr);
	uint8_t R = (wPixel & RGB565_MASK_RED) >> 11; // 取值范围0-31
	uint8_t G = (wPixel & RGB565_MASK_GREEN) >> 5; // 取值范围0-63
	uint8_t B = wPixel & RGB565_MASK_BLUE; // 取值范围0-31
	return  (R * 299 + G * 587 + B * 114 + 500) / 1000;
}

template<FS_COLOR_SPACE COLOR_SPACE>
typename std::enable_if<COLOR_SPACE == FSC_YCbCr || COLOR_SPACE == FSC_YCCK, uint8_t>::type
gray_value(const uint8_t *ptr) {
	return *ptr;
}

template<FS_COLOR_SPACE COLOR_SPACE>
void convert(const uint8_t*src_ptr, uint8_t*dst_ptr, size_t size, size_t src_step) {
	for (size_t y = 0; y < size; ++y, src_ptr += src_step, ++dst_ptr) {
		*dst_ptr = gray_value<COLOR_SPACE>(src_ptr);
	}
}

int fs_to_gray_image_matrix(const fs_image_matrix_ptr matrix, fs_image_matrix_ptr out) {
	if (nullptr == matrix || nullptr == out) {
		return 0;
	}
	if (FSC_GRAYSCALE == matrix->color_space) {
		// 调用复制构造函数
		out = matrix;
	}
	auto new_size = matrix->get_row_stride()*matrix->height;
	fs_make_matrix(out,matrix->width, matrix->height, 1, FSC_GRAYSCALE, matrix->align, nullptr);
	auto dimbuf = fs_color_depth(matrix->color_space);
	auto src_ptr = matrix->pixels;
	auto dst_ptr = out->pixels;
	switch (matrix->color_space) {
	case FSC_RGB:
		convert<FSC_RGB>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_RGBA:
		convert<FSC_EXT_RGBA>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_RGB:
		convert<FSC_EXT_RGB>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_RGBX:
		convert<FSC_EXT_RGBX>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_ARGB:
		convert<FSC_EXT_ARGB>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_XRGB:
		convert<FSC_EXT_XRGB>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_BGRA:
		convert<FSC_EXT_BGRA>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_BGR:
		convert<FSC_EXT_BGR>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_ABGR:
		convert<FSC_EXT_ABGR>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_EXT_XBGR:
		convert<FSC_EXT_XBGR>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_RGB565:
		convert<FSC_RGB565>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_YCbCr:
		convert<FSC_YCbCr>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	case FSC_YCCK:
		convert<FSC_YCCK>(src_ptr, dst_ptr, new_size, dimbuf);
		break;
	default:
		return 0;
	}
	return 1;
}
