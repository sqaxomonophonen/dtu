#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>

#include <png.h>

#include "mud.h"
#include "log.h"

mud_fd mud_open(const char* pathname)
{
	// TODO patch pathname to MUD directory?
	int fd = open(pathname, O_RDONLY);
	if(fd == -1) {
		panicf("open(%s): %s", pathname, strerror(errno));
	}
	return (mud_fd) fd;
}

void mud_readn(mud_fd fd, void* vbuf, size_t n)
{
	char* buf = (char*) vbuf;
	while(n > 0) {
		ssize_t n_read = read((int)fd, buf, n);
		if(n_read == -1) {
			if(errno == EINTR) continue;
			panicf("read: %s", strerror(errno));
		}
		n -= n_read;
		buf += n_read;
	}
}


void mud_close(mud_fd fd) {
	int ret = close((int) fd);
	if(ret == -1) {
		panicf("close: %s", strerror(errno));
	}
}

static void user_error_fn(png_structp png_ptr, png_const_charp error_msg)
{
	panicf("libpng error - %s", error_msg);
}

static void user_warning_fn(png_structp png_ptr, png_const_charp warning_msg)
{
	warnf("libpng warning - %s", warning_msg);
}

static void user_read_data_fn(png_structp png_ptr, png_bytep dest, png_size_t length)
{
	mud_fd fd = *((mud_fd*) png_get_io_ptr(png_ptr));
	mud_readn(fd, dest, length);
}

void mud_load_png_rgba(const char* rel, void** data, int* widthp, int* heightp)
{
	mud_fd fd = mud_open(rel);

	png_byte header[8];
	mud_readn(fd, header, 8);
	if(png_sig_cmp(header, 0, 8) != 0) {
		panicf("mud_load_png_rgba: %s is not a PNG\n", rel);
	}

	png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, (png_voidp)0, user_error_fn, user_warning_fn);
	if(png_ptr == NULL) {
		panicf("png_create_read_struct failed for '%s'", rel);
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if(info_ptr == NULL) {
		panicf("png_create_info_struct failed for '%s'", rel);
	}

	png_set_sig_bytes(png_ptr, 8);
	png_set_read_fn(png_ptr, &fd, user_read_data_fn);
	png_read_info(png_ptr, info_ptr);

	int width = png_get_image_width(png_ptr, info_ptr);
	int height = png_get_image_height(png_ptr, info_ptr);

	if(widthp != NULL) *widthp = width;
	if(heightp != NULL) *heightp = height;

	int channels = png_get_channels(png_ptr, info_ptr);
	int bit_depth = png_get_bit_depth(png_ptr, info_ptr);
	int color_type = png_get_color_type(png_ptr, info_ptr);
	int rowbytes = png_get_rowbytes(png_ptr, info_ptr);

	if(bit_depth != 8) {
		panicf("%s: bit depth != 8", rel);
	}

	if(channels != 4 || color_type != PNG_COLOR_TYPE_RGBA) {
		panicf("%s: not RGBA", rel);
	}

	if(data != NULL) {
		*data = malloc(width * height * 4);

		png_bytep *row_pointers = (png_bytep*) malloc(height * sizeof(png_bytep*));
		for(int i = 0; i < height; i++) {
			row_pointers[i] = ((png_bytep) *data) + rowbytes * i;
		}
		png_read_image(png_ptr, row_pointers);
		free(row_pointers);
	}

	// FIXME clean up?
	mud_close(fd);
}

