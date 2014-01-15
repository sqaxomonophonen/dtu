#ifndef __MUD_H__
#define __MUD_H__

// my useless data

typedef int mud_fd;

mud_fd mud_open(const char* pathname);
void mud_readn(mud_fd fd, void* vbuf, size_t count);
void mud_close(mud_fd fd);

void mud_load_png_rgba(const char* rel, void** data, int* width, int* height);

#endif//__MUD_H__
