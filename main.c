
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <SDL.h>
#include <glew.h>

#include "psys.h"


void sdl_panic(const char* msg)
{
	fprintf(stderr, "%s: %s\n", msg, SDL_GetError());
	SDL_Quit();
	exit(EXIT_FAILURE);
}

void panic(const char* msg)
{
	fprintf(stderr, "%s\n", msg);
	SDL_Quit();
	exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{

	//printf("sizeof(struct particle_hash_bucket) = %zu\n", sizeof(struct particle_hash_bucket));
	//return 0;

	SDL_Init(SDL_INIT_VIDEO);

	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

	SDL_Window* w = SDL_CreateWindow(
		"DESTROY THE UNIVERSE",
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		0, 0,
		SDL_WINDOW_FULLSCREEN_DESKTOP | SDL_WINDOW_OPENGL);

	if(w == NULL) {
		sdl_panic("SDL_CreateWindow");
	}

	SDL_GLContext c = SDL_GL_CreateContext(w);

	SDL_GL_MakeCurrent(w, c);

	GLenum err = glewInit();
	if(GLEW_OK != err) {
		panic((const char*)glewGetErrorString(err));
	}

	int width, height;
	SDL_GetWindowSize(w, &width, &height);
	glViewport(0, 0, width, height);
	printf("window: %dx%d\n", width, height);

	glOrtho(-width/2, width/2, height/2, -height/2, 1, 0);
	glDisable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPointSize(2.0f);

	struct psys ps;
	psys_init(&ps);

	int exiting = 0;
	for(;;) {
		SDL_Event e;
		while(SDL_PollEvent(&e)) {
			if(e.type == SDL_KEYDOWN) {
				if(e.key.keysym.sym == SDLK_ESCAPE) {
					exiting = 1;
				}
			}
		}

		if(exiting) break;

		psys_step(&ps);

		glClearColor(0, 0.05, 0.1, 1);
		glClear(GL_COLOR_BUFFER_BIT);

		glColor4f(1,1,1,1);
		psys_draw(&ps);

		SDL_GL_SwapWindow(w);
	}

	SDL_GL_DeleteContext(c);

	SDL_Quit();

	return EXIT_SUCCESS;
}

