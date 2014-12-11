#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
#define M_PI 3.14159265


#include "glew.h"
#include "glut.h"
#include "vec4.h"
#include "mat4.h"
#include "Program.h"
#include <windows.h>
#include <time.h>



#define WAVE_UNIT 0

int width = 1024;

int height = 768;

Program *render = NULL;

GLuint frameobj;

GLuint vbo;

GLuint vboIndices;

int vboSize = 0;

float sunTheta = 1.0;//1.0

float sunPhi = -0.5;//-0.5

float cameraHeight = 6.0;

float cameraTheta = 0.0;

float gridSize = 8.0;

float nyquistMin = 1.0;

float nyquistMax = 1.5;

float seaColor[4] = { 0 / 255.0, 0 / 255.0, 255 / 255.0, 0.01 };

float hdrExposure = 0.4;

GLuint waveTex;

int nbWaves = 60;

vec4f *gerstnerwaves = NULL;

float lambdaMin = 0.02;

float lambdaMax = 30.0;

float heightMax = 0.12;

float waveDirection = 2.4;

float waveDispersion = 1.25f;

float sigmaXsq = 0.0;

float sigmaYsq = 0.0;

float meanHeight = 0.0;

float heightVariance = 0.0;

float amplitudeMax = 0.0;


void load()
{
	char* files[2];

	files[0] = "radiance.sh";
	files[1] = "ocean.sh";

	if (render != NULL) {
		delete render;
	}
	render = new Program(2, files, NULL);
	glUseProgram(render->program);
	glUniform1i(glGetUniformLocation(render->program, "wavesSampler"), WAVE_UNIT);



}

void createMesh()
{
	if (vboSize != 0) {
		glDeleteBuffers(1, &vbo);
		glDeleteBuffers(1, &vboIndices);
	}
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);


	float s = 0.5;

	float vmargin = 0.1;
	float hmargin = 0.1;

	vec4f *data = new vec4f[int(ceil(height * (s + vmargin) / gridSize) + 5) * int(ceil(width * (1.0 + 2.0 * hmargin) / gridSize) + 5)];

	int n = 0;
	int nx = 0;
	for (float j = height * s - 0.1; j > -height * vmargin - gridSize; j -= gridSize) {
		nx = 0;
		for (float i = -width * hmargin; i < width * (1.0 + hmargin) + gridSize; i += gridSize) {
			data[n++] = vec4f(-1.0 + 2.0 * i / width, -1.0 + 2.0 * j / height, 0.0, 1.0);
			nx++;
		}
	}

	glBufferData(GL_ARRAY_BUFFER, n * 16, data, GL_STATIC_DRAW);
	delete[] data;

	glGenBuffers(1, &vboIndices);
	glBindBuffer(GL_ARRAY_BUFFER, vboIndices);

	vboSize = 0;
	GLuint *indices = new GLuint[6 * int(ceil(height * (s + vmargin) / gridSize) + 4) * int(ceil(width * (1.0 + 2.0 * hmargin) / gridSize) + 4)];

	int nj = 0;
	for (float j = height * s - 0.1; j > -height * vmargin; j -= gridSize) {
		int ni = 0;
		for (float i = -width * hmargin; i < width * (1.0 + hmargin); i += gridSize) {
			indices[vboSize++] = ni + (nj + 1) * nx;
			indices[vboSize++] = (ni + 1) + (nj + 1) * nx;
			indices[vboSize++] = (ni + 1) + nj * nx;
			indices[vboSize++] = (ni + 1) + nj * nx;
			indices[vboSize++] = ni + (nj + 1) * nx;
			indices[vboSize++] = ni + nj * nx;
			ni++;
		}
		nj++;
	}

	glBufferData(GL_ARRAY_BUFFER, vboSize * 4, indices, GL_STATIC_DRAW);
	delete[] indices;
}



void generateWaves()
{

	float ratio = pow(lambdaMax / lambdaMin, (float)1 / (nbWaves - 1));
	sigmaXsq = 0.0;
	sigmaYsq = 0.0;
	meanHeight = 0.0;
	if (gerstnerwaves != NULL) {
		delete[] gerstnerwaves;
	}
	gerstnerwaves = new vec4f[nbWaves];
	for (int i = 0; i < nbWaves; ++i) {
		float wavelen = lambdaMin* pow(ratio, i);
		float theta = rand();
		float k = 2.0f * M_PI / wavelen;
		float omega = sqrt(9.81f * k);
		float amp;
		float omega0 = 9.81 / 10.0;
		amp = (8.1*9.81*9.81 / 1000) / pow(omega, 5) * exp(-0.74*pow(omega0 / omega, 4));
		amp *= 0.5 *sqrt(2 * M_PI*9.81 / wavelen);
		amp = 3 * heightMax*sqrt(amp);
		gerstnerwaves[i].x = amp;
		gerstnerwaves[i].y = omega;
		gerstnerwaves[i].z = k * cos(theta);
		gerstnerwaves[i].w = k * sin(theta);
		sigmaXsq += pow(cos(theta), 2.0f) * (1.0 - sqrt(1.0 - k * k * amp * amp));
		sigmaYsq += pow(sin(theta), 2.0f) * (1.0 - sqrt(1.0 - k * k * amp * amp));
		meanHeight -= k * amp * amp * 0.5f;
	}
	glActiveTexture(GL_TEXTURE0 + WAVE_UNIT);
	glBindTexture(GL_TEXTURE_1D, waveTex);
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA32F_ARB, nbWaves, 0, GL_RGBA, GL_FLOAT, gerstnerwaves);
}



double time()
{

	__int64 time;
	__int64 cpuFrequency;
	QueryPerformanceCounter((LARGE_INTEGER*)&time);
	QueryPerformanceFrequency((LARGE_INTEGER*)&cpuFrequency);
	return time / double(cpuFrequency);

}

void displayCallback()
{

	glClearColor(0.0, 0.5, 1.0, 1.0); // sky color is light blue
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	vec4f sun = vec4f(sin(sunTheta) * cos(sunPhi), sin(sunTheta) * sin(sunPhi), cos(sunTheta), 0.0);
	glPolygonMode(GL_FRONT, GL_FILL);
	glPolygonMode(GL_BACK, GL_FILL);
	glDisable(GL_DEPTH_TEST);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameobj);
	glGenerateMipmapEXT(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	glViewport(0, 0, width, height);
	float ch = cameraHeight - meanHeight;

	mat4f view = mat4f(
		0.0, -1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, -ch,
		-1.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0
		);
	view = mat4f::rotatex(cameraTheta / M_PI * 180.0) * view;

	mat4f proj = mat4f::perspectiveProjection(90.0, float(width) / float(height), 0.1 * ch, 1000000.0 * ch);

	float worldToWind[4];
	worldToWind[0] = cos(waveDirection);
	worldToWind[1] = sin(waveDirection);
	worldToWind[2] = -sin(waveDirection);
	worldToWind[3] = cos(waveDirection);

	float windToWorld[4];
	windToWorld[0] = cos(waveDirection);
	windToWorld[1] = -sin(waveDirection);
	windToWorld[2] = sin(waveDirection);
	windToWorld[3] = cos(waveDirection);
	glUseProgram(render->program);
	glUniformMatrix4fv(glGetUniformLocation(render->program, "screenToCamera"), 1, true, proj.inverse().coefficients());
	glUniformMatrix4fv(glGetUniformLocation(render->program, "cameraToWorld"), 1, true, view.inverse().coefficients());
	glUniformMatrix4fv(glGetUniformLocation(render->program, "worldToScreen"), 1, true, (proj * view).coefficients());
	glUniformMatrix2fv(glGetUniformLocation(render->program, "worldToWind"), 1, true, worldToWind);
	glUniformMatrix2fv(glGetUniformLocation(render->program, "windToWorld"), 1, true, windToWorld);
	glUniform3f(glGetUniformLocation(render->program, "worldCamera"), 0.0, 0.0, ch);
	glUniform3f(glGetUniformLocation(render->program, "worldSunDir"), sun.x, sun.y, sun.z);
	glUniform1f(glGetUniformLocation(render->program, "hdrExposure"), hdrExposure);

	glUniform1f(glGetUniformLocation(render->program, "nbWaves"), nbWaves);
	glUniform1f(glGetUniformLocation(render->program, "heightOffset"), -meanHeight);
	glUniform2f(glGetUniformLocation(render->program, "sigmaSqTotal"), sigmaXsq, sigmaYsq);
	glUniform1f(glGetUniformLocation(render->program, "time"), time());
	glUniform4f(glGetUniformLocation(render->program, "lods"),
		gridSize,
		atan(2.0 / height) * gridSize,
		log(lambdaMin) / log(2.0f),
		(nbWaves - 1.0f) / (log(lambdaMax) / log(2.0f) - log(lambdaMin) / log(2.0f)));
	glUniform1f(glGetUniformLocation(render->program, "nyquistMin"), nyquistMin);
	glUniform1f(glGetUniformLocation(render->program, "nyquistMax"), nyquistMax);
	glUniform3f(glGetUniformLocation(render->program, "seaColor"), seaColor[0], seaColor[1], seaColor[2]);
	glPolygonMode(GL_FRONT, GL_FILL);
	glPolygonMode(GL_BACK, GL_FILL);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIndices);
	glVertexPointer(4, GL_FLOAT, 16, 0);
	glEnableClientState(GL_VERTEX_ARRAY);
	glDrawElements(GL_TRIANGLES, vboSize, GL_UNSIGNED_INT, 0);
	glDisableClientState(GL_VERTEX_ARRAY);
	glutSwapBuffers();
}
int oldx;
int oldy;
bool drag;

void mouseClickFunc(int b, int s, int x, int y)
{
	drag = false;
	if ( b == 0) {
		oldx = x;
		oldy = y;
		drag = true;
	}
}

void mouseMotionFunc(int x, int y)
{
	if (drag) {
		sunPhi += (oldx - x) / 400.0;
		sunTheta += (y - oldy) / 400.0;
		oldx = x;
		oldy = y;
		cout << sunPhi<< " " << sunTheta;

	}
	
}
void reshapeFunc(int x, int y)
{
	width = x;
	height = y;

	glutPostRedisplay();
}



void idleFunc()
{
	glutPostRedisplay();
}
void generateTextures(){

	glActiveTexture(GL_TEXTURE0 + WAVE_UNIT);
	glGenTextures(1, &waveTex);
	glBindTexture(GL_TEXTURE_1D, waveTex);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
}
int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("Ocean Lighting");
	glutCreateMenu(NULL);
	glutDisplayFunc(displayCallback);
	glutReshapeFunc(reshapeFunc);
	glutIdleFunc(idleFunc);
	//glutMouseFunc(mouseClickFunc);
	//glutMotionFunc(mouseMotionFunc);
	glewInit();
	generateTextures();
	generateWaves();
	glGenFramebuffersEXT(1, &frameobj);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameobj);
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	createMesh();
	load();
	glutMainLoop();
}
