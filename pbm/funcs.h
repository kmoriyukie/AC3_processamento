#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>



#define UM_SEC 1000000000L
#define PI 3.14159265359
#define MAX 5

struct Matrix{
  int rows;
  int cols;
  unsigned char * data;
  int widthStep;
};
typedef struct Matrix Image; 

struct MatrixF{
  int rows;
  int cols;
  double * data;
  int widthStep;
};
typedef struct MatrixF ImageF;

Image * loadPBM(char * fname);
void savePBM(char * fname, Image * image);

/* coloque aqui a seguir adeclaração das funções a desenvolver */
ImageF * genlpfmask(int , int);
void dofilt(ImageF * , ImageF * , ImageF * , ImageF * , ImageF * );
void fti(ImageF *, ImageF *, ImageF *, ImageF *, int);
struct timespec SubtracaoTempo(struct timespec Inicio, struct timespec Fim);
void teste(ImageF * , ImageF *);
ImageF newImageF(int rows, int cols, int widthStep);
Image newImage(int rows, int cols, int widthStep);
double * uchar2db(unsigned char *data, int rows, int cols);
unsigned char * db2uchar(double * data, int rows, int cols);