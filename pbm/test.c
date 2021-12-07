#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

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


ImageF * genlpfmask(int rows, int cols){
    double *filter = (double *)malloc(rows*cols*sizeof(double));
    ImageF *filter_img  = (ImageF *)malloc(sizeof(ImageF));
    
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if ((i <= floor(rows/4.0) || i >= floor(3*rows/4.0)) && (j <= floor(cols/4.0) || j > 3*floor(cols/4.0))){
                filter[i*rows+j] = 1.0;
            }
            else
                filter[i*rows+j] = 0.0;
        }
    }
    filter_img->rows=rows;
    filter_img->cols=cols;
    filter_img->data=filter;
    filter_img->widthStep=rows;
    return filter_img;
}

void fti(ImageF * in_re, ImageF * in_img, ImageF * out_re, ImageF * out_img, int inverse){
  int R = in_re->rows;
  int C = in_re->cols;

  int R_img = in_img->rows;
  int C_img = in_img->cols;

  double *transf = (double *) malloc(R*C*sizeof(double));
  double *transf2 = (double *) malloc(R_img*C_img*sizeof(double));
  if (inverse == 1)
  {
    /* Calcular DFT inversa */
    for (int i = 0; i < C; i++)
    {
      for (int j = 0; j < R; j++)
      {
        transf[j*in_re->widthStep+i] += in_re->data[j*in_re->widthStep + i]*exp(1*_Complex_I*2*M_PI*(i*in_re->widthStep/R));
      }
      out_re->(*data)[j][i] += (1/R*C)*transf[j][i]*exp(1*_Complex_I*2*M_PI*(j*j/C));
    }

    for (int i = 0; i < C_img; i++)
      {
        for (int j = 0; j < R_img; j++)
        {
          transf2[j][i] += in_img->(*data)[j][i]*exp(1*_Complex_I*2*M_PI*(i*i/R_img));
        }
        out_img->(*data)[j][i] += (1/R*C)*transf2[j][i]*exp(1*_Complex_I*2*M_PI*(j*j/C_img)); 
      }

    
  }
  else
  {
    /* Calcular DFT */
    for (int i = 0; i < C; i++)
    {
      for (int j = 0; j < R; j++)
      {
        transf[j][i] += in_re->(*data)[j][i]*exp(-1*_Complex_I*2*M_PI*(i*i/R));
      }
      out_re->(*data)[j][i] += transf[j][i]*exp(-1*_Complex_I*2*M_PI*(j*j/C));
    }

    for (int i = 0; i < C_img; i++)
      {
        for (int j = 0; j < R_img; j++)
        {
          transf2[j][i] += in_img->(*data)[j][i]*exp(-1*_Complex_I*2*M_PI*(i*i/R_img));
        }
        out_img->(*data)[j][i] += transf2[j][i]*exp(-1*_Complex_I*2*M_PI*(j*j/C_img)); 
      }
  }
}

int main(){
    ImageF img, *img2;
    Image imgout;
    int i,j;
    
    img.rows=128;
    img.cols=128;
    img.widthStep=128; // Largura da linha + padding (se existir)
    
    img.data=(double *)malloc(img.rows*img.cols*sizeof(double));

    for (i=0;i<img.rows;i++)
	for (j=0;j<img.rows;j++)
	    img.data[img.widthStep*i+j]=128.0;

    img2 = genlpfmask(10, 10);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            printf("%g ",img2->data[i*10+j]);
        }
        printf("rows: %d\n", img2->rows);
    }
    
    return 0;
}