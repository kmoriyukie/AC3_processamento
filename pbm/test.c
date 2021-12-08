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

  double *transf = (double *) malloc(in_re->rows*in_re->cols*sizeof(double));
  double *transf2 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));

    for (int k = 0; k < in_re->cols; k++)
    {
        for (int l = 0; l < in_re->rows; l++)
        {
            if(inverse == 1){
                for (int m = 0; m < in_re->cols; m++)
                {
                    for (int n = 0; n < in_re->rows; n++)
                    {
                        transf[n*in_re->widthStep+k] += in_re->data[n*in_re->widthStep+k]*exp(1*_Complex_I*2*M_PI*(l*n/in_re->rows));
                        transf2[n*in_img->widthStep+k] += in_img->data[n*in_img->widthStep+k]*exp(1*_Complex_I*2*M_PI*(l*n/in_img->rows));                        
                    }
                    transf[l*in_re->widthStep + m] += transf[l*in_re->widthStep+m]*exp(1*_Complex_I*2*M_PI*(k*m/in_re->cols));
                    transf2[l*in_img->widthStep + m] += transf2[l*in_img->widthStep+m]*exp(1*_Complex_I*2*M_PI*(k*m/in_img->cols));
                }
                out_re->data[l*in_re->widthStep + k] = (1.0/(in_re->cols*in_re->rows))*transf[l*in_re->widthStep+k];
                out_img->data[l*in_re->widthStep + k] = (1.0/(in_img->cols*in_img->rows))*transf[l*in_img->widthStep+k];
                // printf("%f ", out_re->data[l*in_re->widthStep + k]);
            }
            else{
                for (int m = 0; m < in_re->cols; m++)
                {
                    for (int n = 0; n < in_re->rows; n++)
                    {
                        transf[n*in_re->widthStep+k] += in_re->data[n*in_re->widthStep+k]*exp(-1*_Complex_I*2*M_PI*(l*n/in_re->rows));
                        transf2[n*in_img->widthStep+k] += in_img->data[n*in_img->widthStep+k]*exp(-1*_Complex_I*2*M_PI*(l*n/in_img->rows));                        
                    }
                    transf[l*in_re->widthStep + m] += transf[l*in_re->widthStep+m]*exp(-1*_Complex_I*2*M_PI*(k*m/in_re->cols));
                    transf2[l*in_img->widthStep + m] += transf2[l*in_img->widthStep+m]*exp(-1*_Complex_I*2*M_PI*(k*m/in_img->cols));
                }
                out_re->data[l*in_re->widthStep + k] = transf[l*in_re->widthStep+k];
                out_img->data[l*in_re->widthStep + k] = transf[l*in_img->widthStep+k];
            }
        }
    }
}

void dofilt(ImageF * in_re, ImageF * in_img, ImageF * mask, ImageF * out_re, ImageF * out_img)
{
  int rows = in_re->rows;
  int cols = in_re->cols;
  double *back_re = (double *)malloc(rows*cols*sizeof(double));
  double *back_img = (double *)malloc(rows*cols*sizeof(double));

  out_re->rows = in_re->rows;
  out_re->cols = in_re->cols;
  out_re->widthStep = in_re->widthStep;

  out_img->rows = in_img->rows;
  out_img->cols = in_img->cols;
  out_img->widthStep = in_img->widthStep;
  
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < cols; j++)
    {
      back_re[in_re->widthStep*i+j] = in_re->data[in_re->widthStep*i+j]*mask->data[mask->widthStep*i+j];
      back_img[in_img->widthStep*i+j] = in_img->data[in_img->widthStep*i+j]*mask->data[mask->widthStep*i+j];
    }
  }
  out_re->data = back_re;
  out_img->data = back_img;
}


int main(){
    ImageF img, *mask;
    ImageF *img2 = (ImageF *)malloc(sizeof(ImageF));
    ImageF *img3  = (ImageF *)malloc(sizeof(ImageF));
    ImageF *img4  = (ImageF *)malloc(sizeof(ImageF));
    Image imgout;
    int i,j;
    
    img.rows=10;
    img.cols=10;
    img.widthStep=10; // Largura da linha + padding (se existir)
    
    img.data=(double *)malloc(img.rows*img.cols*sizeof(double));

    img2->rows=10;
    img2->cols=10;
    img2->widthStep=10; // Largura da linha + padding (se existir)
    
    img2->data=(double *)malloc(img.rows*img.cols*sizeof(double));


    img3->rows=10;
    img3->cols=10;
    img3->widthStep=10; // Largura da linha + padding (se existir)
    
    img3->data=(double *)malloc(img.rows*img.cols*sizeof(double));

    img4->rows=10;
    img4->cols=10;
    img4->widthStep=10; // Largura da linha + padding (se existir)
    
    img4->data=(double *)malloc(img.rows*img.cols*sizeof(double));

    for (i=0;i<img.rows;i++)
        for (j=0;j<img.rows;j++){
            img.data[img.widthStep*i+j]=128.0;
            img2->data[img2->widthStep*i+j] = 7;
        }

    mask = genlpfmask(10, 10);
    dofilt(&img, img2, mask, img3, img4);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            printf("%03g ", img3->data[i*10+j]);
        }
        printf("row: %d\n", i);
    }
    // fti(img3, img4, &img, img2, 1);
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            printf("%03g ", img4->data[i*10+j]);
        }
        printf("rows: %d\n", img.rows);
    }
    return 0;
}