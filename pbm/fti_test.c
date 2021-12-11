#ifdef __APPLE__
#include <netpbm/pam.h>
#else
#include <pam.h>
#endif
#include "funcs.h"
#include "mpi.h"

void fti(ImageF * in_re, ImageF * in_img, ImageF * out_re, ImageF * out_img, int inverse){

  double *transf = (double *) malloc(in_re->rows*in_re->cols*sizeof(double));
  double *transf2 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));
  double *transf3 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));
  double *transf4 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));
  double *transf5 = (double *) malloc(in_re->rows*in_re->cols*sizeof(double));
  double *transf6 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));
  double *transf7 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));
  double *transf8 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));
    for (int k = 0; k < in_re->cols; k++)
    {
        for (int l = 0; l < in_re->rows; l++)
        {
            if(inverse == 1){
                for (int m = 0; m < in_re->cols; m++)
                {
                    for (int n = 0; n < in_re->rows; n++)
                    {
                        transf3[l*in_re->widthStep+k] += in_re->data[n*in_re->widthStep+m]*cexp(1*_Complex_I*2*M_PI*(l*n/in_re->rows))/(in_re->rows);
                        transf4[l*in_img->widthStep+k] += in_img->data[n*in_img->widthStep+m]*cexp(1*_Complex_I*2*M_PI*(l*n/in_img->rows))/(in_re->rows); 
                        transf[l*in_re->widthStep + k] += transf3[n*in_re->widthStep+m]*cexp(1*_Complex_I*2*M_PI*(k*m/in_re->cols))/(in_re->cols);
                        transf2[l*in_img->widthStep + k] += transf4[n*in_img->widthStep+m]*cexp(1*_Complex_I*2*M_PI*(k*m/in_img->cols))/(in_re->cols);                    
                    }
                    
                }
                out_re->data[l*in_re->widthStep + k] = transf[l*in_re->widthStep+k];
                out_img->data[l*in_re->widthStep + k] = transf2[l*in_img->widthStep+k];
            }
            else{
                for (int m = 0; m < in_re->cols; m++)
                {
                    for (int n = 0; n < in_re->rows; n++)
                    {
                        transf3[l*in_re->widthStep+k] += in_re->data[n*in_re->widthStep+k]*cexp(-1*_Complex_I*2*M_PI*(l*n/in_re->rows));
                        transf4[l*in_img->widthStep+k] += in_img->data[n*in_img->widthStep+k]*cexp(-1*_Complex_I*2*M_PI*(l*n/in_img->rows));  
                        transf[l*in_re->widthStep + k] += transf3[l*in_re->widthStep+m]*cexp(-1*_Complex_I*2*M_PI*(k*m/in_re->cols));
                        transf2[l*in_img->widthStep + k] += transf4[l*in_img->widthStep+m]*cexp(-1*_Complex_I*2*M_PI*(k*m/in_img->cols));                    
                    }

                }
                out_re->data[l*in_re->widthStep + k] = transf[l*in_re->widthStep+k];
                out_img->data[l*in_re->widthStep + k] = transf2[l*in_img->widthStep+k];
            }
        }
    }
    
}