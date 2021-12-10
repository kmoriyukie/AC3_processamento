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
  int id, numproc, countn, countm;
  MPI_Init( int* argc , char*** argv);
  MPI_Comm_rank( MPI_COMM_WORLD , &id);
  MPI_Comm_size( MPI_COMM_WORLD , &numproc);
  countn = 0;
  countm = 0;
    for (int k = id; k < in_re->cols; k+=numproc)
    {
        for (int l = id; l < in_re->rows; l+=numproc)
        {
            if(inverse == 1){
                for (int m = id; m < in_re->cols; m+=numproc)
                {
                    for (int n = id; n < in_re->rows; n+=numproc)
                    {
                        transf3[l*in_re->widthStep+k] += in_re->data[n*in_re->widthStep+k]*exp(1*_Complex_I*2*M_PI*(l*n/in_re->rows))/n;
                        transf4[l*in_img->widthStep+k] += in_img->data[n*in_img->widthStep+k]*exp(1*_Complex_I*2*M_PI*(l*n/in_img->rows))/n;  
                         countn++;                    
                    }
                    MPI_Reduce( &transf3[l*in_re->widthStep+k] , &transf3[l*in_re->widthStep+k] , countn , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD); 
                    MPI_Reduce( &transf4[l*in_img->widthStep+k] , &transf4[l*in_img->widthStep+k] , countn , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD); 
                    transf[l*in_re->widthStep + k] += transf3[l*in_re->widthStep+m]*exp(1*_Complex_I*2*M_PI*(k*m/in_re->cols))/m;
                    transf2[l*in_img->widthStep + k] += transf4[l*in_img->widthStep+m]*exp(1*_Complex_I*2*M_PI*(k*m/in_img->cols))/m;
                    countm++;
                }
                MPI_Reduce( &transf[l*in_re->widthStep+k] , &transf[l*in_re->widthStep+k] , countm , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD); 
                MPI_Reduce( &transf2[l*in_img->widthStep+k] , &transf2[l*in_img->widthStep+k] , countm , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD); 
                out_re->data[l*in_re->widthStep + k] = (1.0/(in_re->cols*in_re->rows))*transf[l*in_re->widthStep+k];
                out_img->data[l*in_re->widthStep + k] = (1.0/(in_img->cols*in_img->rows))*transf2[l*in_img->widthStep+k];
            }
            else{
                for (int m = id; m < in_re->cols; m+=numproc)
                {
                    for (int n = 0; n < in_re->rows; n++)
                    {
                        transf3[l*in_re->widthStep+k] += in_re->data[n*in_re->widthStep+k]*exp(-1*_Complex_I*2*M_PI*(l*n/in_re->rows));
                        transf4[l*in_img->widthStep+k] += in_img->data[n*in_img->widthStep+k]*exp(-1*_Complex_I*2*M_PI*(l*n/in_img->rows)); 
                        countn++;                       
                    }
                    MPI_Reduce( &transf3[l*in_re->widthStep+k] , &transf3[l*in_re->widthStep+k] , countn , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD); 
                    MPI_Reduce( &transf4[l*in_img->widthStep+k] , &transf4[l*in_img->widthStep+k] , countn , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD); 
                    transf[l*in_re->widthStep + k] += transf3[l*in_re->widthStep+m]*exp(-1*_Complex_I*2*M_PI*(k*m/in_re->cols));
                    transf2[l*in_img->widthStep + k] += transf4[l*in_img->widthStep+m]*exp(-1*_Complex_I*2*M_PI*(k*m/in_img->cols));
                    countm++;
                }
                MPI_Reduce( &transf[l*in_re->widthStep+k] , &transf[l*in_re->widthStep+k] , countm , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD); 
                MPI_Reduce( &transf2[l*in_img->widthStep+k] , &transf2[l*in_img->widthStep+k] , countm , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD); 
                out_re->data[l*in_re->widthStep + k] = transf[l*in_re->widthStep+k];
                out_img->data[l*in_re->widthStep + k] = transf2[l*in_img->widthStep+k];
            }
        }
    }
    MPI_Finalize();
}