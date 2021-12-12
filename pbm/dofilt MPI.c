#ifdef __APPLE__
#include <netpbm/pam.h>
#else
#include <pam.h>
#endif
#include "funcs.h"
#include "mpi.h"
void dofilt(ImageF * in_re, ImageF * in_img, ImageF * mask, ImageF * out_re, ImageF * out_img)
{
  int rows = in_re->rows;
  int cols = in_re->cols;
  double *back_re = (double *)malloc(rows*cols*sizeof(double));
  double *back_img = (double *)malloc(rows*cols*sizeof(double));
  int id, numproc;
  MPI_Comm_rank( MPI_COMM_WORLD , &id);
  MPI_Comm_size( MPI_COMM_WORLD , &numproc);
  out_re->rows = in_re->rows;
  out_re->cols = in_re->cols;
  out_re->widthStep = in_re->widthStep;

  out_img->rows = in_img->rows;
  out_img->cols = in_img->cols;
  out_img->widthStep = in_img->widthStep;

  for(int j = id; j < cols; j+=numproc)
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