#ifdef __APPLE__
#include <netpbm/pam.h>
#else
#include <pam.h>
#endif
#include "funcs.h"
#include "mpi.h"
ImageF * genlpfmask(int rows, int cols){
    double *filter = (double *)malloc(rows*cols*sizeof(double));
    ImageF *filter_img  = (ImageF *)malloc(sizeof(ImageF));
    int id, numproc;
    MPI_Init( int* argc , char*** argv);
    MPI_Comm_rank( MPI_COMM_WORLD , &id);
    MPI_Comm_size( MPI_COMM_WORLD , &numproc);
    for (int i = id; i < rows; i+=numproc)
    {
        for (int j = id; j < cols; j+=numproc)
        {
            if ((i <= floor(rows/4.0) || i >= floor(3*rows/4.0)) && (j <= floor(cols/4.0) || j > 3*floor(cols/4.0))){
                filter[i*rows+j] = 1.0;
            }
            else
                filter[i*rows+j] = 0.0;
        }
    }
    MPI_Finalize();
    filter_img->rows=rows;
    filter_img->cols=cols;
    filter_img->data=filter;
    filter_img->widthStep=rows;
    return filter_img;
}