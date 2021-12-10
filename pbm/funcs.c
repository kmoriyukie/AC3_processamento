
#ifdef __APPLE__
#include <netpbm/pam.h>
#else
#include <pam.h>
#endif
#include "funcs.h"

void teste(ImageF *in, ImageF *out)
{
  int i, j;
  for (i = 0; i < out->rows; i++)
  {
    // ponteiros para linha i de cada imagem
    double *row1 = (in->data) + i * in->cols;
    double *row2 = (out->data) + i * out->cols;

    for (j = 0; j < out->cols; j++)
    {
      row2[j] = row1[j] + sin(j / 24.0 * M_PI) * 50.0;
    }
  }
}

Image *loadPBM(char *fname)
{
  FILE *file;
  struct pam inpam;
  tuple *tuplerow;
  unsigned int row;
  Image *image;
  int aux;

  file = fopen(fname, "r");
  pnm_readpaminit(file, &inpam, /*PAM_STRUCT_SIZE(tuple_type)*/ sizeof(struct pam));

  printf("Reading image\n");
  printf("width=%d,height=%d,depth=%d\n", inpam.width, inpam.height, inpam.depth);

  /* allocating image*/
  image = (Image *)malloc(sizeof(Image));
  image->cols = inpam.width;
  image->rows = inpam.height;
  image->widthStep = image->cols;
  aux = image->cols & 0x3;
  if (aux != 0)
  {
    image->widthStep += 4 - aux;
  }
  image->data = (unsigned char *)malloc(image->widthStep * image->rows);

  tuplerow = pnm_allocpamrow(&inpam);

  for (row = 0; row < inpam.height; row++)
  {
    unsigned int column;
    pnm_readpamrow(&inpam, tuplerow);
    for (column = 0; column < inpam.width; ++column)
    {
      unsigned int plane;
      for (plane = 0; plane < inpam.depth; ++plane)
      {
        image->data[image->widthStep * row + column] = tuplerow[column][plane];
      }
    }
  }

  pnm_freepamrow(tuplerow);
  fclose(file);
  return image;
}

void savePBM(char *fname, Image *image)
{
  FILE *file;
  struct pam outpam;
  tuple *tuplerow;
  unsigned int row;

  int aux;

  file = fopen(fname, "w");
  outpam.file = file;
  outpam.size = sizeof(struct pam);
  outpam.len = sizeof(struct pam);
  outpam.format = RPGM_FORMAT;
  outpam.plainformat = 0;
  outpam.height = image->rows;
  outpam.width = image->cols;
  outpam.depth = 1;
  outpam.maxval = 255;
  strcpy(outpam.tuple_type, PAM_PGM_TUPLETYPE);
  /*  outpam.allocation_depth=0;
  outpam.comment_p="ficha 4 de Arquitecura de computadores 2010";
  */
  pnm_writepaminit(&outpam);

  printf("Writing image\n");

  tuplerow = pnm_allocpamrow(&outpam);

  for (row = 0; row < outpam.height; row++)
  {
    unsigned int column;
    for (column = 0; column < outpam.width; ++column)
    {
      unsigned int plane;
      for (plane = 0; plane < outpam.depth; ++plane)
      {
        tuplerow[column][plane] = image->data[image->widthStep * row + column];
      }
    }
    pnm_writepamrow(&outpam, tuplerow);
  }

  pnm_freepamrow(tuplerow);
  fclose(file);
}



// --------------------------------------------------------------------------------------------------

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

// --------------------------------------------------------------------------------------------------

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
//----------------------------------------------------------------------------------------------
void dofilt(ImageF * in_re, ImageF * in_img, ImageF * mask, ImageF * out_re, ImageF * out_img)
{
  int rows = in_re->rows;
  int cols = in_re->cols;
  double *back_re = (double *)malloc(rows*cols*sizeof(double));
  double *back_img = (double *)malloc(rows*cols*sizeof(double));
  int id, numproc;
  MPI_Init( int* argc , char*** argv);
  MPI_Comm_rank( MPI_COMM_WORLD , &id);
  MPI_Comm_size( MPI_COMM_WORLD , &numproc);

  for(int j = id; j < cols; j+=numproc)
  {
    for(int j = 0; j < cols; j++)
    {
      back_re[in_re->widthStep*i+j] = in_re->data[in_re->widthStep*i+j]*mask->data[mask->widthStep*i+j];
      back_img[in_img->widthStep*i+j] = in_img->data[in_img->widthStep*i+j]*mask->data[mask->widthStep*i+j];
    }
  }
  MPI_Finalize();
  out_re->data = back_re;
  out_img->data = back_img;
}
//----------------------------------------------------------------------------------------
ImageF newImageF(int rows, int cols, int widthStep){
  ImageF *img  = (ImageF *)malloc(sizeof(ImageF));
  img->cols = cols;
  img->rows = rows;
  img->widthStep = cols;
  img->data = (double*)malloc(cols*rows*sizeof(double));
  return *img;
}
//----------------------------------------------------------------------------------------
Image newImage(int rows, int cols, int widthStep){
  Image *img  = (Image *)malloc(sizeof(Image));
  img->cols = cols;
  img->rows = rows;
  img->widthStep = cols;
  return *img;
}
//----------------------------------------------------------------------------------------
double * uchar2db(unsigned char *data, int rows, int cols){
  double * my_double = (double*) malloc(rows*cols*sizeof(double));

  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols; j++)
    {
      my_double[i + j*cols] = (double) data[i + j*cols];
    } 
  }
  
  return my_double;
}
//----------------------------------------------------------------------------------------
unsigned char * db2uchar(double * data, int rows, int cols){
  unsigned char * my_uchar = (unsigned char *) malloc(rows*cols);

  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols; j++)
    {
      my_uchar[i + j*cols] = (unsigned char) data[i+ j*cols];
    }
    
  }
  return my_uchar;

}