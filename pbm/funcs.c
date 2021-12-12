
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
  MPI_Comm_rank( MPI_COMM_WORLD , &id);
  MPI_Comm_size( MPI_COMM_WORLD , &numproc);
    for (int i = id; i < rows; i+=numproc)
    {
        for (int j = 0; j < cols; j++)
        {
            if ((i < rows/4.0 || i >= 3*rows/4.0) && (j < cols/4.0 || j >= 3*cols/4.0)){
                filter[i*rows+j] = 1.0;
            }
            else
                filter[i*rows+j] = 0.0;
        }
    }
    filter_img->rows=rows;
    filter_img->cols=cols;
    filter_img->data=filter;
    filter_img->widthStep=cols;
    return filter_img;
}

// --------------------------------------------------------------------------------------------------

void fti(ImageF * in_re, ImageF * in_img, ImageF * out_re, ImageF * out_img, int inverse){

  double complex transf = 0;
  double complex transf2 = 0;
  
  double complex transf_re = 0;
  double complex transf2_img = 0;
  double complex theta = 0;
  int rows = in_re->rows;
  int cols = in_re->cols;
  int step = in_re->cols;
  int id, numproc;
  int countm = 0;
  MPI_Comm_rank( MPI_COMM_WORLD , &id);
  MPI_Comm_size( MPI_COMM_WORLD , &numproc);
  for (int k = 0; k < cols; k++)
  {
      for (int l = 0; l < rows; l++)
      {
        if(inverse == 0){
            for (int m = 0; m < cols; m++)
            {
                for (int n = id; n < rows; n+=numproc)
                {
                    theta = 2.0*M_PI*((double)l*n/rows);
                    transf2 += (in_re->data[m*step+n]+ _Complex_I*in_img->data[m*step + n])*cexp(_Complex_I*theta);
                    countm++;                      
                }
                MPI_Reduce( trasnf2 , transf2 , countm, MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
                countm = 0;
                theta = 2*M_PI*((double)k*m/cols);
                transf  += transf2*cexp(_Complex_I*theta);
            }
        }
        else{
            for (int m = 0; m < cols; m++)
            {
              for (int n = id; n < rows; n+=numproc)
                {
                    theta = -2.0*M_PI*((double)l*n/rows);
                    transf2 += (in_re->data[m*step+n]+ _Complex_I*in_img->data[m*step + n])*cexp(_Complex_I*theta);
                    countm++;                      
                }
                MPI_Reduce( trasnf2 , transf2 , countm , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD);
                countm = 0;
                theta = -2.0*M_PI*((double)k*m/cols);
                transf  += transf2*cexp(_Complex_I*theta);
            }
            transf /= cols*rows;
        }
        out_re->data[l + step*k] = creal(transf);
        out_img->data[l + step*k] = cimag(transf);
        transf=0;
        transf2=0;    
      }
  }
}
//----------------------------------------------------------------------------------------------
void dofilt(ImageF * in_re, ImageF * in_img, ImageF * mask, ImageF * out_re, ImageF * out_img)
{
  int rows = in_re->rows;
  int cols = in_re->cols;
  double *back_re = (double *)malloc(rows*cols*sizeof(double));
  double *back_img = (double *)malloc(rows*cols*sizeof(double));
  int id, numproc;
  MPI_Comm_rank( MPI_COMM_WORLD , &id);
  MPI_Comm_size( MPI_COMM_WORLD , &numproc);

  for(int i = id; i < rows; i+=numproc)
  {
    for(int j = 0; j < cols; j++)
    {
      back_re[cols*j+i] = in_re->data[cols*j+i]*mask->data[cols*j+i];
      back_img[cols*j+i] = in_img->data[cols*j+i]*mask->data[cols*j+i];
    }
  }
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