
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
    
    for (int i = 0; i < rows; i++)
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
    filter_img->widthStep=rows;
    return filter_img;
}

// --------------------------------------------------------------------------------------------------

void fti(ImageF * in_re, ImageF * in_img, ImageF * out_re, ImageF * out_img, int inverse){

  double *transf = (double *) malloc(in_re->rows*in_re->cols*sizeof(double));
  double *transf2 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));
  double *transf3 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));
  double *transf4 = (double *) malloc(in_img->rows*in_img->cols*sizeof(double));

    for (int k = 0; k < in_re->cols; k++)
    {
        for (int l = 0; l < in_re->rows; l++)
        {
            if(inverse == 1){
                for (int m = 0; m < in_re->cols; m++)
                {
                    for (int n = 0; n < in_re->rows; n++)
                    {
                        transf3[l*in_re->widthStep+k] += in_re->data[n*in_re->widthStep+k]*exp(1*_Complex_I*2*M_PI*(l*n/in_re->rows))/n;
                        transf4[l*in_img->widthStep+k] += in_img->data[n*in_img->widthStep+k]*exp(1*_Complex_I*2*M_PI*(l*n/in_img->rows))/n;                        
                    }
                    transf[l*in_re->widthStep + k] += transf3[l*in_re->widthStep+m]*exp(1*_Complex_I*2*M_PI*(k*m/in_re->cols))/m;
                    transf2[l*in_img->widthStep + k] += transf4[l*in_img->widthStep+m]*exp(1*_Complex_I*2*M_PI*(k*m/in_img->cols))/m;
                }
                out_re->data[l*in_re->widthStep + k] = (1.0/(in_re->cols*in_re->rows))*transf[l*in_re->widthStep+k];
                out_img->data[l*in_re->widthStep + k] = (1.0/(in_img->cols*in_img->rows))*transf2[l*in_img->widthStep+k];
            }
            else{
                for (int m = 0; m < in_re->cols; m++)
                {
                    for (int n = 0; n < in_re->rows; n++)
                    {
                        transf3[l*in_re->widthStep+k] += in_re->data[n*in_re->widthStep+k]*exp(-1*_Complex_I*2*M_PI*(l*n/in_re->rows));
                        transf4[l*in_img->widthStep+k] += in_img->data[n*in_img->widthStep+k]*exp(-1*_Complex_I*2*M_PI*(l*n/in_img->rows));                       
                    }
                    transf[l*in_re->widthStep + k] += transf3[l*in_re->widthStep+m]*exp(-1*_Complex_I*2*M_PI*(k*m/in_re->cols));
                    transf2[l*in_img->widthStep + k] += transf4[l*in_img->widthStep+m]*exp(-1*_Complex_I*2*M_PI*(k*m/in_img->cols));
                }
                out_re->data[l*in_re->widthStep + k] = transf[l*in_re->widthStep+k];
                out_img->data[l*in_re->widthStep + k] = transf2[l*in_img->widthStep+k];
            }
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
//----------------------------------------------------------------------------------------
ImageF newImageF(int rows, int cols, int widthStep){
  ImageF *img  = (ImageF *)malloc(sizeof(ImageF));
  img->cols = cols;
  img->rows = rows;
  img->widthStep = cols;
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