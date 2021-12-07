
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

// --------------------------------------------------------------------------------------------------

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
        transf[j*R+i] += in_re->(*data)[j][i]*exp(1*_COMPLEX_I*2*M_PI*(i*i/R));
      }
      out_re->(*data)[j][i] += (1/R*C)*transf[j][i]*exp(1*_COMPLEX_I*2*M_PI*(j*j/C));
    }

    for (int i = 0; i < C_img; i++)
      {
        for (int j = 0; j < R_img; j++)
        {
          transf2[j][i] += in_img->(*data)[j][i]*exp(1*_COMPLEX_I*2*M_PI*(i*i/R_img));
        }
        out_img->(*data)[j][i] += (1/R*C)*transf2[j][i]*exp(1*_COMPLEX_I*2*M_PI*(j*j/C_img)); 
      }

    
  }
  else
  {
    /* Calcular DFT */
    for (int i = 0; i < C; i++)
    {
      for (int j = 0; j < R; j++)
      {
        transf[j][i] += in_re->(*data)[j][i]*exp(-1*_COMPLEX_I*2*M_PI*(i*i/R));
      }
      out_re->(*data)[j][i] += transf[j][i]*exp(-1*_COMPLEX_I*2*M_PI*(j*j/C));
    }

    for (int i = 0; i < C_img; i++)
      {
        for (int j = 0; j < R_img; j++)
        {
          transf2[j][i] += in_img->(*data)[j][i]*exp(-1*_COMPLEX_I*2*M_PI*(i*i/R_img));
        }
        out_img->(*data)[j][i] += transf2[j][i]*exp(-1*_COMPLEX_I*2*M_PI*(j*j/C_img)); 
      }
  }
}
//----------------------------------------------------------------------------------------------
void dofilt(ImageF * in_re, ImageF * in_img, ImageF * mask, ImageF * out_re, ImageF * out_img)
{
  int R = in_re->rows;
  int C = in_re->cols;
  double *get = (double *)malloc(rows*cols*sizeof(double)); 
  double *get_img = (double *)malloc(rows*cols*sizeof(double));
  double *get_mask = (double *)malloc(rows*cols*sizeof(double));
  double *back_re = (double *)malloc(rows*cols*sizeof(double));
  double *back_img = (double *)malloc(rows*cols*sizeof(double));
  /*double get[R][C];
  double get_img[R_img][C_img];
  double get_mask[R_mask][C_mask];*/

  for(int i = 0; i<C; i++)
  {
    for(int j = 0; j<R; j++)
    {
      get[R*i+j] = in_re->(*data)[R*i+j];
      get_img[R*i+j] = in_img->(*data)[R*i+j];
      get_mask[R*i+j] = mask->(*data)[R*i+j];
      back_re[R*i+j] = get[R*i+j]*get_mask[R*i+j];
      back_img[R*i+j] = get_img[R*i+j]*get_mask[R*i+j];
    }
  }
  out_re->data = back_re;
  out_img->data = back_img;
}