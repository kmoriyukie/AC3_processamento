

#include "funcs.h"

int main(){
    Image* imgin;
    ImageF img, img2, img3, img4;
    Image imgout;
    int i,j;

    imgin = loadPBM("quokka.pgm");
    int rows = imgin->rows;
    int cols = imgin->cols;
    int step = imgin->widthStep;
    
    img = newImageF(rows, cols, step);
    img2 = newImageF(rows, cols, step);
    img3 = newImageF(rows, cols, step);
    img4 = newImageF(rows, cols, step);
    img.data = uchar2db(imgin->data, rows, cols);
    
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            img2.data[j*img2.widthStep+i] = 0;
        }
    }
    

    img3.data = uchar2db(imgin->data, rows, cols);
    img4.data = uchar2db(imgin->data, rows, cols);
    ImageF *mask = genlpfmask(rows, cols);
    
    unsigned char * data=(unsigned char *) malloc(rows * cols*sizeof(unsigned char));
    printf("\nStarting DFT...\n");
    fti(&img, &img2, &img3, &img4, 1);
    printf("\nFinished DFT! Started Filtering...\n");
    dofilt(&img3, &img4, mask, &img, &img2);
    printf("\n Finished filtering! Starting IDFT!\n");
    fti(&img, &img2, &img3, &img4, 1);
    printf("\nFinished IDFT!\n");

    data = db2uchar(img3.data, rows, cols);
    imgout.rows=rows;
    imgout.cols=cols;
    imgout.widthStep=step;
    imgout.data=(unsigned char *) malloc(rows *cols*sizeof(unsigned char));

    imgout.data = data;
    savePBM("img.pbm",&imgout);
}

    
