

#include "funcs.h"

int main(){
    Image* imgin;
    ImageF img, img2, img3, img4;
    Image imgout;
    int i,j;

    imgin = loadPBM("quokka.pgm");
    img = newImageF(100, 100, 100);
    img2 = newImageF(img.rows, img.cols, img.widthStep);
    img3 = newImageF(img.rows, img.cols, img.widthStep);
    img4 = newImageF(img.rows, img.cols, img.widthStep);

    // img.data = (double*) imgin->data;
    // img2.data = (double*) imgin->data;
    // img3.data = (double*) imgin->data;
    // img4.data = (double*) imgin->data;
    ImageF *mask = genlpfmask(img.rows, img.cols);

    printf("after mask create\n");
     
    // printf("img2.rows = %d, img2.cols = %d\n", img2.rows, img2.cols);
    // imgout.data=(unsigned char *) malloc(mask->rows * mask->cols*sizeof(unsigned char));
    // for (i=0;i<img.rows;i++){
    //     for (j=0;j<img.rows;j++)
    //         printf("%d ", (unsigned char )mask->data[img.widthStep*i+j]);
    //     printf("\n\n");
    // }
        
    // dofilt(&img, &img2, mask, &img3, &img4);
    printf("after filt");
    
    imgout.rows=mask->rows;
    imgout.cols=mask->cols;
    imgout.widthStep=1;
    imgout.data=(unsigned char *) malloc(mask->rows * mask->cols*sizeof(unsigned char));

    imgout.data = (unsigned char*) mask->data;
    savePBM("img.pbm",&imgout);
}

    
