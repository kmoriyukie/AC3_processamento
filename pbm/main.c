

#include "funcs.h"

int main(){
    Image* imgin;
    ImageF img, img2, img3, img4;
    Image imgout;
    int i,j;

    imgin = loadPBM("quokka.pgm");
    img = newImageF(200, 200, 200);
    img2 = newImageF(img.rows, img.cols, img.widthStep);
    img3 = newImageF(img.rows, img.cols, img.widthStep);
    img4 = newImageF(img.rows, img.cols, img.widthStep);

    img.data = (double*) imgin->data;
    img2.data = (double*) imgin->data;
    img3.data = (double*) imgin->data;
    img4.data = (double*) imgin->data;
    ImageF *mask = genlpfmask(img.rows, img.cols);

    printf("after mask create\n");
     
    // printf("img2.rows = %d, img2.cols = %d\n", img2.rows, img2.cols);

    // dofilt(&img, &img2, mask, &img3, &img4);
    printf("after filt");
    
    imgout.rows=img.rows;
    imgout.cols=img.cols;
    imgout.widthStep=img.widthStep;
    imgout.data=(unsigned char *) malloc(imgout.rows *imgout.cols);
    
    printf("c");
    imgout.data = (unsigned char*) mask->data;
    savePBM("img.pbm",&imgout);
}

    
