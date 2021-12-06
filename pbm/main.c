

#include "funcs.h"

int main(){
    ImageF img, img2;
    Image imgout;
    int i,j;
    
    img.rows=128;
    img.cols=128;
    img.widthStep=128;
    
    img.data=(double *)malloc(img.rows*img.cols*sizeof(double));

    for (i=0;i<img.rows;i++)
	for (j=0;j<img.rows;j++)
	    img.data[img.widthStep*i+j]=128.0;


    img2.rows=128;
    img2.cols=128;
    img2.widthStep=128;

    img2.data=(double *)malloc(img.rows*img.cols*sizeof(double));

    teste(&img,&img2);

    imgout.rows=128;
    imgout.cols=128;
    imgout.widthStep=128;
    imgout.data=(unsigned char *) malloc(imgout.rows *imgout.cols);
    savePBM("img.pbm",&imgout);
}

    
