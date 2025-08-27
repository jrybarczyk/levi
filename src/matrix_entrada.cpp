#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List matrix_entrada(NumericMatrix coordPiso,
                    NumericMatrix SignalOut,NumericMatrix signalExp,
                    NumericMatrix signalCtrl,NumericMatrix coord,
                    int resolutionValue, double increase, double zoomValue,
                    int numberCoord) {

    List resultado;
    double p=zoomValue, q= zoomValue;
    int h =0;

    NumericMatrix matrixIn(numberCoord+(resolutionValue*resolutionValue),5);
    NumericMatrix bandcoord (numberCoord, 1);

    for(int i=0; i<resolutionValue; i++){
        for(int j=0; j<resolutionValue; j++){
            for(int m=0; m<numberCoord; m++){
                if ((coord(m,0) < q) && (bandcoord(m,0) == 0)) {
                    matrixIn(h,0)=coord(m,0);
                    matrixIn(h,1)=coord(m,1);
                    matrixIn(h,2)=SignalOut(m,0);
                    matrixIn(h,3)=signalExp(m,0);
                    matrixIn(h,4)=signalCtrl(m,0);
                    h=h+1;
                    bandcoord(m,0) = 10;
                } else {
                    if ((coord(m,0) == q) && (coord(m,1) < p) &
                        (bandcoord(m,0) == 0)) {
                        matrixIn(h,0)=coord(m,0);
                        matrixIn(h,1)=coord(m,1);
                        matrixIn(h,2)=SignalOut(m,0);
                        matrixIn(h,3)=signalExp(m,0);
                        matrixIn(h,4)=signalCtrl(m,0);
                        h=h+1;
                        bandcoord(m,0) = 10;
                    }
                }
            }
            if (coordPiso(i,j) == 0) {
                matrixIn(h,0)=q;
                matrixIn(h,1)=p;
                matrixIn(h,2)=0;
                matrixIn(h,3)=0;
                matrixIn(h,4)=0;
                h=h+1;
            }
            p=p+increase;
        }
        p=zoomValue;
        q=q+increase;
    }
    resultado["m1"] = matrixIn;
    resultado["m3"] = h;

    return resultado;
}

