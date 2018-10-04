#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List matrix_saida(NumericMatrix matrixIn, int resolutionValue,
    double smoothValue, double gamaValue, double increase, double zoomValue,
    int h) {

    List resultado;
    double a=0, b=0, c=0, p=zoomValue, q=zoomValue, t = 0, s1 = 0, s2 = 0, s3 = 0;

    NumericMatrix vetInPos(h+1,1);
    NumericMatrix matrixOut(resolutionValue,resolutionValue);
    NumericMatrix matrixOutExp(resolutionValue,resolutionValue);
    NumericMatrix matrixOutCtrl(resolutionValue,resolutionValue);
    NumericMatrix vetOut(smoothValue+1,1);

    for(int i=0; i<resolutionValue; i++){
        q = zoomValue;
        for(int j=0; j<resolutionValue; j++){
            for(int m=0; m<h; m++){

            a = (p-matrixIn(m,0))*(p-matrixIn(m,0));
            b = (q-matrixIn(m,1))*(q-matrixIn(m,1));
            c = sqrt(a+b);
            vetInPos(m,0) = c;
            }


            for(int m=0; m<smoothValue; m++){
                vetOut(m,0)=smoothValue+1;
                vetInPos(vetOut(m,0),0) = gamaValue;
                for(int n=0; n<h; n++){
                    for(int l=0; l<m; l++){
                        t = 0;
                        if (n == vetOut(l,0)){
                            t = 1;
                        }
                    }
                    if (((vetInPos(n,0) < vetInPos(vetOut(m,0),0)) & (t!=1))){
                        vetOut(m,0) = n;
                    }
                }
            }
            s1 = 0;
            s2 = 0;
            s3 = 0;
            for(int m=0; m<smoothValue; m++){
                s1=s1+matrixIn(vetOut(m,0),2);
                s2=s2+matrixIn(vetOut(m,0),3);
                s3=s3+matrixIn(vetOut(m,0),4);

            }

            //Performs smoothing by joining values in vetOut
            //divided by the smoothing value

            matrixOut(i,j) = s1/smoothValue;
            matrixOutExp(i,j) = s2/smoothValue;
            matrixOutCtrl(i,j) = s3/smoothValue;

            q = q+increase;
        }
        p = p+increase;
    }
    resultado["m1"] = matrixOut;
    resultado["m2"] = matrixOutExp;
    resultado["m3"] = matrixOutCtrl;

    return resultado;
}
