#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix SigCoordPiso(NumericMatrix coord, int
    resolutionValue, double gamaValue, double
    increase, double contrastValue,
    double zoomValue, int numberCoord) {

    double r1 = gamaValue, a=0, b=0, c=0, p=zoomValue, q= zoomValue, t = 0;

    NumericMatrix coordPiso(resolutionValue,resolutionValue);

    for(int i=0; i<resolutionValue; i++){
        q = zoomValue;
        for(int j=0; j<resolutionValue; j++){
            for(int m=0; m<numberCoord; m++){
                a = (double) (p-coord(m,0))*(p-coord(m,0));
                b = (double) (q-coord(m,1))*(q-coord(m,1));
                t = (double) a + (double) b;
                c = (double) sqrt((double) t);

                if ((double) r1 > (double) c) {r1 = c;}
            }
            if((double) r1 > (double) contrastValue) {coordPiso(i,j) = 0;} else
            {coordPiso(i,j) = 10;}
            r1 = gamaValue;
            q = q+increase;
        }
        p = p+increase;
    }
    return coordPiso;
}





