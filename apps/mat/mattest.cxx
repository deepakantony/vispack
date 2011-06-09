#include <stream.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>

#include "mat/matrix.h"

//test functions det, svd, QR, inversesvd


main(char** argv, int argc){
    VISMatrix  tmp;
    char filename[20];
    cout << "VISMatrix data filename [input]:  ";
    cin >> filename;
    ifstream  datafile(filename);
    if( !datafile ){
        cout << "Unable to open matrix data file: " << filename << endl;
        exit(0);
    }
    datafile >> tmp;
    cout << "tmp" << tmp << endl;
   
    VISMatrix  C, U,W,V;
    C = tmp;
    cout << "Testing LU <called by det()>: " << endl;
    cout << "det = " << C.det() << endl << endl;  //test LU called by det

    C.svd(U, W, V);  //test svd
    cout << "Testing SVD: " << endl << flush;
    cout << "U" << U << endl;
    cout << "W" << W << endl;
    cout << "V" << V << endl;
    cout << "U*W*~V" << U*W*~V << endl;

    VISMatrix Q, R;
    C.QR(Q, R);  //test QR
    cout << "Testing QR: " << endl;
    cout << "Q" << Q << endl;
    cout << "R" << R << endl;
    cout << "Q*R" << Q*R << endl;

    VISMatrix Cinv;
    if(abs(C.det())> 1e-4)  // |det|>0, depends on precision
    {
	cout << "Testing inversesvd: " << endl;
	Cinv = C.inv();   //test inversesvd
	cout << "Cinv" << Cinv << endl;
	cout << "C*Cinv" << C*Cinv << endl;
    }

}




