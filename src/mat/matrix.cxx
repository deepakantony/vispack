// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: matrix.cxx,v 1.2 2003/02/12 22:49:25 whitaker Exp $

#include <fstream>
#include "mat/matrix.h"
#include "util/mathutil.h"
#define debug 0

#define XWDEBUG
#define USELAPACKLIB

#ifdef USELAPACKLIB
//#include "f2c.h"
#ifdef XWDEBUG
#include "lapack.h"
#else
#define integer int
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#endif

// for inverse SVD
#define __MAT_EPSILON (1.0e-12)

//extern "C" {
//#include "clapack.h"
//}
#endif




VISMatrixData::VISMatrixData(const datatype* data, int size){
  _refcnt = 0;
  if (size==0){
    cout << "Created a matrix of size zero - stupid" << *((char*)NULL) << endl;
    _data=NULL;
  }
  else{
	_data = new datatype[size];
	_size = size;
	for (int n=0;n<size;n++){
	  _data[n]=data[n];
	}
  }
}


VISMatrix::VISMatrix(const VISVector& other)
    :_cols(other._cols),
     _rows(other._rows), _matrixdata(other._matrixdata)
{
    refCntInc();
}

VISMatrix::VISMatrix(){
  _rows = 0; _cols = 0;
  if (debug)
	cout << "VISMatrix::VISMatrix()" << endl;
  _matrixdata = new VISMatrixData();
  refCntInc();
  if (debug)
	cout << "refCnt()" << refCnt() << endl;
}

VISMatrix::VISMatrix(int row, int col){
    initialize(row, col);
}

VISMatrix::VISMatrix(char* filename){
    ifstream file(filename);
    if (!file){
	cout << "VISMatrix::VISMatrix(char* filename)" << endl;
	cout << "File not found:  " << filename << endl;
    }
    int r,c;
    file >> r >> c;
    initialize(r,c);
    int i, j;
    for (i = 0; i < _rows; i++) {
	for (j = 0; j < _cols; j++) {
	    if (!file)
		cerr << "error reading matrix: ran out of data";
	    file >> poke(i,j);
	}
    }
}

void VISMatrix::initialize(int r,int c){
    if ((r==0)||(c==0)){
      cout << "Created a matrix of size zero -- stupid" << *((char*)NULL) << endl;
      _rows=0;
      _cols=0;
      _matrixdata = new VISMatrixData();
    }
    else{
	_rows=r;
	_cols=c;
	_matrixdata = new VISMatrixData(r*c);
	refCntInc();
    }
}

VISMatrix::~VISMatrix()
{
    _rows = 0;
    _cols = 0;
    refCntDec();
    if (refCnt() == 0)
      delete _matrixdata;
}

VISMatrix& VISMatrix::operator=(const VISVector& other)
{
  if ((VISMatrix*)&other == this)
    return(*this);
  if (debug)
	cout << "VISMatrix& VISMatrix::operator=(const VISMatrix& other)" << endl;
  //if (refCnt() > 0)  //if this has data first delete it
  //delete _matrixdata;
  //cout << refCnt() << endl;
  refCntDec();
  if(refCnt() == 0)
    {
      delete _matrixdata;
    }
  _cols = other._cols;
  _rows = other._rows;
  _matrixdata = other._matrixdata;
  if (debug)
	cout << _matrixdata->_data << "=" << other._matrixdata->_data << endl;
  refCntInc();
  return *this;
}


VISMatrix& VISMatrix::operator=(const VISMatrix& other)
{
  if (&other == this)
    return(*this);
  if (debug)
	cout << 
	  "VISMatrix& VISMatrix::operator=(const VISMatrix& other)" << endl;
  //if (refCnt() > 0)  //if this has data first delete it
  //delete _matrixdata;
  //cout << refCnt() << endl;
  refCntDec();
  if(refCnt() == 0)
    {
      delete _matrixdata;
    }
  _cols = other._cols;
  _rows = other._rows;
  _matrixdata = other._matrixdata;
  if (debug)
	cout << _matrixdata->_data << "=" << other._matrixdata->_data << endl;
  refCntInc();
  return *this;
}

VISMatrix& VISMatrix::operator=(datatype s){
  if (debug){
	cout << "VISMatrix& VISMatrix::operator=(datatype s)" << endl;
	cout << _matrixdata->_data << endl;}
  for (int r=0;r<_rows;r++)
	for (int c=0;c<_cols;c++)
	  poke(r,c) = s;
  return *this;
} 

    
const datatype* VISMatrix::operator[](int row) const
{
  if ((row < _rows)&&(row >= 0)) return dataPtr(row*_cols);
    else{
	  cout << "const datatype* VISMatrix::operator[](int row) const" << endl;
	  cout << "You have attempted to querry row: " << row
		   << " which is out of bounds." << endl;
	  exit(0);
    }
  return 0;  // error condition
}


VISMatrix VISMatrix::operator+(const VISMatrix& other) const{
    if (debug){
	  cout << "VISMatrix VISMatrix::operator+(const VISMatrix& other)" << endl;
	  cout << _matrixdata->_data << " " << other._matrixdata->_data << endl;}
    if ((_rows!=other._rows)||(_cols!=other._cols)){
	  cout << "Invalid matrix addition. (Sizes: " << _rows << "X" << _cols;
	  cout << " + " << other._rows << "X" << other._cols << ")" << endl;
	  //	  return 0; // return null
	  return VISMatrix(); // return null
    }
    else{
	  VISMatrix ans(_rows,_cols);
	  for (int i=0; i<_rows; i++)
	    for (int j=0; j<_cols; j++){
		  ans.poke(i,j)=peek(i,j)+other.peek(i,j);
		  //cout << i << " " << j << endl;
		  //cout << peek(i,j)+other.peek(i,j);
		  //cout << ans.peek(i,j) << " ";
	    }
	  //cout << endl;
	  return ans;
    }
}

VISMatrix VISMatrix::operator+(datatype s) const{
    VISMatrix ans(_rows,_cols);

    for (int i=0; i<_rows; i++)
	for (int j=0; j<_cols; j++)
	    ans.poke(i,j)=peek(i,j)+s;
    return ans;
}

VISMatrix VISMatrix::operator-(const VISMatrix& other) const
{
  if (debug)
	{
	  cout << "VISMatrix VISMatrix::operator-(const VISMatrix& other)" << endl;
	  cout << _matrixdata->_data << " " << other._matrixdata->_data << endl;
	}
  if ((_rows!=other._rows)||(_cols!=other._cols))
	{
	  cout << "Invalid matrix subtraction.(Sizes: " << _rows << "X" << _cols;
	  cout << " - " << other._rows << "X" << other._cols << ")" << endl;
	  exit(0);
	  return 0; 
    }
  else
	{
	  VISMatrix ans(_rows,_cols);
	  for (int i=0; i<_rows; i++)
	    for (int j=0; j<_cols; j++)
		  {
			ans.poke(i,j)=peek(i,j)-other.peek(i,j);
		  }
	  return ans;
    }
}

VISMatrix VISMatrix::operator-(datatype s) const
{
  VISMatrix ans(_rows,_cols);
  for (int i=0; i<_rows; i++)
	for (int j=0; j<_cols; j++)
	  ans.poke(i,j)=peek(i,j)-s;
  return ans;
}

VISMatrix VISMatrix::operator-() const
{
  VISMatrix ans(_rows,_cols);
  for (int i=0; i<_rows; i++)
	for (int j=0; j<_cols; j++)
	  ans.poke(i,j)= -1.0f * peek(i,j);
  return ans;
}

VISMatrix VISMatrix::operator*(datatype s) const
{
  VISMatrix ans(_rows,_cols);
  
  for (int i=0; i<_rows; i++)
	for (int j=0; j<_cols; j++)
	  ans.poke(i,j)=peek(i,j)*s;
  return ans;
}

VISMatrix VISMatrix::operator/(datatype s) const
{
  VISMatrix ans(_rows,_cols);
  
  for (int i=0; i<_rows; i++)
	for (int j=0; j<_cols; j++)
	  ans.poke(i,j)=peek(i,j)/s;
  return ans;
}

VISMatrix VISMatrix::operator*(const VISMatrix& mat2) const
{
  if (_cols!=mat2._rows){
	cout << "Invalid multiplication:  " << _rows << "X" << _cols;
	cout << " times " << mat2._rows << "X" << mat2._cols << endl;
	//	exit(0);
	// 	return(*((VISMatrix*)NULL));
	return VISMatrix();
  }
  const VISMatrix mat1 = *this;
  int J = mat1.r(), K = mat1.c(), L = mat2.c();
  int size = J*L;
  if (!(K == mat2.r()))
    return(VISMatrix());
  VISMatrix ret(J, L);
  datatype *ret_buf; 
  int j, k, l;
  const datatype *mat1_buf, *mat2_buf;
  
  ret_buf = ret.dataBufRef();
  for (j = 0; j < size; j++)
    *(ret_buf++) = 0.0f;

  ret_buf = ret.dataBufRef();
  mat1_buf = mat1.dataBuf();
  mat2_buf = mat2.dataBuf();
  for (l = 0; l < L; l++)
    {
      mat1_buf = mat1.dataBuf();
      for (k = 0; k < K; k++)      
	{
	  ret_buf = ret.dataBufRef() + l*J;
	  for (j = 0; j < J; j++)      
	    *(ret_buf++) += (*(mat1_buf++))*(*mat2_buf);
	  mat2_buf++;
	}
    }
  return(ret);

//   else{
// 	VISMatrix ans(_rows,other._cols);
// 	ans=0.0f;
// 	for (int i=0; i<_rows; i++){
// 	  for (int j=0; j<other._cols; j++){
// 		for (int k=0; k<other._rows; k++){
// 		  //cout << query(i,j+k) << "*" << other.query(i+k,j) << "=";
// 		  //cout << ans[i][j] << " ";	    
// 		  ans.poke(i,j)=ans.peek(i,j)+ this->peek(i,k)*other.peek(k,j);
// 		}
// 		//cout << endl;
// 	  }
// 	}
// 	return ans;
//   }
}


VISVector VISMatrix::operator*(const VISVector& other) const{
    if (c()!=other.n()){
	  cout << "Invalid multiplication: " << r() << "X" << c();
	  cout << " times " << other.n() << "X" << 1 << endl;
	  //	exit(0);
	  // 
	  // put this in so that I could catch it in the debugger.  We should throw an exception here.
	  // RTW - 10-28-99
	  // We should never just "exit" for such an error, it makes it impossible to find
	  //
	  //return(*((VISVector*)NULL));
	  return VISVector();  // Return null
    }
    else{
	  VISVector ans(r());
	  ans=0.0f;
	  for (int i=0; i<r(); i++){
	    for (int k=0; k<c(); k++){
		  ans.poke(i)=ans.peek(i)+peek(i,k)*other.peek(k);
		}
	  }
	  return ans;
    }
}

VISMatrix VISMatrix::operator~() const{
  VISMatrix ans(_cols,_rows);
  int i, j;
    for (i = 0; i < _rows; i++)
	  for (j = 0; j < _cols; j++)
	    ans.poke(j,i)=peek(i,j);
    return ans;
}//~transpose

VISMatrix VISMatrix::operator==(datatype s) const{
  VISMatrix ans(_rows,_cols);
  int i, j;
  for (i = 0; i < _rows; i++)
	for (j = 0; j < _cols; j++)
	  if (peek(i,j) == s)
		ans.poke(i,j)=1.0f;
	  else
		ans.poke(i,j)=0.0f;
  return ans;
}

VISMatrix VISMatrix::operator==(const VISMatrix& other) const{
  //problem here?
  int r,c;
  if ((_rows!=other.r()) || (_cols!=other.c())){
	cout << "VISMatrix VISMatrix::operator==(const VISMatrix& other) const" << endl;
	cout << "You attempted to compare matricies of differing sizes" << 
	  endl;
  }
  VISMatrix ans(_rows,_cols);
  for (r=0;r<_rows;r++)
	for (c=0;c<_cols;c++)
	  if (peek(r,c)==other.peek(r,c))
		ans.poke(r,c)=1.0f;
	  else
		ans.poke(r,c)=0.0f;
  return ans;
}

// this prints out the matrix is a nice way (as test) to an ostream object
void VISMatrix::print(ostream *os) const{
#ifdef SGI
  *os << setw(3) << _rows << " " << setw(3) << _cols << endl;
#else //ifdef LINUX
  *os << _rows <<" " <<_cols<<endl;
#endif
  int i, j;
  for (i = 0; i < _rows; i++){
	for (j = 0; j < _cols; j++)
#ifdef SGI
	  *os << setw(7) << setprecision(7)  << peek(i,j) << " ";
#else //ifdef LINUX
	*os<< peek(i,j)<<" ";
#endif
	*os << "\n";
  }
}

void VISMatrix::read(istream *is){
    //initialize this
    int r,c;
    *is >> r >> c;
    initialize(r,c);
    int i, j;
    for (i = 0; i < _rows; i++) {
	for (j = 0; j < _cols; j++) {
	    if (!is)
		cerr << "error reading matrix: ran out of data";
	    *is >> poke(i,j);
	}
    }
}


//Returns a 1-D pointer to matrix data as needed by LAPACK
datatype* VISMatrix::createLAPK() const{
  datatype* lapk;
  lapk = new datatype[_rows*_cols];
  
  const datatype *d_buf = dataBuf();
  for (int i=0; i<_rows*_cols; i++){   //matrix => 2-D array to 1-D
    lapk[i] = (*d_buf++);  //Column-major! order in LAPACK
    //    lapk[i] = peek(i%_rows, i/_rows);
  } 
  return lapk;
}//datatype* VISMatrix::createLAPK()

//The matrix becomes the one corresponding to the LAPACK data
void VISMatrix::becomeLAPK(datatype* lapk, int r, int c)
{
  refCntDec();
  if (refCnt() == 0)
    {
      delete _matrixdata;
    }
  initialize(r,c);
  datatype *d_buf = dataBufRef();
  for (int i=0; i<_rows*_cols; i++)
    {   //1-D array to 2-D => matrix
      //      poke(i%_rows, i/_rows) = lapk[i];  //Column-major! order in LAPACK
      // REALLY?
      *(d_buf++) = lapk[i];  //Column-major! order in LAPACK  ??????? 
    }
}//void VISMatrix::becomeLAPK(datatype*, int, int)

void printLAPK(datatype* lapk, int r, int c){
    for (int i=0; i<r*c; i++){
	cout << lapk[i] << "  ";
    }
    cout << endl;
}//printLAPK(datatype*, int, int)

//Clean up the data declared by createLAPK
void destroyLAPK(datatype* lapk){
 
  if(lapk!=NULL)
    delete[] lapk;
 
}//destroyLAPK(datatype*)


//void VISMatrix::destroy(){
//    _cols=0;
//    _rows=0;
//    if (_matrixdata!=NULL)
//	delete [] _data;
//}


VISVector VISMatrix::vec(int col) const{
    int cnt;
    VISVector temp(_rows);
    for(cnt=0;cnt<_rows;cnt++){
      temp.poke(cnt) = peek(cnt,col);
    }
    return temp;
}

VISVector VISMatrix::vecFromRow(int row) const{
    int cnt;
    VISVector temp(_cols);
    for(cnt=0;cnt<_cols;cnt++){
	temp.poke(cnt) = peek(row,cnt);
    }
    return temp;
}

VISVector VISMatrix::vecFromDiag() const{
    VISVector temp(_cols);
    for(int cnt=0;cnt<_cols;cnt++)
	temp.poke(cnt) = peek(cnt,cnt);
    return temp;
}

VISMatrix VISMatrix::col(int col) const {
    
    int cnt;
    VISMatrix temp(_rows,1);
    for(cnt=0;cnt<_rows;cnt++){
	temp.poke(cnt,0) = peek(cnt,col);
    }
    return temp;
}

VISMatrix VISMatrix::row(int row) const {
    
    int cnt;
    VISMatrix temp(1,_cols);
    for(cnt=0;cnt<_cols;cnt++){
	temp.poke(0,cnt) = peek(row,cnt);
    }
    return temp;
}

VISMatrix VISMatrix::peekROI(int startrow,int endrow,int startcol,int endcol) const{
    if (((startrow >= _rows)||(endrow >= _rows))|| 
        ((startcol >= _cols)||(endcol >= _cols))){
	cout << "Error in VISMatrix VISMatrix::peekROI(int startrow,int endrow," <<
	    "int startcol,int endcol)" << endl;
	cout << "out of bounds: matrix size is (" << _rows << "," << _cols
	     << ")" << endl;
	cout << "startrow=" << startrow << ", endrow=" << endrow << ", "
	     << "startcol=" << startcol << ", endcol=" << endcol << endl;
	exit(0);
    }
	
    VISMatrix ans(endrow-startrow+1,endcol-startcol+1);
    for (int r=startrow;r<=endrow;r++)
	for (int c=startcol;c<=endcol;c++){
	    ans.poke(r-startrow,c-startcol) = peek(r,c);
	}
    return ans;
}

VISMatrix&  VISMatrix::pokeROI(int row,int col, const VISMatrix& roi){
    //pokeROI does not change the size of the matrix being operated on
    int roirows;
    int roicols;
    roirows = roi.r();
    roicols = roi.c();

    if (((row+roirows)>_rows)||((col+roicols)>_cols)){
	cout << "VISMatrix& VISMatrix::pokeROI(int row,int col,VISMatrix roi)"
	     << endl;
	cout << "You have attempted to poke a matrix out of bounds "
	     << "of the matrix being operated on." << endl;
	cout << "Matrix being operated on is unchanged." << endl;
	return *this;
    }

    if (refCnt() > 1){
	refCntDec();
	_matrixdata = new VISMatrixData(dataPtr(),_rows*_cols);
	refCntInc();
    }
    
    int r,c;
    for (r=0;r<roirows;r++)
	for (c=0;c<roicols;c++)
	    poke(row+r,col+c) = roi.peek(r,c);
    return *this ;
}

VISMatrix& VISMatrix::pokeROI(int row,int col, const VISVector& roi){
    //pokeROI does not change the size of the matrix being operated on
    int roirows;
    roirows = roi.n();

    if ((row+roirows)>_rows){
	cout << "VISMatrix& VISMatrix::pokeROI(int row,int col,VISVector roi)"
	     << endl;
	cout << "You have attempted to poke a vector out of bounds "
	     << "of the matrix being operated on." << endl;
	cout << "Matrix being operated on is unchanged." << endl;
	return *this;
    }

    if (refCnt() > 1){
	refCntDec();
	_matrixdata = new VISMatrixData(dataPtr(),_rows*_cols);
	refCntInc();
    }

    int r;
    for (r=0;r<roirows;r++)
      poke(row+r,col) = roi.peek(r);

    return(*this);
}

datatype VISMatrix::normOfCol(int col) const{
    double total=0.0;
    datatype ans;
    for (int row=0;row<r();row++)
	total = total + (peek(row,col)*peek(row,col));
    ans = (datatype)sqrt(total); //_sqrt_d
    return ans;
}

datatype VISMatrix::norm() const{
    double total=0.0;
    datatype ans;
    int row, col;
    for (row=0; row < _rows; row++)
	for (col=0;col < _cols; col++)
	    total = total + pow(peek(row,col),2);
    ans = (datatype)sqrt(total); //_sqrt_d
    return ans;
}

datatype VISMatrix::norm(datatype e) const{
    double total=e*e;
    datatype ans;
    int row, col;
    for (row=0; row < _rows; row++)
	for (col=0;col < _cols; col++)
	    total = total + pow(peek(row,col),2);
    ans = (datatype)sqrt(total); //_sqrt_d
    return ans;
}


datatype VISMatrix::normOfRow(int row) const{
    double total=0.0;
    datatype ans;
    for (int col=0;col<_cols;col++)
	total = total + (peek(row,col)*peek(row,col));
    ans = (datatype)sqrt(total);  //_sqrt_d
    return ans;
}

datatype VISMatrix::max(int &r,int &c) const{
    datatype ans;
    ans = peek(0,0);
    r=0;
    c=0;
    for (int rows=0;rows<_rows;rows++)
	  for (int cols=0;cols<_cols;cols++)
	    if (peek(rows,cols)>ans){
		  ans = peek(rows,cols);
		  r=rows;
		  c=cols;
	    }
    return ans;
}

datatype VISMatrix::min(int &r,int &c) const{
  datatype ans;
  ans = peek(0,0);
  r=0;
  c=0;
  for (int rows=0;rows<_rows;rows++)
	for (int cols=0;cols<_cols;cols++)
	  if (peek(rows,cols)<ans){
		ans = peek(rows,cols);
		r=rows;
		c=cols;
	  }
  return ans;
}

VISVector VISMatrix::mean() const
{
  int rows,cols;
  rows=r();
  cols=c();
  VISVector ans(rows);
  double total;
  for(int row=0;row<rows;row++)
    {
	total=0;
	for(int col=0;col<cols;col++)
	  total=total+peek(row,col);
	ans.poke(row) = (datatype)(total/(datatype)cols);
    }
  return ans;
}

VISVector VISMatrix::meanOfCols() const{
  VISVector ans(c());
  int rows=r();
    int cols=c();
    double total;
    for(int col=0;col<cols;col++){
	total=0;
	for(int row=0;row<rows;row++)
	    total=total+peek(row,col);
	ans.poke(col) = (datatype)(total/(double)rows);
    }
    return ans;
}


datatype& VISMatrix::poke(int r,int c)
{
  if (refCnt() > 1)
	{
	  refCntDec();
	  _matrixdata = new VISMatrixData(dataPtr(),_cols*_rows);
	  refCntInc();
	}
  if((r < _rows)&&(c < _cols))
	{
	  return dataRef(r+c*_rows);
	}
  else
	{
	  cout << "You have attempted to poke row: " << r << " and column: "
		   << c << " which does not exist." << endl;
	  //return(*((datatype*)NULL));
	  return dataRef();
	}
}

VISMatrix VISMatrix::dot(const VISMatrix &other) const{
  //dot(column wise) two matricies that hold
  //column vectors.
  //result(ans) is a scalar for each column
  int row,col;
  VISMatrix ans(_rows,1);
  for (row=0;row<_rows;row++){
	for (col=0;col<_cols;col++){
	  ans.poke(row,0) = ans.peek(1,col) + 
		(peek(row,col) * other.peek(row,col));
	}
  }
  return ans;
}

VISMatrix VISMatrix::cov() const{
    //this routine assumes a the input is a matrix
    //that as sets of numerical observations where
    //each column of *this comprises a set of observations
    //if you have the oposite situation (numerical observations
    //along the rows) simply use (~M).cov() -- i.e. cov of transpose
    int rows=r(),cols=c();
    VISMatrix Y;
    Y = *this-((mean()).mat(cols));
    VISMatrix ans;
    ans = (1.0f/((datatype)cols-1.0f))*Y*~Y;
    return ans;
}

datatype VISMatrix::sum() const 
{
  
  datatype total = 0;
  for (int i = 0; i < _cols; i++)
    for (int j = 0; j < _rows; j++)    
      {
	total += peek(j, i);
      }

  return(total);

}

datatype VISMatrix::trace() const
{

  datatype total = 0;
  int n = VISmin(_rows, _cols);

  for (int i = 0; i < n; i++)
    total += peek(i, i);

  return total;
}

#ifdef USELAPACKLIB

//Compute the determinant of a matrix
datatype VISMatrix::det() const{
 
//    datatype *datatype_tmp;
//    datatype_tmp = new datatype[5*5];

 
    integer  M = _rows;
    integer  N = _cols;
    integer  lda = M;  //leading dimension of A: >=max(1,M)
    integer  ret;      //return value: 0=successful exit
    integer  *ipiv;    //pivot indices array
    datatype *a;      //input/output array



 
//    datatype *datatype_tmp2;
//    datatype_tmp2 = new datatype[5*5];

 
    a = createLAPK();
    integer ldipiv = MIN(M, N);
    ipiv = new integer [ldipiv];

#ifdef VISMATRIXDOUBLE
    dgetrf_(&M, &N, a, &lda, ipiv, &ret); //LU from LAPACK
 #else
    sgetrf_(&M, &N, a, &lda, ipiv, &ret); //LU from LAPACK
#endif

    VISMatrix A(M, M);
    A = 0.0f;

    A.becomeLAPK(a, M, M);
    datatype ans=1.0f;
    for (int cnt=0; cnt<M; cnt++)
	ans = ans*A[cnt][cnt];
    int isign=1; 
    for(int i=0; i<ldipiv; i++){ 
        if(ipiv[i] != i) //interchanges
	    isign = -1 * isign;
    }
    ans = ans*isign;
    delete[] ipiv;
    destroyLAPK(a);
    return ans;
}//det

//Perform the SVD decomposition of a matrix
void VISMatrix::svd(VISMatrix &U, VISMatrix &W, VISMatrix &V) const{
    //A=U*W*transpose(V)
    //The following stuff are needed by SVD of LAPACK//
    integer  M = _rows;
    integer  N = _cols;
    integer S = VISmax(M, N);
    char  ju  = 'A';   //'All' columns of U are returned
    char  jvt = 'A';   //'All' columns of VT are returned
    integer  la = M;   //leading dimension of A: >=max(1,M)
    integer  lu = M;   //leading dimension of U: >= M
    integer  lvt = N;  //leading dimension of VT: >= N
    integer  lwk;      //dimension of the array work
    integer  ret;      //return value: 0=successful exit
    datatype *a, *work;     //input array and work space
    datatype *u, *w, *vt;   //output arrays
    // from the lapack documentation...
    lwk =VISmax(3*VISmin(M, N) + VISmax(M, N), 5*VISmin(M, N));
    work = new datatype [lwk];   
    a = createLAPK();
    w = new datatype [N]; 
    //    w = new datatype [S]; // this should be N but lapack seems to be broken
    vt = new datatype [N*N];
    //vt = new datatype [S*S]; // this should be NxN but lapack seems to be broken
    //  u = new datatype [S*S];
    //u = new datatype [N*N];
    u = new datatype [M*M];

    //printLAPK(a, M, N);
#ifdef VISMATRIXDOUBLE


    //    cout << "about to do svd" << endl;
    dgesvd_(&ju, &jvt, &M, &N, a, &la, w, u, &lu, vt, &lvt, work, &lwk, &ret);
#else
    sgesvd_(&ju, &jvt, &M, &N, a, &la, w, u, &lu, vt, &lvt, work, &lwk, &ret);
#endif
    //    cout << "done svd" << endl;

    VISMatrix um(M, N);
    um = 0.0f;

    um.becomeLAPK(u, M, N);  //u from lapack
    U = um;  //U is to be returned
    //printLAPK(u, M, N);
    VISMatrix vm(N, N);
    vm = 0.0f;

    vm.becomeLAPK(vt, N, N);  //NOTE: vt here!

    V = vm.t();  //V is to be returned

    //printLAPK(vt, N, N);
    VISMatrix wm(N, N);
    int j;
    wm = 0.0f;
    for (j=0; j<VISmin(N, M); j++)
      wm.poke(j, j) = w[j];  //w from lapack
    for (; j<_cols; j++)
      wm.poke(j, j) = 0.0f;

    W = wm;  //W is to be returned
    //printLAPK(w, M, N);

    destroyLAPK(u);  destroyLAPK(vt);
    destroyLAPK(w);  destroyLAPK(work);
    destroyLAPK(a);

    if (ret == 0){
//	cout << "SVD called successfully." << endl << endl;
    }
    else{
      if(ret < 0)
        cout << "Illegal value for " << -(ret) << "th argument!" << endl;
      else
	{
        cout << "Error: " << ret << " bidiagonals not zero!" << endl;
	cout << "input is " << *this << endl;
	}
      //      exit(0);
      U = V = W = VISMatrix();
    }
}//svd

//Perform the QR decomposition of a matrix -- 9/02/99
void VISMatrix::QR(VISMatrix &Q, VISMatrix &R) const{
    //A=Q*R, A:MXN Q:MXmin(M,N) R: min(M,N)XN
    //The following stuff are needed by QRD of LAPACK//
    integer  M = _rows;
    integer  N = _cols;
    integer  lda = M;   //leading dimension of A: >=max(1,M)
    integer  lwork;     //dimension of work space: >=max(1,N).
    integer  ret;  //return value: 0=successful exit
    datatype  *a;     //input/output array for LAPACK
    //On entry, the M-by-N matrix A.
    //On exit, the elements on and above the diagonal of the array
    //contain the min(M,N)-by-N upper trapezoidal matrix R ( R is
    //upper triangular if m>= n); the elements below the diagonal, 
    //with the array TAU, represent the orthogonal matrix Q as a
    //product of min(m,n) elementary reflectors. 
    datatype  *tau;   //scalar factors of the elementary reflectors
    datatype  *work;  //workspace array, dimension lwork
    integer  minmn = MIN(M, N);  //dimension of array tau
    tau = new datatype [minmn];
    lwork = MAX(M, N) + 1;
    work = new datatype [lwork];

    a = createLAPK();
#ifdef VISMATRIXDOUBLE
    dgeqrf_(&M, &N, a, &lda, tau, work, &lwork, &ret); //QR from LAPACK
#else
    sgeqrf_(&M, &N, a, &lda, tau, work, &lwork, &ret); //QR from LAPACK
#endif

    datatype *r;
    char  jr = 'U';  //'U': Upper triangular or trapezoid
    VISMatrix rm(minmn, N);
    rm = 0.0f;
    r = rm.createLAPK();

#ifdef VISMATRIXDOUBLE
    dlacpy_(&jr, &M, &N, a, &lda, r, &minmn); //LAPACK copy
#else
    slacpy_(&jr, &M, &N, a, &lda, r, &minmn); //LAPACK copy
#endif
    rm.becomeLAPK(r, minmn, N);
    R = rm;  //R is to be returned 

    datatype *q;
    char  jq = 'L';  //'L': Lower triangular or trapezoid
    VISMatrix qm(M, minmn);
    qm = 0.0f;
    q = qm.createLAPK();
#ifdef VISMATRIXDOUBLE
    dlacpy_(&jq, &M, &minmn, a, &lda, q, &lda); //LAPACK copy
#else
    slacpy_(&jq, &M, &minmn, a, &lda, q, &lda); //LAPACK copy
#endif
#ifdef VISMATRIXDOUBLE
    dorgqr_(&M, &minmn, &minmn, q, &lda, tau, work, &lwork, &ret); //Q maker
#else
    sorgqr_(&M, &minmn, &minmn, q, &lda, tau, work, &lwork, &ret); //Q maker
#endif
    qm.becomeLAPK(q, M, minmn);
    Q = qm;  //Q is to be returned

    delete[] tau;  delete[] work;
    destroyLAPK(q), destroyLAPK(r);
    if (ret == 0){
//        cout << "QRD called successfully." << endl << endl;
    }
    else{
	cerr << "Error: QRD failed." << endl << endl;
        exit(0);
    }
}//QR

VISMatrix VISMatrix::inverseSVD() const {
    //A=U*W*transpose(V) ==>> A_inverse=V*W_inverse*U_transpose
    //The following stuff are needed by SVD of LAPACK//
    integer  M = _rows;
    integer  N = _cols;
    integer S = VISmax(M, N);
    char  ju  = 'A';   //'All' columns of U are returned
    char  jvt = 'A';   //'All' columns of VT are returned
    integer  la = M;   //leading dimension of A: >=max(1,M)
    integer  lu = M;   //leading dimension of U: >= M
    integer  lvt = N;  //leading dimension of VT: >= N
    integer  lwk;      //dimension of the array work
    integer  ret;      //return value: 0=successful exit
    datatype *a, *work;     //input array and work space
    datatype *u, *w, *vt;   //output arrays
    lwk = M*N + N*N + N;
    work = new datatype [lwk];   
//      u = new datatype [M*N];
//      w = new datatype [N];
//      vt = new datatype [N*N];
// the sizes above are supposed to be correct, but lapack appears to be broken
      u = new datatype [S*S];
      w = new datatype [S];
      vt = new datatype [S*S];

 
    a = createLAPK();
    //printLAPK(a, M, N);
#ifdef VISMATRIXDOUBLE
    dgesvd_(&ju, &jvt, &M, &N, a, &la, w, u, &lu, vt, &lvt, work, &lwk, &ret);
#else
    sgesvd_(&ju, &jvt, &M, &N, a, &la, w, u, &lu, vt, &lvt, work, &lwk, &ret);
#endif

    VISMatrix um(M, N);
    um = 0.0f;
    um.becomeLAPK(u, M, N);  //u from lapack
    um = ~um;  //U_transpose
    //printLAPK(u, M, N);
  
    VISMatrix vm(N, N);
    vm = 0.0f;
    vm.becomeLAPK(vt, N, N);  //NOTE: vt here!
    vm = ~vm;  //V
    //printLAPK(vt, N, N);
 
    VISMatrix wm(N, N);
    int rank = VISmin(_cols, _rows);
    wm = 0.0f;
    for (int j=0; j<rank; j++){ //W_inverse
      if( fabs(w[j]) > __MAT_EPSILON)  //depends on precison
	wm.poke(j, j) = 1.0f/w[j];  //w from lapack
    }
    //printLAPK(w, M, N);

    destroyLAPK(u);  destroyLAPK(vt);
    destroyLAPK(w);  destroyLAPK(work);
    destroyLAPK(a);

    VISMatrix inverse(N, M);
    inverse = 0.0f;
    inverse = vm*wm*um; //A_inverse=V*W_inverse*U_transpose
    return inverse;
}//inverseSVD

VISMatrix VISMatrix::inverseGJ() const{
    //A*X=B 
    //The following stuff are needed by LAPACK//
    integer  NA = _rows;
    integer  NB = NA;
    integer  lda = NA;  //leading dimension of A: >=max(1,M)
    integer  ldb = NA;  //leading dimension of b: >=max(1,M)
    integer  ret;      //return value: 0=successful exit
    integer  *ipiv;    //pivot indices array of dimension M

    if(_cols!=_rows){
	cout << "inverseGJ(): matrix must be square" << endl;
        exit(0);
    }
    datatype *a, *b;
    a = createLAPK();
    VISMatrix B(NA, NB);
    B = 0.0f;
    for(int cnt=0;cnt<NB;cnt++)
	B.poke(cnt, cnt) = 1.0f;
    b = B.createLAPK();
    ipiv = new integer [NA];

#ifdef VISMATRIXDOUBLE
    dgesv_(&NA, &NB, a, &lda, ipiv, b, &ldb, &ret);
#else
    sgesv_(&NA, &NB, a, &lda, ipiv, b, &ldb, &ret);
#endif
    VISMatrix a_inv;
    a_inv.becomeLAPK(b, NA, NA); //depends on a or b
    delete[] ipiv; destroyLAPK(a); destroyLAPK(b);
    return a_inv;
}//inverseGJ


VISVector solvLinearEqn(const VISMatrix& A,const VISVector& b){
    // solve the linear system of eqns:
    // Ax=b, simply return x=inv(A)b
    return ((A.inverseSVD())*b);
}

#else

datatype VISMatrix::det() const{
    cout << "datatype VISMatrix::det() const" << endl;
    cout << "LAPACK is needed to execute this function." << endl;
    return(0.0f);
}
void VISMatrix::svd(VISMatrix &U, VISMatrix &W, VISMatrix &V) const{
    cout << "void VISMatrix::svd(VISMatrix &U, VISMatrix &W, VISMatrix &V) const" << endl;
    cout << "LAPACK is needed to execute this function." << endl;
    cout << "done svd stub not quite" << endl;
    U=V=W=VISMatrix();
    cout << "done svd stub" << endl;
}
VISMatrix VISMatrix::inverseSVD() const {
    cout << "VISMatrix VISMatrix::inverseSVD() const" << endl;
    cout << "LAPACK is needed to execute this function." << endl;
    return VISMatrix();
}
VISMatrix VISMatrix::inverseGJ() const{
    cout << "VISMatrix VISMatrix::inverseGJ() const" << endl;
    cout << "LAPACK is needed to execute this function." << endl;
    return VISMatrix();
}

VISVector solvLinearEqn(const VISMatrix& A,const VISVector& b){
    cout << "VISVector solvLinearEqn(const VISMatrix& A,const VISVector& b)" << endl;
    cout << "LAPACK is needed to execute this function." << endl;
    return VISVector();
}
void VISMatrix::QR(VISMatrix &Q, VISMatrix &R) const{
    cout << "void VISMatrix::QR(VISMatrix &Q, VISMatrix &R) const" << endl;
    cout << "LAPACK is needed to execute this function." << endl;
}
#endif



VISMatrix VISIdentity(int size){
    VISMatrix ans(size,size);
    ans=0.0f;
    for (int cnt=0;cnt<size;cnt++)
	ans.poke(cnt,cnt)=1.0f;

    return(ans);
}

VISMatrix homogeneousTransform(VISMatrix &mother, VISVector &vother)
{
    if((mother.c()==3) && (mother.r()==3) && (vother.n() == 3) )
    {
        VISMatrix tmp(1,4);
        tmp = 0.0;
        tmp.poke(0,3) = 1.0;

        return((mother.concat(vother)).concatRow(tmp));     
    }
    else
    {
        cout << "VISMatrix(const VISMatrix& mother, const VISVector& vother)" << endl;
        cout << "You have attempted to create an HVISMatrix from "
             << "a rotation matrix that is not 3x3 and/or " << endl;
        cout << "a translation vector that is not 3x1." << endl;
        exit(0);
		return 0; // keeps compiler happy
    }
}//homogeneousTransform


//*******CONCATINATION******* the below looks ugly but works.
//note you can not concatRow a matrix and a vector because a
//vector is ALWAYS one column.  if you want to stick a vector
//onto a matrix as a row use concatRow(v,m) where v is vector
//and m is matrix
VISMatrix VISMatrix::concat(const VISMatrix& second)
{
  //concatinate cols i.e. a=[1 2 3  b=[ 7  8  9
  //                         4 5 6]    10 11 12]
  // a = a.concatcols(b)   c=[1 2 3  7  8  9
  //                          4 5 6 10 11 12]
  if (((c()==0)&&(r()==0)) &&
	  ((second.c()==0)&&(second.r()==0)))
	{
	  cout << "VISMatrix concat(const VISMatrix& "
		   << "first,const VISMatrix& second)" << endl;
	  cout << "You have attempted to row concatinate matricies both"
		   << "of which have zero size." << endl;
    }
  else if ((c()==0)&&(r()==0))
	{
	  VISMatrix ans;
	  ans=second;
	  return ans;
    }
  else if ((second.c()==0)&&(second.r()==0))
	{
	  VISMatrix ans;
	  ans=*this;
	  return ans;
	}	
  else if (r()!=second.r())
	{
	  cout << "VISMatrix concat(const VISMatrix& "
		   << "first,const VISMatrix& second)" << endl;
	  cout << "You have attempted to column concatinate matricies "
		   << "of different numbers of rows." << endl;
	  cout << "r()  " << r() << "     "
		   << "second.r()  " << second.r() << endl;
	}
  else
	{
	  int numrows = r();
	  int numcols = c()+second.c();
	  int firstcols = c();
	  VISMatrix ans(numrows,numcols);
	  for (int row=0;row<numrows;row++)
	    for (int col=0;col<numcols;col++)
		  if(col<firstcols)
		    ans.poke(row,col) = peek(row,col);
		  else
		    ans.poke(row,col) = second(row,col-firstcols);
	  return ans;
    }
  return 0; // Error condition
}


VISMatrix VISMatrix::concat(const VISVector& second)
{
  //concatinate cols i.e. a=[1 2 3  b=[ 7 
  //                         4 5 6]    10]
  // c=concatcols(a,b)   c=[1 2 3  7 
  //                        4 5 6 10]
  if ((second.n()==0) && ((c()==0)&&(r()==0)))
	{
	  cout << "VISMatrix concat(const VISMatrix& "
		   << "first,const VISVector& second)" << endl;
	  cout << "You have attempted to row concatinate a vector and "
		   << "a matrix, both of which have zero size." << endl;
    }
  else if ((c()==0)&&(r()==0))
	{
	  VISMatrix ans;
	  ans=second.mat();
	  return ans;
	}
  else if (second.n()==0)
	{
	  VISMatrix ans;
	  ans=*this;
	  return ans;
	}	
  else if (r()!=second.n())
	{
	  cout << "VISMatrix concat(const VISMatrix& "
		   << "first,const VISVector& second)" << endl;
	  cout << "You have attempted to column concatinate a matrix "
		   << "and a vector that have different numbers of rows/elements."
		   << endl;
	  cout << "r()  " << r() << "     "
		   << "second.n()  " << second.n() << endl;
	  exit(0);
    }
  else
	{
	  int numrows = second.n();
	  int numcols = c()+1;
	  int firstcols = c();
	  VISMatrix ans(numrows,numcols);
	  for (int row=0;row<numrows;row++)
	    for (int col=0;col<numcols;col++)
		  if(col<firstcols)
		    ans.poke(row,col) = peek(row,col);
		  else
		    ans.poke(row,col) = second[row];
	  return ans;
    }
  return 0; // Error condition
}

VISMatrix VISMatrix::concatRow(const VISMatrix& second){
  //concatinate rows i.e. a=[1 2 3  b=[ 7  8  9
  //                         4 5 6]    10 11 12]
  // c=concatRows(a,b)   c=[1  2  3
  //                        4  5  6
  //                        7  8  9
  //                       10 11 12]
  if (((c()==0)&&(r()==0)) &&
	  ((second.c()==0)&&(second.r()==0))){
	cout << "VISMatrix concatRows(const VISMatrix& "
	     << "first,const VISMatrix& second)" << endl;
	cout << "You have attempted to row concatinate matricies of zero size "
	     << endl;
	exit(0);
  }
  else if ((c()==0)&&(r()==0)){
	VISMatrix ans;
	ans=second;
	return ans;
  }
  else if ((second.c()==0)&&(second.r()==0)){
	VISMatrix ans;
	ans=*this;
	return ans;
  }	
  else if (c()!=second.c()){
	cout << "VISMatrix concatRows(const VISMatrix& "
	     << "first,const VISMatrix& second)" << endl;
	cout << "You have attempted to row concatinate matricies of different "
	     << "numbers of columns." << endl;
	cout << "c()  " << c() << "     "
	     << "second.c()  " << second.c() << endl;
	exit(0);
  }
  else{
	int numrows = r()+second.r();
	int numcols = c();
	int firstrows = r();
	VISMatrix ans(numrows,numcols);
	for (int row=0;row<numrows;row++)
	  for (int col=0;col<numcols;col++)
		if(row<firstrows)
		  ans.poke(row,col) = peek(row,col);
		else
		  ans.poke(row,col) = second[row-firstrows][col];
	return ans;
  }
  return 0; // error condition
}

