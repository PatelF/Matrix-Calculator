#include <iostream>
#include "SafeMatrix.h"
#include <stdlib.h>

using namespace std;

//------------------Solve/Inverse---------------//
bool SafeMatrix::invert(){

  if(_numRows != _numCols){
    return false;
  }

  if(_numRows < 0 || _numCols < 0){
    return false;
  }

  if(_data == NULL){
    return false;
  }

  int i,j,k,n;
  float a[100][200],t;

  n = _numRows;

  k = 0;

  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      a[i][j] = _data[k];
      k++;
    }
  }

  k = 0;

  for(i = 0; i < n; i++){
    for(j = n; j < 2*n; j++){
      if(i == j - n)
        a[i][j] = 1;
      else
        a[i][j] = 0;
    }
  }

  for(i = 0; i < n; i++){
    t = a[i][i];
    if(t == 0)
      return false;

    for( j = i; j < 2*n; j++){
      if(t == 0)
        return false;

      a[i][j] = a[i][j] /t;
    }

    for(j = 0; j < n; j++){
      if(i != j){
        t = a[j][i];
        for(k = 0; k < 2*n; k++){
          a[j][k] = a[j][k] - t*a[i][k];
        }
      }
    }
  }

  delete[] _data;
  _data = new (std::nothrow) float [_dataSpaceAllocated];
  int count = 0;

  for(i = 0; i < n; i++){

    for(j = n; j < 2*n; j++){
      _data[count] = a[i][j];
      count++;
    }
  }

  return true;
}

SafeMatrix SafeMatrix::solve(const SafeMatrix& A){
  INVARIANT_TEST(_numRows < -1, "solve");

  if(_numRows != _numCols){
    return *this;
  }

  SafeMatrix s(A._numRows, 1);

  SafeMatrix invertible(*this);

  invertible.invert();

  s = invertible*A;

  return s;
  
}

//---------------------Bools---------------------//
bool SafeMatrix::append(const SafeMatrix& m){

  if(m._numCols < 0 || m._numRows < 0 || _numCols < 0 || _numRows < 0 || _dataSpaceAllocated < 0 || m._dataSpaceAllocated < 0){
    return false;
  }

  if(_numCols == m._numCols){

    int d = _numRows + m._numRows;

    int i = 0;
    float data[_numCols];

    while(i < m._dataSpaceAllocated){

      for(int j = 0; j < _numCols; j++){

        data[j] = m._data[i];
        i++;
      }

      SafeMatrix::appendRow(m._numCols, data);
    }

    if(_numRows != d)
      _numRows = d;
    return true;
  }

  else if (_numRows == m._numRows) {

    int d = _numCols + m._numCols;
    int j = 0;
    int k = 0;

    float data[_numRows];

    while(j < m._dataSpaceAllocated){

      for(int i = 0; i < _numRows;i++){

        int start = i *m._numCols + k;
        data[i] = m._data[start];
        j++;
      }
      k++;
      SafeMatrix::appendColumn(m._numRows, data);
    }

    if(_numCols != d)
      _numRows = d;

    return true;
  }

  return false;
}


bool SafeMatrix::isNaM()const{

  if(_numRows < 0 || _numCols < 0 || _dataSpaceAllocated < 0){
    return true;
  }
  return false;
}

bool SafeMatrix::swapRow(const int row1, const int row2){

  if(row1 > _numRows || row2 > _numRows){
    return false;
  }

  float firstR[_numCols];
  float secondR[_numCols];

  int j = 0;
  int k = 0;

  for(int i = 0; i < _dataSpaceAllocated; i++){
    if(i >= (row1*_numCols) && i < ((row1*_numCols) + _numCols)){
      firstR[j] = _data[i];
      j++;
    }
    else if(i >= (row2*_numCols) && i < ((row2*_numCols) + _numCols)){
      secondR[k] = _data[i];
      k++;
    }
  }

  j = 0;
  k = 0;
  int rowNum = 0;

  for(int i = 0; i < _dataSpaceAllocated; i++){
    if(i != 0 && i % _numCols == 0){
      rowNum++;
    }

    if (rowNum == row1){
      _data[i] = secondR[j];
      j++;
    }
    else if (rowNum == row2){
      _data[i] = firstR[k];
      k++;
    }
  }
  return true;
}

bool SafeMatrix::swapColumn(const int column1, const int column2){
  if(column1 > _numCols || column2 > _numCols){
    return false;
  }

  (*this).transpose();
  (*this).swapRow(column1,column2);
  (*this).transpose();

  return true;
}

bool SafeMatrix::deleteRow(const int rowNumber){

  if(rowNumber > _numRows - 1){
    return false;
  }
  _dataSpaceAllocated = (_numRows-1)*(_numCols);

  float* tmp = new (std::nothrow) float[_dataSpaceAllocated];

  int front = (rowNumber)* _numCols;
  int end = (rowNumber*+_numCols) + (_numCols - 1);

  int j = 0;
  int k = 0;
  for(int i = 0; i < _dataSpaceAllocated + _numCols; i++){
    if(i >= front && i <= end){
      k += 1;
    }
    else{
      tmp[j] = _data[i];
      j++;
    }

  }

  delete[] _data;

  _numRows--;
  _data = new (std::nothrow) float[_dataSpaceAllocated];

  _data = tmp;
  return true;


}

bool SafeMatrix::deleteColumn(const int columnNumber){

  if(columnNumber > _numCols - 1){
    return false;
  }

  (*this).transpose();
  (*this).deleteRow(columnNumber);
  (*this).transpose();

  return true;
}

bool SafeMatrix::appendColumn(const int rows, const float data[]) {
  if(rows != _numRows){
    return false;
  }

  (*this).transpose();
  (*this).appendRow(rows,data);
  (*this).transpose();

  return true;
}

bool SafeMatrix::appendRow(const int cols, const float data[]){
  if(cols != _numCols){
    return false;
  }
  _dataSpaceAllocated += _numCols;
  _numRows++;

  float* tmp = new (std::nothrow) float[_dataSpaceAllocated];

  for(int i = 0; i < _dataSpaceAllocated - _numCols; i++){
    tmp[i] = _data[i];
  }

  for(int i = _dataSpaceAllocated - _numCols; i < _dataSpaceAllocated; i++){
    tmp[i] = data[i - _dataSpaceAllocated+_numCols];
  }

  delete[] _data;
  _data = tmp;

  return true;
}

void SafeMatrix::transpose() {

  INVARIANT_TEST(_numRows < -1, "row");

  SafeMatrix s(_numCols, _numRows, 0);

  int pos = 0;

  for(int i = 0; i < s._numRows; i++){
    for(int j = 0; j < s._numCols; j++){
      s._data[pos] = (*this)(j,i);
      pos++;
    }
  }

  for(int i = 0; i < _dataSpaceAllocated; i++){
    _data[i] = s._data[i];
  }

  _numRows = s._numRows;
  _numCols = s._numCols;
}

//-----------------Find Row/Col -----------------------//

SafeMatrix SafeMatrix::row(const int rowNumber)const{

  INVARIANT_TEST(_numRows < -1, "row");

  if(rowNumber < 0 || rowNumber > _numRows - 1){
    SafeMatrix bad (-1,-1);

    return bad;
  }

  SafeMatrix s (1, _numCols);

  int start = (rowNumber*_numCols);
  int j = 0;
  for(int i = start; i < start + _numCols; i++){
    s._data[j] = _data[i];
    j++;
  }

  return s;
}

SafeMatrix SafeMatrix::column(const int columnNumber)const{

  INVARIANT_TEST(_numRows < - 1, "column");

  if(columnNumber < 0 || columnNumber > _numCols - 1){
    SafeMatrix bad(-1,-1);
    return bad;
  }

  SafeMatrix s (1,_numRows);

  SafeMatrix trans(*this);
  trans.transpose();

  s = trans.row(columnNumber);

  s.transpose();

  return s;
}

//-----------------OPERATORS---------------------------//
float& SafeMatrix::operator()(int row, int col) {

  INVARIANT_TEST(_numRows < -1, "operator()");

  if(row < 0 || row > _numRows-1 || col < 0 || col > _numCols-1){
    return _nan;
  }

  int location = (_numCols)*(row) + (col);
  return _data[location];
}

SafeMatrix SafeMatrix::operator+(const SafeMatrix& m) const {

  INVARIANT_TEST(_numRows < -1, "operator+");
  INVARIANT_TEST(m._numRows < -1, "operator+(m)");

  if(_dataSpaceAllocated == 0 || m._dataSpaceAllocated == 0){
    SafeMatrix bad(-1,-1);
    return bad;
  }

  if(_numRows != m._numRows || _numCols != m._numCols || _dataSpaceAllocated != m._dataSpaceAllocated){
    SafeMatrix bad(-1,-1);
    return bad;
  }

  if(_data == 0 || m._data == 0){
    SafeMatrix bad(-1,-1);
    return bad;
  }

  SafeMatrix s(_numRows, _numCols);
  if(s._dataSpaceAllocated == 0){
    SafeMatrix bad(-1,-1);
    return bad;
  }
  
  for(int i = 0; i < s._dataSpaceAllocated; i++){
    s._data[i] = _data[i] + m._data[i];
  }

  return s;
}

SafeMatrix SafeMatrix::operator-(const SafeMatrix& m)const{
  INVARIANT_TEST(_numRows < -1, "operator-");
  INVARIANT_TEST(m._numRows < -1, "operator-(m)");

  if(_dataSpaceAllocated == 0 || m._dataSpaceAllocated == 0){
    SafeMatrix bad(-1,-1);
    return bad;
  }

  if(_numRows != m._numRows || _numCols != m._numCols || _dataSpaceAllocated != m._dataSpaceAllocated){
    SafeMatrix badMatrix(-1,0);
    return badMatrix;
  }

  SafeMatrix s(_numRows, _numCols);
  if(s._dataSpaceAllocated == 0){
    SafeMatrix bad(-1,-1);
    return bad;
  }

  for(int i = 0; i < s._dataSpaceAllocated; i++){
    s._data[i] = _data[i] - m._data[i];
  }

  return s;
}

SafeMatrix SafeMatrix::operator*(const SafeMatrix& m) const {
  INVARIANT_TEST(_numRows < -1, "operator*");
  INVARIANT_TEST(m._numRows < -1, "operator*(m)");

  if(m._numRows < 0 || m._numCols < 0){
    SafeMatrix badMatrix(-1,-1);
    return badMatrix;
  }

  if(_numCols != m._numRows){
    SafeMatrix badMatrix(-1,-1);
    return badMatrix;
  }

  SafeMatrix s(_numRows, m._numCols,0);
  if(s._dataSpaceAllocated == 0){
    SafeMatrix bad(-1,-1);
    return bad;
  }

  for(int k = 0; k < _numRows; k++){

    for(int j = 0; j < m._numCols; j++){

      for(int i = 0; i < _numCols; i++){

        s._data[k*m._numCols+j] += _data[k*_numCols + i]*m._data[i*m._numCols+j];
      }
    }
  }
  return s;
}

void SafeMatrix::operator=(const SafeMatrix& m) {
  INVARIANT_TEST(_numRows < -1, "operator=");
  INVARIANT_TEST(m._numRows < -1, "operator=(m)");
  
  _numRows = m._numRows;
  _numCols = m._numCols;

  float* tmp = new (std::nothrow) float [m._dataSpaceAllocated];

  for(int i = 0; i < m._dataSpaceAllocated; i++){
    tmp[i] = m._data[i];
  }

  delete[] _data;
  _data = tmp;

  return;
}

bool SafeMatrix::operator==(const SafeMatrix& m)const{
  INVARIANT_TEST(_numRows < -1, "operator=");
  INVARIANT_TEST(m._numRows < -1, "operator=(m)");

  if(_numRows != m._numRows || _numCols != m._numCols|| _dataSpaceAllocated != m._dataSpaceAllocated||m._numRows < 0 || m._numCols < 0){
    return false;
  }

  for(int i = 0; i < _dataSpaceAllocated; i++){
    if (_data[i] != m._data[i]){
      return false;
    }
  }
  return true;
}

bool SafeMatrix::operator!=(const SafeMatrix& m)const{
  INVARIANT_TEST(_numRows < -1, "operator=");
  INVARIANT_TEST(m._numRows < -1, "operator=(m)");

  if(_numRows != m._numRows || _numCols != m._numCols|| _dataSpaceAllocated != m._dataSpaceAllocated){
    return true;
  }

  for(int i = 0; i < _dataSpaceAllocated; i++){
    if (_data[i] != m._data[i]){
      return true;
    }
  }
  return false; 
}

//-----------------DIMENSIONS-----------------//
Dimensions SafeMatrix::dimensions() const {

  Dimensions dim;
  dim.rows = _numRows;
  dim.cols = _numCols;

  return dim;
}
//------------------------------------------//
SafeMatrix SafeMatrix::identity(const int n){

  int tmp[n][n];

  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i == j){
        tmp[i][j] = 1;
      }
      else{
        tmp[i][j] = 0;
      }
    }
  }

  SafeMatrix s(n,n);
  int k = 0;

  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      s._data[k] = tmp[i][j];
      k++;
    }
  }

  // SafeMatrix identity(n,n,0);

  // for (int i = 0; i < _numCols; i++) {
  //   identity(i,i) = 1;
  // }

  // return identity;
  return s;
}

//----------------CONSTRUCTORS--------------//
SafeMatrix::SafeMatrix(){
  _numRows = 0;
  _numCols = 0;
  _dataSpaceAllocated = 0;
  _data = new(std::nothrow) float[1];
}

SafeMatrix::SafeMatrix(const SafeMatrix& m){

  if(m._numRows < 0|| m._numCols < 0){
    _numRows = -1;
    _numCols = -1;
    _dataSpaceAllocated = 0;
    _data = NULL;
    return;
  }
  
  _numRows = m._numRows;
  _numCols = m._numCols;
  _dataSpaceAllocated = m._dataSpaceAllocated;

  _data = new (std::nothrow) float[_dataSpaceAllocated];

  if(_data == 0){
    return;
  }

  for(int i = 0; i < _dataSpaceAllocated; i++){
    _data[i] = m._data[i];
  } 
}

SafeMatrix::SafeMatrix(const Dimensions& d){

  if(d.rows < 0 || d.cols < 0){
    _numRows = -1;
    _numCols = -1;
    _dataSpaceAllocated = 0;
    _data = NULL;
    return;
  }

  _numRows = d.rows;
  _numCols = d.cols;

  _dataSpaceAllocated = d.rows*d.cols;


  _data = new (std::nothrow) float[_dataSpaceAllocated];

  if(_data == 0){
    return;
  }

}

SafeMatrix::SafeMatrix(const int rows, const int cols){

  if(rows < 0 || cols < 0){
    _numRows = -1;
    _numCols = -1;
    _dataSpaceAllocated = 0;
    _data = NULL;
    return;
  }

  _numRows = rows;
  _numCols = cols;

  _dataSpaceAllocated = rows*cols;


  _data = new (std::nothrow) float[_dataSpaceAllocated];

  if(_data == 0){
    return;
  }


}

SafeMatrix::SafeMatrix(const Dimensions& d, const float initVal){
  _numCols = d.cols;
  _numRows = d.rows;

  if(_numCols < 0 || _numRows < 0){
    _numRows = -1;
    _numCols = -1;
    _dataSpaceAllocated = 0;
    _data = NULL;
    return;
  }

  _dataSpaceAllocated = _numCols*_numRows;

  _data = new (std::nothrow) float[_dataSpaceAllocated];

  if(_data == 0){
    return;
  }

  for(int i = 0; i < _dataSpaceAllocated; i++){
    _data[i] = initVal;
  }
}

SafeMatrix::SafeMatrix(const int rows, const int cols, const float initVal){
  if(rows < 0 || cols < 0){
    _numRows = -1;
    _numCols = -1;
    _dataSpaceAllocated = 0;
    _data = NULL;
    return;
  }

  _numRows = rows;
  _numCols = cols;
  _dataSpaceAllocated = rows*cols;

  _data = new (std::nothrow) float [_dataSpaceAllocated];

  if(_data == 0){
    return;
  }

  for(int i = 0; i < _dataSpaceAllocated; i++){
    _data[i] = initVal;
  }
}

SafeMatrix::SafeMatrix(const int rows, const float data[]){

  if(rows < 0){
    _numRows = 1;
    _numCols = -1;
    _dataSpaceAllocated = 0;
    _data = NULL;
    return;
  }

  _numRows = rows;
  _numCols = 1;
  _dataSpaceAllocated = _numRows;

  _data = new (std::nothrow) float[_dataSpaceAllocated];

  if(_data == 0){
    return;
  }

  for(int i = 0; i < _dataSpaceAllocated; i++){
    _data[i] = data[i];
  }

}

SafeMatrix::~SafeMatrix(){
  _numRows = -2;
  delete[] _data;
}

std::ostream& operator<<(std::ostream& os, const SafeMatrix& m) {
  INVARIANT_TEST(m._numRows < NOT_A_MATRIX, "ostream::operator<<()");
  if (m._numRows < NOT_A_MATRIX) {
    m.errorMsg("Method ostream::operator<<: attempting to output deleted SafeMatrix");
    abort();
  }
  if (m._numRows == NOT_A_MATRIX) {
    os << "Not-a-Matrix";
    return os;
  }
  if (m._numRows == 0)
    os << "[]";
  for (int i = 0; i < m._numRows; ++i) {
    os << "[";
    for (int j = 0; j < m._numCols; ++j) {
      os << MATRIX(m,i,j);
      if (j < (m._numCols - 1))
        os << ", ";
    }
    os << "]";
    if (i < (m._numRows - 1))
      os << std::endl;
  }
  return os;
}

void SafeMatrix::errorMsg(const char msg[]) const {
  std::cerr << msg << std::endl;
}

int main(const int argc, const char* const argv[]){

  cout << "Choose Number of Rows" << endl;
  int rows;
  cin >>rows;

  cout << "\nChoose Number of Columns" << endl;
  int cols;
  cin >> cols;

  int dataSpace = cols*rows;

  float data[dataSpace];

  int k = 0;

  cout << "\nEnter the numbers in each row and press enter" << endl;

  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      int num;
      cin >> num;
      data[k] = num;
      k++;      
    }
  }

  SafeMatrix s(1, cols, 0);

  for(int i = 0; i < rows; i++){
    if(i == 0){
      s.appendRow(cols, data);
    }
    else{
      s.appendRow(cols, data + cols*i);
    }
  }
  s.deleteRow(0);

  cout << "\nThis is your Matrix\n" << s << endl;

  cout << "\nWhat would you like to do?" << endl;
  cout << "Choose an option and press enter" << endl;

  int input;
  cout << "1. Invert" << endl;
  cout << "2. Transpose" << endl;
  cout << "3. Solve" << endl;
  cout << "4. Add Another Matrix" << endl;
  cout << "5. Subtract Another Matrix" << endl;
  cout << "6. Multiply Another Matrix" << endl;

  cin >> input;

  if(input == 1){

    if(!s.invert())
      cout << "Unfortunately your matrix is not invertible" << endl;
    else{
      cout << "\nInverting Your Matrix" << endl;
      cout << s << endl;
    }
  }

  else if (input == 2){
    cout << "\nTranposed Matrix" << endl;
    s.transpose();
    cout << s << endl;
  }

  else if (input == 3){
    cout << "\nEnter the number in each row and press enter" << endl;

    for(int i = 0; i < rows; i++){
      int num;
      cin >> num;
      data[i] = num;     
    }

    SafeMatrix b(1, 1, 0);

    for(int i = 0; i < rows; i++){
      if(i == 0){
        b.appendRow(1, data);
      }
      else{
        b.appendRow(1, data + i);
      }
    }
    b.deleteRow(0);

    cout << "\nThe solution to your matrix is\n" << s.solve(b) << endl;
  }

  else if (input == 4){

    cout << "Enter another matrix to add" << endl;

    float data[dataSpace];

    int k = 0;

    cout << "\nEnter the numbers in each row and press enter" << endl;

    for(int i = 0; i < rows; i++){
      for(int j = 0; j < cols; j++){
        int num;
        cin >> num;
        data[k] = num;
        k++;      
      }
    }

    SafeMatrix s2(1, cols, 0);

    for(int i = 0; i < rows; i++){
      if(i == 0){
        s2.appendRow(cols, data);
      }
      else{
        s2.appendRow(cols, data + cols*i);
      }
    }
    s2.deleteRow(0);

    cout <<"The added matrix is\n" << s.operator+(s2) << endl;
  }

  else if (input == 5){

    cout << "Enter another matrix to subtract" << endl;

    float data[dataSpace];

    int k = 0;

    cout << "\nEnter the numbers in each row and press enter" << endl;

    for(int i = 0; i < rows; i++){
      for(int j = 0; j < cols; j++){
        int num;
        cin >> num;
        data[k] = num;
        k++;      
      }
    }

    SafeMatrix s2(1, cols, 0);

    for(int i = 0; i < rows; i++){
      if(i == 0){
        s2.appendRow(cols, data);
      }
      else{
        s2.appendRow(cols, data + cols*i);
      }
    }
    s2.deleteRow(0);

    cout <<"The subtracted matrix is\n" << s.operator-(s2) << endl;
  }

  else if (input == 6){

    cout << "Enter another matrix to multiply" << endl;

    float data[dataSpace];

    int k = 0;

    cout << "\nEnter the numbers in each row and press enter" << endl;

    for(int i = 0; i < rows; i++){
      for(int j = 0; j < cols; j++){
        int num;
        cin >> num;
        data[k] = num;
        k++;      
      }
    }

    SafeMatrix s2(1, cols, 0);

    for(int i = 0; i < rows; i++){
      if(i == 0){
        s2.appendRow(cols, data);
      }
      else{
        s2.appendRow(cols, data + cols*i);
      }
    }
    s2.deleteRow(0);

    cout <<"The multiplied matrix is\n" << s.operator*(s2) << endl;
  }

  return 0;
}
// //--------------MAIN-----------//
// #ifndef MARMOSET_TESTING
// int main() {
//   SafeMatrix m(3,3,1);
//   cout << "Matrix m:\n" << m << endl;

//   float data[] = {1,0,0,0,1,0,0,0,1};
//   SafeMatrix data2(1,3,1);

//   data2.appendRow(3,data);
//   data2.appendRow(3, data+3);
//   data2.appendRow(3,data+6);
//   data2.deleteRow(0);

//   SafeMatrix s1(2,3,1);
//   SafeMatrix s2(2,5,2);

//   s1.append(s2);

//   s1.invert();

//   cout << "Appending s1 and s2\n" << s1 << endl;

//   cout << "After adding rows\n" << data2 << endl;

//   data2.invert();
//   cout << "Invert \n" << data2 << endl;
//   float data3[] = {2,3};

//   SafeMatrix data4(1,1,1);
//   data4.appendRow(1,data3);
//   data4.appendRow(1,data3 + 1);

//   cout << data4 << endl;
//   data2.solve(data4);

// }
// #endif