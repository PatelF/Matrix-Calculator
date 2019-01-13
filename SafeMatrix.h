////////////////////////////////////////////////////////////////////////////////
//
// SafeMatrix.h
// Defines class SafeMatrix
// The SafeMatrix class should provide safe access to 2D arrays
// as well as the elemntary Matrix operations +, * and transpose
//

#ifndef SAFEMATRIX_H
#define SAFEMATRIX_H

// Invariant defintions

#define INVARIANT_FAILURE(DATA) do {					\
    std::cerr << "Method " << DATA << ": Invariant Failure; Aborting" << std::endl; \
    abort();								\
  } while (0)

#define INVARIANT_TEST(TEST, DATA) do {		\
    if ((TEST)) {				\
      INVARIANT_FAILURE(DATA);			\
    }						\
  } while(0)

// Use _numRows for these two

#define NOT_A_MATRIX  -1
#define MATRIX_DELETED -2

#define isNaN(X) (X != X)

// Reset NaN is neccessary in case someone changed it

#define RESET_NAN do {							\
    if (!isNaN(_nan)) {							\
      errorMsg("Warning: _nan was set to something other than NaN; resetting");	\
      _nan = 0.0/0.0;							\
    }									\
  } while (0)

// Given SafeMatrix M, the access to _valid_ locations M(X,Y) is given by this macro:

#define MATRIX(M,X,Y) (M)._data[X*(M)._numCols+Y]

////////////////////////////////
//
// Dimensions
//

struct Dimensions {
  int rows;
  int cols;
};

/////////////////////////////////
//
// SafeMatrix Class
//

class SafeMatrix {
public:
  // Methods
  void       transpose();
  bool       invert();                       // Invert Matrix
  SafeMatrix solve(const SafeMatrix& A);     // Solve MX = A
  
  bool       appendRow(const int cols, const float data[]);
  bool       appendColumn(const int rows, const float data[]);
  bool       deleteRow(const int rowNumber);
  bool       deleteColumn(const int rowNumber);
  bool       swapRow(const int row1, const int row2);
  bool       swapColumn(const int column1, const int column2);
  bool       append(const SafeMatrix& m);

  bool       isNaM() const;                  // Is Not-a-Matrix

  SafeMatrix row(const int rowNumber) const;
  SafeMatrix column(const int columnNumber) const;
  
  Dimensions dimensions() const;

  static SafeMatrix identity(const int n);
  
  // Operators
  float&     operator()(int row, int col);
  SafeMatrix operator+ (const SafeMatrix& m) const;
  SafeMatrix operator- (const SafeMatrix& m) const;
  SafeMatrix operator* (const SafeMatrix& m) const;
  void       operator= (const SafeMatrix& m);
  bool       operator==(const SafeMatrix& m) const;
  bool       operator!=(const SafeMatrix& m) const;

  // Constructors/Destructor
  SafeMatrix();                                                         // Size zero matrix
  SafeMatrix(const SafeMatrix& m);
  SafeMatrix(const Dimensions& d);                                      // Uninitialized matrix
  SafeMatrix(const int rows, const int cols);                           // Uninitialized matrix
  SafeMatrix(const Dimensions& d, const float initVal);                 // Initialized matrix
  SafeMatrix(const int rows, const int cols, const float initVal);      // Initialized matrix
  SafeMatrix(const int rows, const float data[]);                       // Initialized vector (R,1)
 ~SafeMatrix();

//#include "myPublicAdditions.h"
  
private:
  friend std::ostream& operator<<(std::ostream& os, const SafeMatrix& m);
  void errorMsg(const char msg[]) const;             // Display error message to std::cerr

  float _nan = 0.0/0.0;                              // Return this for out-of-bounds access
                                                     // Make sure to reset it if it changes
  int     _numRows;
  int     _numCols;
  int     _dataSpaceAllocated;
  float*  _data;

//#include "myPrivateAdditions.h"

};

#endif
