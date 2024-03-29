#include <cstddef>
#include <iostream>
#include <math.h>
#include <algorithm>


class Matris {
public:
   // constructors
   Matris(){
      m_rows =  0;
      m_cols =  0;
      m_capacity = 0;
      m_vec = new double[0];
    }
    Matris(int x, int y){
      m_rows = y;
      m_cols = x;
      m_capacity =  x * y;
      m_vec = new double[x*y]();
    }
    Matris(Matris const & m){
      m_rows = m.m_rows;
      m_cols = m.m_cols;
      m_capacity = m.m_capacity;
      m_vec = new double[m.rows() * m.cols()]();
      for (int i = 0; i < rows() * cols(); i++) {
        m_vec[i] = m.m_vec[i];
      }
    }

    Matris(int cols, int rows, std::initializer_list<double> l){
      m_rows = rows;
      m_cols = cols;
      m_capacity = rows * cols;
      m_vec = new double[l.size()]();
      std::copy(l.begin(), l.end(), m_vec);


      int valid = (*this).validate();
      if (valid != -1) {
        std::cerr << "Row error at row: " <<  valid << '\n';
      }
    }


    ~Matris(){
    // std::cerr << "Before Free" << '\n';
     delete [] m_vec;
     // std::cerr << "after Free" << '\n';

   }

   // operators
   Matris& operator= (Matris const& m){
      m_rows = m.m_rows;
      m_cols = m.m_cols;
      m_capacity = m.m_capacity;
      double * temp_vec = new double[m.cols()*m.rows()]();

      for (int i = 0; i < rows() * cols(); i++) {
        temp_vec[i] = m.m_vec[i];
      }
      // std::cerr << "Before free 1" << '\n';
      delete [] m_vec;
      // std::cerr << "After free 1" << '\n';
      // double * m_vec;
      m_vec = temp_vec;

      return *this;
    }

   Matris& operator*=(Matris m){
      if (m_cols != m.m_rows) {
        throw std::out_of_range("Multiplication: dimension error.");
      }
      Matris temp = *this;
      Matris tempRes(m.cols(), rows());

      for (int y = 0; y < rows(); y++) {
        for (int x = 0; x < m.cols(); x++) {
          tempRes(x,y) = 0;
          for (int x2 = 0; x2 < cols(); x2++) {
            tempRes(x,y) += temp(x2,y) * m(x, x2);

          }
        }
      }
      *this = tempRes;
      return *this;
    }

   // methods
   int cols() const {
      return m_cols ;
    }

    int rows() const {
      return m_rows ;
    }
    double& operator()(int x, int y) {
      return m_vec[x + y * cols()];
    }
    const double& operator()(int x, int y) const {
     return m_vec[x + y * m_cols];
   }

   int validate(){
     for (int row = 0; row < rows(); row++) {
       double row_sum = 0;
       for (int col = 0; col < cols(); col++) {
         row_sum += (*this)(col, row);
       }
       if (abs(row_sum - 1.0) > 0.000001) {
         std::cerr << "Error sum: " << row_sum << '\n';
         return row;
       }
     }
     return -1;
   }
private:
   size_t m_rows;
   size_t m_cols;
   size_t m_capacity;
   double * m_vec;
};

Matris operator* (Matris m1, Matris m2){
  return m1 *= m2;
}

std::ostream &operator<<(std::ostream &os, const Matris &m){
  for (int y = 0; y < m.rows(); y++) {
    for (int x = 0; x < m.cols(); x++) {
      os << m(x,y);
      if (x + 1 != m.cols()) {
        os << " ";
      }
    }
    if (y != m.rows() - 1) {
      os << std::endl;
    }
  }
  return os;
}


std::istream &operator>>(std::istream &in, Matris &m) {
  double value_holder = double();
  int rows;
  int cols;
  in >> rows;
  in.ignore();
  in >> cols;
  in.ignore();
  m = Matris(cols, rows);
  for (int y = 0; y < rows; y++) {
    for (int x = 0; x < cols; x++) {
      in >> value_holder;
      m(x, y) = value_holder;
      in.ignore();
    }
  }
  return in;
}
