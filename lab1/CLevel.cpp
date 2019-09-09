#include <cstddef>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <iterator>


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
    ~Matris(){
     delete [] m_vec;
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
      delete [] m_vec;
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

void printMatrix(Matris m) {
  std::cout << m.rows() << ' ';
  std::cout << m.cols() << ' ';
  for (int i = 0; i < m.rows(); i++) {
    for (int j = 0; j < m.cols(); j++) {
      std::cout << m(j, i) << ' ';
    }
  }
  std::cout << "" << '\n';
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


int main(int argc, char const *argv[]) {
  /* code */

  Matris a;
  Matris b;
  Matris pi;

  std::cin >> a;
  std::cin >> b;
  std::cin >> pi;

  int total_observations;
  std::cin >> total_observations;
  int observations[total_observations];
  for (auto &num: observations) {
    std::cin >> num;
  }


  Matris a_original(a);
  Matris b_original(b);
  Matris pi_original(pi);


  // std::cout << a << "\n\n";
  // std::cout << b << "\n\n";
  // std::cout << pi << "\n\n";

  // std::cout.precision(17);


  int total_observations_original = total_observations;
  for (int observation_sample = 1; observation_sample <= total_observations_original; observation_sample++) {
    a = a_original;
    b = b_original;
    pi = pi_original;
    total_observations = observation_sample;

    // double alpha[a.rows()];
    Matris allAlpha(total_observations, a.rows());

    double oldLogProb = -std::numeric_limits<double>::infinity();
    double logProb = -std::numeric_limits<double>::infinity();
    int maxIter = 10000;
    int itertation = 0;

    // while ((std::abs(logProb - oldLogProb) > 0.0000000000001 || itertation == 0) && itertation < maxIter) {
    while (((logProb - oldLogProb) > 0 || itertation == 0) && itertation < maxIter) {

      //compute alpha zero
      double c[total_observations];
      c[0] = 0;
      for (int i = 0; i < a.rows(); i++) {
        allAlpha(0, i) = pi(i,0)*b(observations[0],i);
        c[0] += allAlpha(0, i);
      }
      // scale the alpha zero
      c[0] = 1.0 / c[0];
      for (int i = 0; i < a.rows(); i++) {
        allAlpha(0, i) *= c[0];
      }


      //compute alpha t
      for (int t = 1; t < total_observations; t++) {
        c[t] = 0;
        for (int i = 0; i < a.rows(); i++) {
          allAlpha(t, i) = 0;
          for (int j = 0; j < a.rows(); j++) {
            allAlpha(t, i) += allAlpha(t-1, j) * a(i,j);
          }
          allAlpha(t, i) = allAlpha(t, i) * b(observations[t],i);
          c[t] += allAlpha(t, i);
        }

        // std::cout << "Time step: " << t << '\n';
        // std::cout << "c-value: " << c[t] << '\n';

        // scale the alpha 1 to T
        c[t] = 1.0 / c[t];

        for (int i = 0; i < a.rows(); i++) {
          allAlpha(t, i) *= c[t];
        }
      }



      double resultProb = 0;

      for (int i = 0; i < a.rows(); i++) {
        resultProb += allAlpha(total_observations - 1, i);
      }

      // std::cout << resultProb << '\n';

      Matris allBeta(total_observations, a.rows());
      // calculate Beta pass T
      for (int state = 0; state < a.rows(); state++) {
        allBeta(total_observations - 1, state) = 1 * c[total_observations - 1];
      }

      // std::cout << "Allbeta: " << '\n';
      // std::cout << allBeta << '\n';

      for (int time_step = total_observations - 2; time_step >= 0; time_step--) {
        for (int i = 0; i < a.rows(); i++) {
          for (int j = 0; j < a.rows(); j++) {
            allBeta(time_step, i) += a(j, i) * b(observations[time_step + 1], j) * allBeta(time_step + 1, j);
          }
          // std::cout << "Beta: " << allBeta(time_step, i) << '\n';
          // Scale beta t with same factor as alpha t
          allBeta(time_step, i) *= c[time_step];
        }
      }


      //Compute sum of alpha T
      double denom = 0;
      for (int state = 0; state < a.rows(); state++) {
        denom += allAlpha(total_observations-1, state);
      }

      //compute gamma i,j and gamma i
      Matris digamma[total_observations];
      Matris gamma(total_observations,a.rows());

      // std::cout << "allAlpha: " << allAlpha << '\n';
      // std::cout << "allBeta: " << allBeta << '\n';
      // for (auto i: c) {
      //   std::cout <<  i << '\n';
      // }

      for (int time_step = 0; time_step < total_observations - 1; time_step++) {
        digamma[time_step] = Matris(a.rows(), a.rows());
        for (int i = 0; i < a.rows(); i++) {
          for (int j = 0; j < a.rows(); j++) {
            digamma[time_step](j,i) = allAlpha(time_step, i)*a(j,i)*b(observations[time_step+1], j)*allBeta(time_step+1, j)/denom;
            gamma(time_step, i) += digamma[time_step](j,i);
          }
        }
      }


      // special case for gamma T
      denom = 0;
      for (int i = 0; i < a.rows(); i++) {
        denom += allAlpha(total_observations-1, i);
      }
      for (int i = 0; i < a.rows(); i++) {
        gamma(total_observations-1, i) = allAlpha(total_observations-1, i)/denom;
      }

      // Re-estimate A, B and pi

      // Re-estimate pi
      for (int i = 0; i < a.rows(); i++) {
        pi(i,0) = gamma(0, i);
      }

      // Re-estimate A
      for (int i = 0; i < a.rows(); i++) {
        for (int j = 0; j < a.rows(); j++) {
          double numer = 0;
          double denom = 0;
          for (int time_step = 0; time_step < total_observations - 1; time_step++) {
            numer = numer + digamma[time_step](j, i);
            denom = denom + gamma(time_step, i);
          }
          // std::cout << "Estimate A step("<< i << " " << j <<  "): denom: " << denom << " numer: " << numer << '\n';
          // if (denom < 0.0000000000001 && denom > -0.0000000000001 ) {
          //  a(j,i) = 0;
          // } else {
            a(j,i) = numer/denom;
          // }
        }
      }

      // Re-estimate B
      for (int j = 0; j < a.rows(); j++) {
        for (int k = 0; k < b.cols(); k++) {
          double numer = 0;
          double denom = 0;
          for (int time_step = 0; time_step < total_observations - 1; time_step++) {
            if (observations[time_step] == k) {
              numer = numer + gamma(time_step, j);
            }
            denom = denom + gamma(time_step, j);
          }
          // if (denom < 0.000000000000001 && denom > -0.0000000000001 ) {
          //   b(k,j) = 0;
          // } else {
            b(k,j) = numer/denom;
          // }
        }
      }

      // Compute log[P(O|λ)]
      oldLogProb = logProb;
      logProb = 0;
      for (int time_step = 0; time_step < total_observations; time_step++) {
        logProb -= log(c[time_step]);
      }
      itertation++;
      // std::cout << "Iteration: " << itertation << '\n';
      // std::cout << "LogProb:   " << logProb << '\n';
    }

    // std::cout << "/* message */" << '\n';
    // for (int i = 0; i < total_observations; i++) {
    //   std::cout << c[i] << ' ';
    // }
    // std::cout << "/* message */" << '\n';

    std::cout << "" << '\n';
    if (itertation == maxIter) {
      std::cout << "#################################################" << '\n';
    }
    std::cout << "Observations " << observation_sample << '\n';
    std::cout << "Iterations: " << itertation << '\n';
    std::cout << "A: \n" << a << '\n';
    std::cout << "B: \n"<< b << '\n';
    std::cout << "" << '\n';
  }

  // printMatrix(a);
  // printMatrix(b);

  return 0;
}
