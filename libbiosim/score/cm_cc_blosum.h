#ifndef score_cm_cc_blosum_h
#define score_cm_cc_blosum_h

#include "score/function.h"
#include "che/cc.h"
#include "boost/numeric/ublas/triangular.hpp"

namespace biosim {
  namespace score {
    // comparison function for comparing two cc based on a blosum scoring matrix; blosum62 is default constructed,
    // any other blosum matrix can be constructed either from its frequency or its bitscore matrix; compare uses the
    // bitscore matrix, frequency_compare uses the frequency matrix.
    // note: the frequency matrix is empty if constructed from a bitscore matrix.
    class cm_cc_blosum : public cm_function<che::cc> {
    public:
      using lower = boost::numeric::ublas::lower;
      using row_major = boost::numeric::ublas::row_major;
      using dbl_triangular_matrix =
          boost::numeric::ublas::triangular_matrix<double, lower, row_major, std::vector<double>>;
      using int_triangular_matrix = boost::numeric::ublas::triangular_matrix<int, lower, row_major, std::vector<int>>;

      // default ctor, constructs a blosum62 scoring matrix
      cm_cc_blosum();
      // ctor from frequency matrix
      cm_cc_blosum(std::string const &__identifier, dbl_triangular_matrix const &__frequencies);
      // ctor from bitscore matrix
      cm_cc_blosum(std::string const &__identifier, int_triangular_matrix const &__bitscores,
                   double const &__bitscore_fraction);
      // returns identifier
      std::string get_identifier() const;
      // returns the frequency matrix
      dbl_triangular_matrix const &get_frequency_matrix() const;
      // returns the bitscore matrix
      int_triangular_matrix const &get_bitscore_matrix() const;
      // returns the bitscore fraction
      double const &get_bitscore_fraction() const;
      // compares the given two instances of cc
      double compare(che::cc const &__first, che::cc const &__second) const;
      // compares the given two instances of cc using the frequency matrix
      double frequency_compare(che::cc const &__first, che::cc const &__second) const;

    private:
      std::string _identifier; // id
      bool _frequency_matrix_empty; // if _frequency_matrix contains values or not
      dbl_triangular_matrix _frequency_matrix; // frequency matrix
      std::pair<int_triangular_matrix, double> _bitscore_matrix; // pair of bitscore matrix and bitscore fraction

      // get the order in which data for cc is stored in the matrices as string of identifier chars
      static std::string get_cc_order();
      // get the blosum62 frequency matrix
      static dbl_triangular_matrix get_frequency_matrix_blosum62();
      // get the pair of blosum62 bitscore integer matrix and bitfraction for this matrix
      static std::pair<int_triangular_matrix, double> get_bitscore_matrix_blosum62();

      // converts a frequency matrix in a bitscore matrix
      static dbl_triangular_matrix to_bitscore_matrix(dbl_triangular_matrix const &__frequencies);
      // converts a floating point bitscore matrix in a integer bitscore matrix
      static int_triangular_matrix to_int_bitscore_matrix(dbl_triangular_matrix const &__bitscores);
    }; // class cm_cc_blosum
  } // namespace score
} // namespace biosim

#endif // score_cm_cc_blosum_h
