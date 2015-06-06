#ifndef che_score_cm_cc_blosum_h
#define che_score_cm_cc_blosum_h

#include "math/function.h"
#include "che/cc.h"
#include "boost/numeric/ublas/triangular.hpp"

namespace biosim {
  namespace che {
    namespace score {
      // comparison function for comparing two cc based on a blosum scoring matrix; blosum62 is default constructed,
      // any other blosum matrix can be constructed either from its frequency or its bitscore matrix; compare uses the
      // bitscore matrix, frequency_compare uses the frequency matrix.
      // note: the frequency matrix is empty if constructed from a bitscore matrix.
      class cm_cc_blosum : public math::cm_function<che::cc> {
      public:
        using lower = boost::numeric::ublas::lower;
        using row_major = boost::numeric::ublas::row_major;
        using d_triangular_matrix =
            boost::numeric::ublas::triangular_matrix<double, lower, row_major, std::vector<double>>;
        using i_triangular_matrix = boost::numeric::ublas::triangular_matrix<int, lower, row_major, std::vector<int>>;

        // default ctor, constructs a blosum62 scoring matrix
        explicit cm_cc_blosum(bool __use_dbl_bitscore = false);
        // ctor from frequency matrix
        cm_cc_blosum(std::string const &__id, d_triangular_matrix const &__frequencies,
                     bool __use_dbl_bitscore = false);
        // ctor from bitscore matrix
        cm_cc_blosum(std::string const &__id, i_triangular_matrix const &__bitscores, double __bitscore_fraction);
        // returns identifier
        std::string get_identifier() const;
        // returns the frequency matrix
        d_triangular_matrix const &get_frequency_matrix() const;
        // returns if double bitscores are used
        bool get_use_dbl_bitscore() const;
        // returns the bitscore fraction
        double const &get_bitscore_fraction() const;
        // returns the double bitscore matrix
        d_triangular_matrix const &get_dbl_bitscore_matrix() const;
        // returns the int bitscore matrix
        i_triangular_matrix const &get_int_bitscore_matrix() const;
        // compares the given two instances of cc
        double compare(che::cc const &__first, che::cc const &__second) const;

      private:
        std::string _identifier; // id
        d_triangular_matrix _frequency_matrix; // frequency matrix
        bool _use_dbl_bitscore; // if double or int bitscores should be used in compare()
        double _bitscore_fraction; // bitscore fraction
        d_triangular_matrix _dbl_bitscore_matrix; // floating point bitscore matrix
        i_triangular_matrix _int_bitscore_matrix; // int bitscore matrix

        // get the order in which data for cc is stored in the matrices as string of identifier chars
        static std::string get_cc_order();
        // get the blosum62 frequency matrix
        static d_triangular_matrix get_frequency_matrix_blosum62();
        // get the pair of blosum62 bitscore integer matrix and bitfraction for this matrix
        static std::pair<i_triangular_matrix, double> get_bitscore_matrix_blosum62();

        // converts a frequency matrix in a bitscore matrix
        static d_triangular_matrix to_bitscore_matrix(d_triangular_matrix const &__frequencies);
        // converts a floating point bitscore matrix in a integer bitscore matrix
        static i_triangular_matrix to_int_bitscore_matrix(d_triangular_matrix const &__bitscores);
      }; // class cm_cc_blosum
    } // namespace score
  } // namespace che
} // namespace biosim

#endif // che_score_cm_cc_blosum_h
