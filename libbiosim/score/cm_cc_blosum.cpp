#include "score/cm_cc_blosum.h"
#include "tools/log.h"

namespace biosim {
  namespace score {
    // default ctor, constructs a blosum62 scoring matrix
    cm_cc_blosum::cm_cc_blosum()
        : _identifier("blosum62"),
          _frequency_matrix_empty(false),
          _frequency_matrix(get_frequency_matrix_blosum62()),
          _bitscore_matrix(get_bitscore_matrix_blosum62()) {}
    // ctor from frequency matrix
    cm_cc_blosum::cm_cc_blosum(std::string const &__identifier, dbl_triangular_matrix const &__frequencies)
        : _identifier(__identifier),
          _frequency_matrix_empty(false),
          _frequency_matrix(__frequencies),
          _bitscore_matrix() {
      int_triangular_matrix bits = to_int_bitscore_matrix(to_bitscore_matrix(_frequency_matrix));
      _bitscore_matrix = std::pair<int_triangular_matrix, double>(bits, 1);
    } // ctor
    // ctor from bitscore matrix
    cm_cc_blosum::cm_cc_blosum(std::string const &__identifier, int_triangular_matrix const &__bitscores,
                               double const &__bitscore_fraction)
        : _identifier(__identifier),
          _frequency_matrix_empty(true),
          _frequency_matrix(get_cc_order().length(), get_cc_order().length()),
          _bitscore_matrix(__bitscores, __bitscore_fraction) {}
    // returns identifier
    std::string cm_cc_blosum::get_identifier() const { return _identifier; }
    // returns the frequency matrix
    cm_cc_blosum::dbl_triangular_matrix const &cm_cc_blosum::get_frequency_matrix() const { return _frequency_matrix; }
    // returns the bitscore matrix
    cm_cc_blosum::int_triangular_matrix const &cm_cc_blosum::get_bitscore_matrix() const {
      return _bitscore_matrix.first;
    } // get_bitscore_matrix()
    // returns the bitscore fraction
    double const &cm_cc_blosum::get_bitscore_fraction() const { return _bitscore_matrix.second; }
    // compares the given two instances of cc
    double cm_cc_blosum::compare(che::cc const &__first, che::cc const &__second) const {
      size_t const first(get_cc_order().find(__first.get_identifier_char()));
      size_t const second(get_cc_order().find(__second.get_identifier_char()));
      return first >= second ? _bitscore_matrix.first(first, second) : _bitscore_matrix.first(second, first);
    } // compare()
    // compares the given two instances of cc using the frequency matrix
    double cm_cc_blosum::frequency_compare(che::cc const &__first, che::cc const &__second) const {
      if(_frequency_matrix_empty) {
        throw compare_error(get_identifier() + ": No frequency matrix data available.");
      } // if

      size_t const first(get_cc_order().find(__first.get_identifier_char()));
      size_t const second(get_cc_order().find(__second.get_identifier_char()));
      return first >= second ? _frequency_matrix(first, second) : _frequency_matrix(second, first);
    } // frequency_compare()

    // get the order in which data for cc is stored in the matrices as string of identifier chars
    std::string cm_cc_blosum::get_cc_order() {
      static std::string cc_order("ARNDCQEGHILKMFPSTWYV");
      return cc_order;
    } // get_cc_order()
    // get the blosum62 frequency matrix
    cm_cc_blosum::dbl_triangular_matrix cm_cc_blosum::get_frequency_matrix_blosum62() {
      static std::vector<double> frequency_data = {
          26782.81, //
          5848.09,  22123.18, //
          4857.14,  4930.09,  17616.72, //
          5400.73,  3953.21,  9269.63,  26504.98, //
          3962.65,  980.88,   1092.50,  994.42,   14864.71, //
          4794.70,  6194.03,  3813.41,  4106.64,  770.86,   9131.13, //
          7444.98,  6711.13,  5505.96,  12248.73, 955.29,   8817.09,  20100.50, //
          14491.39, 4291.00,  7124.25,  6284.36,  1917.42,  3409.41,  4829.16,  47098.64, //
          2759.96,  3091.55,  3563.41,  2376.71,  572.50,   2613.56,  3405.63,  2387.39,  11561.82, //
          7943.72,  3098.65,  2477.68,  3076.61,  2730.06,  2220.18,  3037.90,  3450.79,  1447.42,  22975.42, //
          11009.91, 6028.42,  3411.38,  3787.90,  3907.74,  4030.19,  4990.96,  5198.87,  2459.23,  28361.76, 46272.78,
          8338.88,  15533.06, 6080.29,  6092.99,  1248.97,  7716.37,  10296.37, 6327.07,  2958.56,  3900.98,  6138.15,
          20074.90, //
          3341.90,  2001.08,  1319.20,  1156.88,  939.80,   1843.69,  1692.11,  1825.96,  953.44,   6249.52,  12282.50,
          2264.38,  5042.88, //
          4076.59,  2321.82,  1868.88,  1894.29,  1280.49,  1351.88,  2122.42,  2984.24,  2019.30,  7589.45,  13492.91,
          2364.01,  2965.57,  22771.30, //
          5374.35,  2386.63,  2143.28,  3083.07,  899.80,   2109.68,  3542.26,  3398.98,  1190.29,  2508.69,  3524.79,
          3930.17,  1017.22,  1308.84,  23753.47, //
          15579.38, 5646.43,  7840.40,  6985.53,  2599.50,  4717.08,  7360.29,  9553.90,  2753.84,  4291.98,  6049.14,
          7728.29,  2132.99,  2975.00,  4151.82,  15680.23, //
          9264.29,  4435.93,  5571.65,  4724.65,  2318.32,  3437.77,  5106.27,  5446.48,  1853.12,  6716.05,  8281.94,
          5847.29,  2515.49,  2896.40,  3366.54,  11712.22, 15591.62, //
          1003.70,  662.10,   402.65,   404.11,   360.67,   566.56,   660.02,   1015.14,  377.94,   901.59,   1824.09,
          677.71,   495.17,   2115.89,  352.63,   715.91,   712.02,   8060.58, //
          3239.21,  2308.18,  1745.40,  1491.09,  862.22,   1684.00,  2168.88,  2079.76,  3790.82,  3443.82,  5505.85,
          2489.43,  1423.89,  10562.84, 1126.86,  2566.42,  2346.22,  2211.24,  12765.13, //
          12627.73, 3939.43,  2993.54,  3278.68,  3390.37,  2905.57,  4232.85,  4539.45,  1616.69,  29832.35, 23617.94,
          4824.05,  5761.66,  6419.39,  3102.50,  5877.26,  9070.39,  886.51,   3859.61,  24458.47
          // A      R         N         D         C         Q         E         G         H         I         L
          // K      M         F         P         S         T         W         Y         V
      };

      return dbl_triangular_matrix(20, 20, frequency_data);
    } // get_frequency_matrix_blosum62()
    // get the pair of blosum62 bitscore integer matrix and bitfraction for this matrix
    std::pair<cm_cc_blosum::int_triangular_matrix, double> cm_cc_blosum::get_bitscore_matrix_blosum62() {
      static std::vector<int> bitscore_data = {
          4, //
          -1, 5, //
          -2, 0,  6, //
          -2, -2, 1,  6, //
          0,  -3, -3, -3, 9, //
          -1, 1,  0,  0,  -3, 5, //
          -1, 0,  0,  2,  -4, 2,  5, //
          0,  -2, 0,  -1, -3, -2, -2, 6, //
          -2, 0,  1,  -1, -3, 0,  0,  -2, 8, //
          -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, //
          -1, -2, -3, -4, -1, -2, -3, -4, -3, 2,  4, //
          -1, 2,  0,  -1, -3, 1,  1,  -2, -1, -3, -2, 5, //
          -1, -1, -2, -3, -1, 0,  -2, -3, -2, 1,  2,  -1, 5, //
          -2, -3, -3, -3, -2, -3, -3, -3, -1, 0,  0,  -3, 0,  6, //
          -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, //
          1,  -1, 1,  0,  -1, 0,  0,  0,  -1, -2, -2, 0,  -1, -2, -1, 4, //
          0,  -1, 0,  -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1,  5, //
          -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1,  -4, -3, -2, 11, //
          -2, -2, -2, -3, -2, -1, -2, -3, 2,  -1, -1, -2, -1, 3,  -3, -2, -2, 2,  7, //
          0,  -3, -3, -3, -1, -2, -2, -3, -3, 3,  1,  -2, 1,  -1, -2, -2, 0,  -3, -1,
          4
          // A R  N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y
          // V
      };

      return std::pair<cm_cc_blosum::int_triangular_matrix, double>(int_triangular_matrix(20, 20, bitscore_data), 0.5);
    } // get_bitscore_matrix_blosum62()

    // converts a frequency matrix in a bitscore matrix
    cm_cc_blosum::dbl_triangular_matrix cm_cc_blosum::to_bitscore_matrix(dbl_triangular_matrix const &__frequencies) {
      // total frequencies
      double f_total(0.0);
      for(size_t i = 0; i < __frequencies.size1(); ++i) {
        for(size_t j = 0; j < __frequencies.size2(); ++j) {
          f_total += __frequencies(i, j);
        } // for
      } // for
      // marginal probabilities p_i
      std::vector<double> p(20);
      for(size_t i = 0; i < __frequencies.size1(); ++i) {
        for(size_t j = 0; j < __frequencies.size2(); ++j) {
          p[i] += (__frequencies(i, j) + __frequencies(j, i)) / 2 / f_total;
        } // for
      } // for

      // bit score (in bits)
      boost::numeric::ublas::matrix<double> s(20, 20);
      for(size_t i = 0; i < s.size1(); ++i) {
        for(size_t j = 0; j < s.size2(); ++j) {
          double r((__frequencies(j, i) + __frequencies(i, j)) / 2 / f_total / p[i] / p[j]); // odds ratio
          s(i, j) = std::log(r) / std::log(2); // log odds
        } // for
      } // for

      return boost::numeric::ublas::triangular_adaptor<boost::numeric::ublas::matrix<double>, lower>(s);
    } // to_bitscore_matrix()
    // converts a floating point bitscore matrix in a integer bitscore matrix
    cm_cc_blosum::int_triangular_matrix cm_cc_blosum::to_int_bitscore_matrix(dbl_triangular_matrix const &__bitscores) {
      int_triangular_matrix s(__bitscores.size1(), __bitscores.size2());
      for(size_t i = 0; i < __bitscores.size1(); ++i) {
        for(size_t j = 0; j <= i && j < __bitscores.size2(); ++j) {
          s(i, j) = std::round(__bitscores(i, j));
        } // for
      } // for

      return s;
    } // to_int_bitscore_matrix()
  } // namespace score
} // namespace biosim
