#include <Rcpp.h>
using namespace Rcpp;
using std::sort;


//' Filter a matrix based on a threshold value
//'
//' This function filters an input matrix based on a threshold value. Values below
//' the threshold are set to 0. For each row this function calculates \code{max(row)}.
//' The actual threshold used for filtering is \code{min(threshold, max(row))}. This guarantees that
//' there is always at least one non-zero value in each row. WARNING: the input matrix is modified
//' in-place.
//'
//' @param m The input matrix
//' @param treshold The threhsold used for filtering. See Description for details
//' @return This function does not return any value. The input matrix is modified in-place.
//'
//'
// [[Rcpp::export]]

void filter_matrix(NumericMatrix m, double threshold) {
    unsigned int nrow = m.nrow();
    unsigned int ncol = m.ncol();

    for(unsigned int i = 0; i < nrow; i++) {
        double row_max = m(i, 0);

        for(unsigned int j = 1; j < ncol; j++)
            if(m(i, j) > row_max)
                row_max = m(i, j);

        double t = row_max < threshold ? row_max : threshold;

        for(unsigned int j = 0; j < ncol; j++)
            if(m(i, j) < t)
                m(i, j) = 0.0;

    }
}


class Comparator {

    private:
        const NumericVector& v;

    public:
        Comparator(const NumericVector &v_) : v(v_) {}


        bool operator() (const int lhs, const int rhs) {return v[lhs] < v[rhs];}
};




//' Filter a matrix based on the rank of the values in each row
//'
//' This function filters an input matrix, by ranking the data in each row from largest to smallest,
//' and setting any element whose rank is greater than the input threshold to 0. In other words,
//' if the treshold is X, only the X greatest values in each row will be kept, and the rest will
//' be set to 0. WARNING: the input matrix is modified in-place.
//'
//' @param m The input matrix
//' @param threshold The threshold rank. Values with rank greater than the threhsold will be set to
//'   0. Note that the rank is 1-based (i.e. the largest observation has rank 1)
//' @return This function does not return any value. The input matrix is modified in-place.
//' @export
//'
// [[Rcpp::export]]
void filter_matrix_by_rank(NumericMatrix m, unsigned int threshold) {
    unsigned int nrow = m.nrow();
    IntegerVector v = seq(0, m.ncol() - 1);


    for(unsigned int i = 0; i < nrow; i++) {
        NumericVector row_v = m(i, _);

        sort(v.begin(), v.end(), Comparator(row_v));

        for(unsigned int idx = 0; idx < v.size(); idx++) {
            unsigned int rank = v.size() - idx;
            if(rank > threshold)
                m(i, v[idx]) = 0.0;
        }
    }
}

/*

m <- matrix(nrow = 9000, ncol = 9000, data = rnorm(9000 ^ 2))

kk <- vite:::filter_similarity_matrix_by_rank(m, 2)

vite:::filter_similarity_matrix_by_rank_cpp(m, 2)


all(kk == m)

*/



