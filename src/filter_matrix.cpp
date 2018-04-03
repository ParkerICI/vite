#include <Rcpp.h>
using namespace Rcpp;
using std::sort;

// [[Rcpp::export]]

void filter_similarity_matrix_cpp(NumericMatrix m, double threshold) {
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

// [[Rcpp::export]]

void filter_similarity_matrix_by_rank_cpp(NumericMatrix m, unsigned int threshold) {
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

kk <- scgraphs:::filter_similarity_matrix_by_rank(m, 2)

scgraphs:::filter_similarity_matrix_by_rank_cpp(m, 2)


all(kk == m)

*/



