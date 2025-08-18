#include <Rcpp.h>
using namespace Rcpp;

//' Test function returning 3
//'
//' This is a simple example of an exported C++ function
//' made available to R.
//'
//' @return An integer equal to 3.
//' @examples
//' test0()
//' @export
// [[Rcpp::export]]
int test0()
{
  return 3;
}
