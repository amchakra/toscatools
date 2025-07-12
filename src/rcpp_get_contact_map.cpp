#include <Rcpp.h>
using namespace Rcpp;

//' Get contact map using Rcpp
//'
//' This function gets the contact map from a given DataFrame of hybrids. It uses Rcpp for faster computation.
//' It loops through each hybrid, gets the start and end points of the left and right genes, and increments the corresponding cells in the contact map.
//'
//' @param hybrids The DataFrame of hybrids.
//' @param gene_A_size The size of gene A.
//' @param gene_B_size The size of gene B.
//' @return An IntegerMatrix representing the contact map.
//'
//' @export
// [[Rcpp::export]]

 IntegerMatrix rcpp_get_contact_map(DataFrame hybrids, int gene_A_size, int gene_B_size) {

   IntegerVector L_start_v = hybrids["L_start"];
   IntegerVector L_end_v = hybrids["L_end"];
   IntegerVector R_start_v = hybrids["R_start"];
   IntegerVector R_end_v = hybrids["R_end"];

   IntegerMatrix contact_map(gene_A_size, gene_B_size);

   int n_hybrids = hybrids.nrows();
   for (int i = 0; i < n_hybrids; i++) {

     int L_start = L_start_v[i] - 1;
     int L_end = L_end_v[i] - 1;
     int R_start = R_start_v[i] - 1;
     int R_end = R_end_v[i] - 1;

     for (int x = L_start; x <= L_end; x++) {
       for (int y = R_start; y <= R_end; y++) {
         contact_map(x, y) ++;
       }
     }

   }

   return(contact_map);

 }
