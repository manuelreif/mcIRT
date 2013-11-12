#include <Rcpp.h>
//#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//using namespace arma;



//NumericVector de1nrmC(List ErgE, List PITEMLL, List NODW) {


// [[Rcpp::export]]
Rcpp::NumericVector de1nrmC(List PITEML, Rcpp::NumericVector nodes,NumericMatrix FIQ, NumericMatrix RIQ) {

    int lno = nodes.size(); //number of nodes

    int itzahl = PITEML.size(); // number of items
    int gescat = 0; 
    
    for(int coc = 0; coc < itzahl; coc++) // anzahl der kategorien zählen über all items
      {
      NumericVector wurscht = PITEML[coc];
      gescat += wurscht.size()/2;
      }
    
    //NumericMatrix ZQsternbig(lno,gescat);
//    NumericMatrix Rfo(lno,gescat);
//    NumericMatrix Rfo2(lno,gescat);
      //NumericVector gamderiv(gescat); // vector for 1st derivates of gammas
      //NumericVector xideriv(gescat); // vector for 1st derivates of xis
      NumericVector derivs(gescat*2);

    int endE =0;
    int indexdrivs = 0;
    
    for(int its = 0; its < itzahl; its++)
      { // loops items
        
      Rcpp::NumericVector PITEM = PITEML[its]; // extracts pars for each item
      
      
      int lpi = PITEM.size()/2; // number of categories
      
       Rcpp::NumericMatrix x(lno,lpi);
       
       for(int o = 0; o < lno; o++)
         {
         double gessum = 0;
           
           for(int q = 0; q < lpi; q++)
             { 
             int lpi2 = q+lpi;
             x(o,q) = exp(PITEM(q) + nodes(o)*PITEM(lpi2));
             gessum += x(o,q);
             }

             x(o,_) =  x(o,_) /gessum;
             
         }
       
       
//       for(int blab = 0; blab < lpi; blab++)
//           {
//          ZQsternbig(_,endE)  = x(_,blab);
//           endE += 1;  
//           }

//       for(int blab = 0; blab < lpi; blab++)
//           {
//           Rfo(_,endE)  = RIQ(endE,_) - x(_,blab) * FIQ(_,its);
//           Rfo2(_,endE)  = Rfo(_,endE) * weights;
//           endE += 1;  
//           }

        
       for(int blab = 0; blab < lpi; blab++)
           {
           //gamderiv(endE)  = sum(RIQ(endE,_) - x(_,blab) * FIQ(_,its));
           //xideriv(endE)  = sum((RIQ(endE,_) - x(_,blab) * FIQ(_,its)) * nodes);
           derivs(indexdrivs) = sum(RIQ(endE,_) - x(_,blab) * FIQ(_,its));
           int woxi = indexdrivs + lpi; // position of the xi
           derivs(woxi) = sum((RIQ(endE,_) - x(_,blab) * FIQ(_,its)) * nodes);
           endE += 1;
           indexdrivs += 1;
           }

      Rcout << "The value is " << endE << std::endl;
      Rcout << "The value is " << derivs.size() << std::endl;

        indexdrivs = indexdrivs + lpi;

      }
      
      
      

  //return gamderiv;
//return List::create(_["gamderiv"] = gamderiv, _["xideriv"] = xideriv);

return derivs;

}
