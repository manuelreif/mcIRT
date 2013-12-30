//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

List EnelmC(List PITEMLL, List NODW, List Yl, List NU1) {
   
   // PITEMLL - item parameter estimates
   // NODW - quadrature nodes and weights
   // Yl - response matrices for each group in List
   // NU1 - null1 matrices for each group in List
   
   int ngru = PITEMLL.size();
   
   // create return list
   List riqv_querG(ngru);
   List fiqG(ngru);
   
   // the group loop
   for(int gru = 0; gru < ngru; gru++)
   {
    // extract everything out of the Lists 
    List PITEML = PITEMLL[gru];
    List nodeswei = NODW[gru];
    Rcpp::IntegerMatrix Y = Yl[gru];
    Rcpp::IntegerMatrix nu1m = NU1[gru]; 
     
    NumericVector nodes =  nodeswei[0];
    NumericVector weights =  nodeswei[1];
    
   int listlength = PITEML.size(); // number of items
   int ysi = Y.nrow(); // number of observations per item (including NA's)
   int lno = nodes.size(); // number of quadrature nodes
   
   // create matrix outside the loop
   Rcpp::NumericMatrix ENDm(lno,ysi); // matrix with proper dimensions for multiplication
   ENDm.fill(1); // write 1 in each cell
   
   for(int l = 0; l < listlength; l++)
     { // loops all items
    
    Rcpp::NumericVector PITEM = PITEML[l]; // take out parameters for the l-th item
    IntegerVector y = Y(_,l); // response vector of l-th items
      
    int lpi = PITEM.size()/2; // number of categories
    int lpim1 = lpi - 1; 
     
     //arma::mat x(lno,lpi);
     Rcpp::NumericMatrix x(lno,lpi);
     
     // das hier ist zeile 40 bis 44 des nrm Estep
     for(int o = 0; o < lno; o++)
       {
        double gessum = 0;
        double z2plv = 0;
        double tplf = 0;
        
       for(int q = 0; q < lpi; q++)
         {
         int lpi2 = q+lpi;
         if(q == 0)
         { // hier muss jetzt das 2pl Modell rein
         // exp(Km %*% abpar) / (1 + exp(Km %*% abpar)) 
        // x(o,q) = exp(PITEM(q) + nodes(o)*PITEM(lpi2)) / ( 1+ exp(PITEM(q) + nodes(o)*PITEM(lpi2))); // 2pl
        z2plv = exp(PITEM(q) + nodes(o)*PITEM(lpi2)) / ( 1+ exp(PITEM(q) + nodes(o)*PITEM(lpi2)));
        tplf = 1 - z2plv; // 1-P
        //std::cout << "Return" << tplf << " \n ";
         } else {
                x(o,q) = exp(PITEM(q) + nodes(o)*PITEM(lpi2)) * tplf; // mit 1-P multiplizieren
                }
         gessum += exp(PITEM(q) + nodes(o)*PITEM(lpi2)); // new
         }
          
        x(o,_) = x(o,_) / gessum;
        x(o,0) = z2plv;
        
        //std::cout << "z2plv:" << z2plv << " \n ";
        // std::cout << "inmattrix" << x(o,0) << " \n ";
        //arma::mat zwischen = x(o,arma::span(1,lpim1)) / gessum;
        //x(o,arma::span(1,lpim1)) = zwischen; // nrm cats divided by the gessum
         
       }
  
     
     Rcpp::NumericMatrix z(lno,ysi);
     
     for(int i = 0; i < ysi; i++)
       {
        int whichE = y(i);
        
        // if there is NOT a missing value, make standard procedure
        // else (missing value is there) multiply with 1 - that means make a copy of what was there before
        //if(!NumericVector::is_na(whichE))
        if(!IntegerVector::is_na(whichE))
          {
            z(_,i) = x(_,whichE) * ENDm(_,i); // at the end ENDm will be the product again
          } else {
                  z(_,i) =  ENDm(_,i);
                 }
       }
         
      ENDm = z;       
       
     }
    
    NumericVector colmw; 
    
     for(int col = 0; col < ysi; col++)
       {
       colmw = ENDm(_,col) * weights;
       ENDm(_,col) = colmw / sum(colmw); // normalize
       } // das muss ja fiq sein


    /////// 
  
    arma::mat Anu1m = Rcpp::as<arma::mat>(nu1m);
    arma::mat AENDm = Rcpp::as<arma::mat>(ENDm);
    
    
    //arma::mat riqv_quer = Anu1m * trans(AENDm); 
    arma::mat riqv_quer = trans(Anu1m) * trans(AENDm); 
    
    riqv_querG[gru] = riqv_quer; // save in list
    
    
    //NumericVector fiq(lno);
    // calculate fiq which is the expected number of persons on each node
    
    
    IntegerMatrix fivor(ysi,listlength);
    
    // write 0 if missing value, 1 if valid response
    for(int ww = 0; ww < listlength; ww++)
      {
      fivor(_,ww) = ifelse(is_na(Y(_,ww)),0,1);  
      }
    
    arma::mat Afivor = Rcpp::as<arma::mat>(fivor);
    arma::mat fiq = AENDm * Afivor;
    fiqG[gru] = fiq;
    
    

//     
  } // end of group loop
   
   
    return List::create(_["riqv_querG"] = riqv_querG, _["fiqG"] = fiqG);
     //return riqv_querG;
}
