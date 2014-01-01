//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

List mue_nelmC(List PITEMLL, List NODW, List Yl, List NU1, int sigmaest, double endest) {
   
   // PITEMLL - item parameter estimates
   // NODW - quadrature nodes and weights
   // Yl - response matrices for each group in List
   // NU1 - null1 matrices for each group in List
   
   int ngru = PITEMLL.size();
   
   
   // create return list
   List riqv_querG(ngru);
   List fiqG(ngru);
   Rcpp::NumericVector mues(ngru);
   Rcpp::NumericVector sigmas(ngru);
   Rcpp::NumericMatrix errmatB(2,ngru); // matrix for the standard errors
   List Xqi_querL(ngru);
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
   double ENDsk;
   
   for(int l = 0; l < listlength; l++)
     { // loops all items
    
    Rcpp::NumericVector PITEM = PITEML[l]; // take out parameters for the l-th item
    IntegerVector y = Y(_,l); // response vector of l-th items
      
    int lpi = PITEM.size()/2; // number of categories
     
     Rcpp::NumericMatrix x(lno,lpi);
     
   for(int o = 0; o < lno; o++)
       {
        double gessum = 0;
        double z2plv = 0;
        double tplf = 0;
        
       for(int q = 0; q < lpi; q++)
         {
         int lpi2 = q+lpi;
         if(q == 0)
         {
        z2plv = exp(PITEM(q) + nodes(o)*PITEM(lpi2)) / ( 1+ exp(PITEM(q) + nodes(o)*PITEM(lpi2)));
        tplf = 1 - z2plv; // 1-P
        //std::cout << "Return" << tplf << " \n ";
         } else {
                x(o,q) = exp(PITEM(q) + nodes(o)*PITEM(lpi2)) * tplf; // mit 1-P multiplizieren
                gessum += exp(PITEM(q) + nodes(o)*PITEM(lpi2)); // new
                }
         }
        x(o,_) = x(o,_) / gessum;
        x(o,0) = z2plv;
       }
  
  
  
    // Zeile 46 des nrm Estep
     //int ysi = y.size();
     
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

    
    // ----- computation of the mean -----
    NumericVector LjmalAmalN(ysi);
    NumericVector Pjischlange(ysi);
    NumericVector Xqi_quer(ysi);
    //double sucolmw;
    
//     for(int col = 0; col < ysi; col++)
//       {
//       sucolmw = sum(ENDm(_,col) * weights);
//       LjmalAmalN = ENDm(_,col) * weights * nodes;
//       //LjmalAmalN = sum(ENDm(_,col) * weights * nodes);
//       ENDsk += LjmalAmalN / sucolmw * 1/ysi; // mean Xqi_qq 
//       } 
//    
//    mues(gru) = ENDsk; // mean of group



    for(int col = 0; col < ysi; col++)
       {
       Pjischlange(col) = sum(ENDm(_,col));
       LjmalAmalN(col)  = sum(ENDm(_,col) * nodes); //Xquer_insum
       Xqi_quer(col)    = LjmalAmalN(col) / Pjischlange(col); // Xqi_quer 
       } 

      mues(gru) = mean(Xqi_quer); //Xqi_qq
      Xqi_querL[gru] = Xqi_quer; // save EAP

  if(sigmaest == 1 & ngru > 1) // if more than one group and an sigma estimation is desired
  {
    
    // ----- computation of sigma -----
    
    NumericVector nodMINmean(lno);
    
    // si_term1
    for(int sma = 0; sma < lno; sma++)
      {
      double zww1 = nodes(sma) - ENDsk;
      nodMINmean(sma) = pow(zww1,2); // zww1 to the power of 2
      }
    
    
    double sigma_xugiZ;
    double sigma_xugi;
    double SIGhatsq = 0;
    
     for(int col = 0; col < ysi; col++)
       {
       sigma_xugiZ = sum(ENDm(_,col) * nodMINmean); // sigma_xugiZ
       sigma_xugi  = sigma_xugiZ / Pjischlange(col);
       SIGhatsq += (sigma_xugi + pow(mues(gru) - Xqi_quer(col),2)) * 1/ysi;
       } 
       
    sigmas(gru) = SIGhatsq; //Xqi_qq
    
  } else 
      {
      sigmas(gru) = 1; // set sigma to 1  
      }

    
// ----- computation of SEs -----
    
  NumericVector serror(2);  
    
  if(endest == 1 & sigmaest == 1)
  //if(endest == 1)
    {
      
      NumericVector mulfac1(lno);
      NumericVector mulfac2(lno);
      
      for(int ean = 0; ean < lno; ean++)
        {
        mulfac1(ean) = (nodes(ean) - mues(gru))/sigmas(gru);
        mulfac2(ean) = (pow(nodes(ean) - mues(gru),2) - sigmas(gru)) / (pow(2*sigmas(gru),2)) ;
        }
    
    double muemue = 0;
    double sigsig = 0;
    double billy = 0;
    double corgan = 0;
    double muesig = 0;
    
    NumericMatrix errmat(2,2);
    

    
    for(int col = 0; col < ysi; col++)
       {
       billy = pow(sum(ENDm(_,col) * mulfac1),2);
       corgan = pow(sum(ENDm(_,col) * mulfac2),2);
       muemue += billy;
       sigsig += corgan;
       muesig += billy * corgan;
       } 
     
  errmat(0,0) = muemue;
  errmat(0,1) = muesig; 
  errmat(1,0) = muesig; 
  errmat(1,1) = sigsig;  
  
  arma::vec serror(2);

  arma::mat Aerrmat = Rcpp::as<arma::mat>(errmat);

  serror = sqrt(diagvec(inv(Aerrmat)));
 
  //Rcout << "Pos 2: " << endest << std::endl; 
  
    } else if(endest >= 1 & sigmaest == 0){
    NumericVector mulfac1(lno);
    double billy = 0;

    
      for(int ean = 0; ean < lno; ean++)
        {
        mulfac1(ean) = (nodes(ean) - mues(gru))/sigmas(gru);
        }
        

    double muemue = 0;
    
    for(int col = 0; col < ysi; col++)
       {
       billy = pow(sum(ENDm(_,col) * mulfac1),2);
       muemue += billy;
       } 
   
   serror(0) = sqrt(1/muemue);
   serror(1) = NA_REAL;
   
    } else 
        {
        //NumericVector serror(2);  
          
         serror(0) = NA_REAL;
         serror(1) = NA_REAL;  
            
        }
        
   errmatB(_,gru) = serror;
   
 // Rcout << "Pos 3: " << gru << std::endl;    

  
// ----- computation of Estepthings -----



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
  
  mues = mues - mues(0);
  sigmas = sigmas + (1-sigmas(0)); //first variance is always 1
  

   
   if(endest == 0)
       {
     return List::create(_["riqv_querG"] = riqv_querG, _["fiqG"] = fiqG,
        _["mean_est"]=mues, _["sig_est"] = sigmas,_["errmat"] = errmatB); 
         
       } else 
           { // for the last estimate return the EAPs as well
         return List::create(_["riqv_querG"] = riqv_querG, _["fiqG"] = fiqG,
            _["mean_est"]=mues, _["sig_est"] = sigmas,_["errmat"] = errmatB,  _["thetas"] = Xqi_querL); 
           }
           

}
