
(function(global)  {
    
  "use strict";
  var ComplexArray = global.ComplexArray;

  var Mdfa = function() {};

  Mdfa.solve = function solve(gramm, reg_mat) {
    
    var rex = gramm.rex.slice(0);
    var imx = gramm.imx.slice(0);
    var des_mat = reg_mat.des_mat.slice(0);
    var Qsmooth = reg_mat.Q_smooth.slice(0);
          
    var des_matT = numeric.transpose(des_mat); 
    var reXX = numeric.dot(rex,des_matT);
    
    var imXX = numeric.dot(imx,des_matT);
    var XtX = numeric.dot(numeric.transpose(reXX), reXX);
    var imXtX = numeric.dot(numeric.transpose(imXX), imXX);
    XtX = numeric.add(XtX, imXtX);
    
    var dev = numeric.sum(numeric.getDiag(XtX))/reg_mat.w_eight.length;
       
    var des = numeric.dot(Qsmooth, des_matT);
    var regularize = numeric.dot(des_mat, des);
    var btemp = numeric.dot(Qsmooth, reg_mat.w_eight);
    var reg_xtxy = numeric.dot(des_mat, btemp);
    
    if(reg_mat.smooth > 0) {
        var disentangle = numeric.sum(numeric.getDiag(Qsmooth));
        var denom = numeric.sum(numeric.getDiag(regularize));
        
        disentangle = disentangle/denom;
        for(var i=0; i < reg_mat.ncols2; i++) {
                for(var j=0; j < reg_mat.ncols2; j++) {
                  regularize[i][j] *= regularize[i][j]*disentangle*dev;  
              }
        }
        
        for(var i=0; i < reg_mat.ncols; i++) {
                for(var j=0; j < reg_mat.ncols; j++) {
                  Qsmooth[i][j] *= Qsmooth[i][j]*disentangle;
              }
        }
        
    }
        
    var btemp2 = numeric.dot(rex, reg_mat.w_eight);
    var temp2 = numeric.dot(numeric.transpose(reXX), btemp2);
    var btemp3 = numeric.dot(imx, reg_mat.w_eight);
    var xtxy = numeric.dot(numeric.transpose(imXX), btemp3);
    xtxy = numeric.add(xtxy, temp2); 
   
    var b = numeric.dot(numeric.transpose(rex), gramm.gamma);
    b = numeric.sub(b,xtxy);
    
    for(var i=0; i < reg_xtxy.length; i++) {
        reg_xtxy[i] *= reg_xtxy[i]*dev;
    }
    
    b = numeric.sub(b, reg_xtxy);
    XtX = numeric.add(XtX, regularize);

    var lu = numeric.LU(XtX);
    var bh = numeric.LUsolve(lu, b);
    
    return numeric.dot(des_matT, bh);

 };


 Mdfa.applyFilter = function applyFilter(x, coeffs) {
    
    if (x instanceof ComplexArray) {
        var yt = x.real.slice(0);       
    }
    else {
        var yt = x.slice(0); 
    }
    
    
    const L = coeffs.length;
    const n_obs = yt.length;
    var signal = numeric.rep([n_obs - (L-1)],0);
    var xt = numeric.rep([n_obs - (L-1)],0);
    
    var sum = 0;
    for(var i = L-1; i < n_obs; i++) {
       
       sum = 0; 
       for(var l = 0; l < L; l++) { 
           sum = sum + coeffs[l]*yt[i-l]; 
       }   
       signal[i - (L-1)] = sum;
       xt[i - (L-1)] = yt[i];
    }
    
    return {st:signal, xt:xt};
 };
 
 
 Mdfa.computeFilter = function computeFilter(args) {
     
     
     
     
 } 
 
 
 global.Mdfa = Mdfa;
 
})(mdfajs);