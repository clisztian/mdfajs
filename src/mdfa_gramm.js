(function(global)  {
    
    "use strict";
    

    var Gramm = function(L, N, cutoff, lambda, alpha, Lag, eweight) {
      
      var K = N/2;
      var mdfaWeights = new global.Weights(N, cutoff, alpha); 
                  
            
      this.gamma = mdfaWeights.gamma;
      this.rex = numeric.rep([K+1,L],0);
      this.imx = numeric.rep([K+1,L],0);
      this.rows = K+1;
      this.cols = L;
        
      for(var l = 0; l < L; l++) {
       for(var j = 0; j <= K; j++) {       
           
            var arg = ((l-Lag)*Math.PI*j/K); 
            var weight = new Complex(eweight.real[j]*mdfaWeights.expweight[j], eweight.imag[j]*mdfaWeights.expweight[j]);
            var basis = new Complex(0,arg).exp();
            var w = (basis).mul(weight);
            
            this.rex[j][l] = w.re;
            this.imx[j][l] = Math.sqrt(1.0 + mdfaWeights.gamma[j]*lambda)*w.im;
        }
      }  
    };
    
    global.Gramm = Gramm;
    
})(mdfajs);

