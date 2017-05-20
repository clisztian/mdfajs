(function(global) { 
    
    "use strict";  
    
    var Periodogram = function(data) {
        
        var N = data.length;
        var rab = 0;
        var iab = 0;
        var K = N/2;
        
        this.real = new Array(K+1);
        this.imag = new Array(K+1);
        
        for(var j = 0; j <= K; j++) {
            
         rab = 0; iab = 0;
         for(var i = 0; i < N; i++)  {
            
            var arg = ((i+1)*Math.PI*j/K); 
            rab += Math.cos(arg)*data[i];
            iab += Math.sin(arg)*data[i];
         }
         this.real[j] = rab; this.imag[j] = iab;
        }
    };
    
    
    
global.Periodogram = Periodogram;    
    
    
})(mdfajs);