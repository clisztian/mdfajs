/* 
 * The MIT License
 *
 * Copyright 2017 clisztian.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

var mdfajs = mdfajs || { REVISION: 'ALPHA' };
/* 
 * The MIT License
 *
 * Copyright 2017 clisztian.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

(function(global) {
  
    
  "use strict";
  const PI = Math.PI;
  const SQRT1_2 = Math.SQRT1_2;


  var ComplexArray = function(other) {
 
    if (other instanceof ComplexArray) {
      // Copy constuctor.
      this.real = new Array(other.real);
      this.imag = new Array(other.imag);
    } else {
      // other can be either an array or a number.
      this.real = new Array(other);
      this.imag = new Array(this.real.length);
    }

    this.length = this.real.length;
  };

  
  ComplexArray.prototype = {
       
   
  
  FFT: function() {
    return ComplexArray.fft(this, false);
  },
  
  InvFFT: function() {
    return ComplexArray.fft(this, true);
  }
 
 
 };


 ComplexArray.fft = function(input, inverse) {
  const n = input.length;

  if (n & (n - 1)) {
    return ComplexArray.FFT_Recursive(input, inverse);
  } else {
    return ComplexArray.FFT_2_Iterative(input, inverse);
  }
 };

  ComplexArray.FFT_Recursive = function(input, inverse) {
   const n = input.length;

   if (n === 1) {
    return input;
   }

   const output = new ComplexArray(n);

   // Use the lowest odd factor, so we are able to use FFT_2_Iterative in the
   // recursive transforms optimally.
   const p = LowestOddFactor(n);
   const m = n / p;
   const normalisation = 1 / Math.sqrt(p);
   var recursive_result = new ComplexArray(m, input.ArrayType);

   // Loops go like O(n Σ p_i), where p_i are the prime factors of n.
   // for a power of a prime, p, this reduces to O(n p log_p n)
   for(var j = 0; j < p; j++) {
    for(var i = 0; i < m; i++) {
      recursive_result.real[i] = input.real[i * p + j];
      recursive_result.imag[i] = input.imag[i * p + j];
    }
    // Don't go deeper unless necessary to save allocs.
    if (m > 1) {
      recursive_result = ComplexArray.fft(recursive_result, inverse);
    }

    const del_f_r = Math.cos(2*PI*j/n);
    const del_f_i = (inverse ? -1 : 1) * Math.sin(2*PI*j/n);
    var f_r = 1;
    var f_i = 0;

    for(var i = 0; i < n; i++) {
      const _real = recursive_result.real[i % m];
      const _imag = recursive_result.imag[i % m];

      output.real[i] += f_r * _real - f_i * _imag;
      output.imag[i] += f_r * _imag + f_i * _real;
      
      f_r = f_r * del_f_r - f_i * del_f_i;
      f_i = f_r * del_f_i + f_i * del_f_r;

    }
  }

  // Copy back to input to match FFT_2_Iterative in-placeness
  // TODO: faster way of making this in-place?
  for(var i = 0; i < n; i++) {
    input.real[i] = normalisation * output.real[i];
    input.imag[i] = normalisation * output.imag[i];
  }

  return input;
};


ComplexArray.FFT_2_Iterative = function(input, inverse) {
  const n = input.length;

  const output = BitReverseComplexArray(input);
  const output_r = output.real;
  const output_i = output.imag;
  // Loops go like O(n log n):
  //   width ~ log n; i,j ~ n
  var width = 1;
  while (width < n) {
    const del_f_r = Math.cos(PI/width);
    const del_f_i = (inverse ? -1 : 1) * Math.sin(PI/width);
    for (var i = 0; i < n/(2*width); i++) {
      var f_r = 1;
      var f_i = 0;
      for (var j = 0; j < width; j++) {
        const l_index = 2*i*width + j;
        const r_index = l_index + width;

        const left_r = output_r[l_index];
        const left_i = output_i[l_index];
        const right_r = f_r * output_r[r_index] - f_i * output_i[r_index];
        const right_i = f_i * output_r[r_index] + f_r * output_i[r_index];

        output_r[l_index] = SQRT1_2 * (left_r + right_r);
        output_i[l_index] = SQRT1_2 * (left_i + right_i);
        output_r[r_index] = SQRT1_2 * (left_r - right_r);
        output_i[r_index] = SQRT1_2 * (left_i - right_i);
        
        f_r = f_r * del_f_r - f_i * del_f_i;
        f_i = f_r * del_f_i + f_i * del_f_r;

      }
    }
    width <<= 1;
  }

  return output;
};

function BitReverseIndex(index, n) {
  var bitreversed_index = 0;

  while (n > 1) {
    bitreversed_index <<= 1;
    bitreversed_index += index & 1;
    index >>= 1;
    n >>= 1;
  }
  return bitreversed_index;
}

function BitReverseComplexArray(array) {
  const n = array.length;
  const flips = new Set();

  for(var i = 0; i < n; i++) {
    const r_i = BitReverseIndex(i, n);

    if (flips.has(i)) continue;
    
    var reali = array.real[i];
    array.real[i] = array.real[r_i];
    array.real[r_i] = reali;
    
    var imagi = array.imag[i];
    array.imag[i] = array.imag[r_i];
    array.imag[r_i] = imagi;

    flips.add(r_i);
  }

  return array;
}

function LowestOddFactor(n) {
  const sqrt_n = Math.sqrt(n);
  var factor = 3;

  while(factor <= sqrt_n) {
    if (n % factor === 0) return factor;
    factor += 2;
  }
  return n;
}

global.ComplexArray = ComplexArray;

})(mdfajs);

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
/* 
 * The MIT License
 *
 * Copyright 2017 clisztian.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


(function(global) { 
    
    "use strict";
    
    
    var Parameters = function(L, cutoff, i1, lag, lambda, expweight, decay, decay2, smooth) {
        
        this.L = L; 
        this.cutoff = cutoff; 
        this.i1 = i1; 
        this.lag = lag; 
        this.lambda = lambda; 
        this.expweight = expweight; 
        this.decay = decay; 
        this.decay2 = decay2; 
        this.smooth = smooth; 
        
        this.reg_mat = new global.Regularization(L, lag, i1, decay, decay2, smooth); 
    };
    
    Parameters.changeRegularization = function(L, lag, i1, decay, decay2, smooth) {
        
        this.L = L;
        this.i1 = i1; 
        this.lag = lag; 
        this.decay = decay; 
        this.decay2 = decay2; 
        this.smooth = smooth;
        
        this.reg_mat = new global.Regularization(L, lag, i1, decay, decay2, smooth);
            
    };
    
    
    global.Parameters = Parameters;
    
})(mdfajs);
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

/* 
 * The MIT License
 *
 * Copyright 2017 clisztian.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * The MIT License
 *
 * Copyright 2017 clisztian.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

(function(global) {

"use strict";

  var Weights = function(N, cutoff, alpha) {
        
        const K = N/2;
        this.cutoff = cutoff;
        this.N = N;
        this.gamma = new Array(K+1);
        this.expweight = new Array(K+1);
        
        
        for(var i = 0; i < this.gamma.length; i++) {
            var loc = i*Math.PI/K;
            this.gamma[i] = (loc < cutoff) ? 1 : 0;
            this.expweight[i] = (loc < cutoff) ? 1 : Math.pow(loc + 1, alpha/10);
        }
 };
 global.Weights = Weights;

})(mdfajs);



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
/* 
 * The MIT License
 *
 * Copyright 2017 clisztian.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

(function(global)  {


  "use strict";

  
  var Regularization = function(L, lag, i1, decay, decay2, smooth) {
        
        this.ncols = L;
        this.ncols2 = (i1 === 1) ? L-1 : L; 
        
        var ncols = this.ncols; 
        var ncols2 = this.ncols2;
        
        const lagi = Math.floor(lag);
        
        var Q_decay = numeric.rep([ncols,ncols],0);
        this.des_mat = numeric.rep([ncols2,ncols],0);
        this.Q_smooth = numeric.rep([ncols,ncols],0);
        this.w_eight = numeric.rep([ncols],0);

        decay =  Math.tan(Math.min(decay,0.999999)*Math.PI/2.0);
        decay2 =  100.0*Math.tan(Math.min(decay2,0.999999)*Math.PI/2.0);
        smooth = 100.0*Math.tan(Math.min(smooth,0.999999)*Math.PI/2.0);
                       
        this.Q_smooth[0][0] = smooth;  
        this.Q_smooth[0][1] = -2.0*smooth;
        this.Q_smooth[0][2] = smooth;  
        Q_decay[0][0]  = decay2*Math.pow(1 + decay,  2.0*Math.abs(-lag)); 
 
        this.Q_smooth[1][0] = -2.0*smooth; 
        this.Q_smooth[1][1] = 5.0*smooth;
        this.Q_smooth[1][2] = -4.0*smooth;
        this.Q_smooth[1][3] = smooth;
        Q_decay[1][1] = decay2*Math.pow(1 + decay,  2.0*Math.abs(1.0-lag)); 

        var l1 = L-1;
        this.Q_smooth[l1][l1-2] =  smooth;       
        this.Q_smooth[l1][l1-1] = -2.0*smooth;     
        this.Q_smooth[l1][l1] =  1.0*smooth;    
        Q_decay[l1][l1] =  decay2*Math.pow(1 + decay, 2.0*Math.abs(1.0*l1-lag));        
        
        l1 = L-2;
        this.Q_smooth[l1][l1-2] =  smooth;       
        this.Q_smooth[l1][l1-1] = -4.0*smooth;     
        this.Q_smooth[l1][l1] =    5.0*smooth;    
        this.Q_smooth[l1][l1+1] =   -2.0*smooth;    
        Q_decay[l1][l1] =  decay2*Math.pow(1 + decay, 2.0*Math.abs(1.0*l1-lag));      
        
        for(var i=2;i < L-2; i++) {
            
            Q_decay[i][i] = decay2*Math.pow(1 + decay, 2.0*Math.abs(1.0*i-lag));
            this.Q_smooth[i][i-2] =  smooth;
            this.Q_smooth[i][i-1] = -4.0*smooth;
            this.Q_smooth[i][i]   =  6.0*smooth;
            this.Q_smooth[i][i+1] = -4.0*smooth; 
            this.Q_smooth[i][i+2] =  smooth; 
        }
        
        
        var trace = 0.000000000000001; 
        var strace = 0.000000000000001; 
        for(var i=0; i < L; i++) {
            trace = trace + Q_decay[i][i];  
            strace = strace + this.Q_smooth[i][i];
        }
    
   
        if(decay2 > 0 || smooth > 0) {            
            for(var i=0; i < L; i++) {
                for(var j=0; j < L; j++) {
                  Q_decay[i][j] *= Q_decay[i][j]*decay2/trace;
                  this.Q_smooth[i][j] *= this.Q_smooth[i][j]*smooth/strace;
                  this.Q_smooth[i][j] += Q_decay[i][j];
              }
            }
        }
               
        if(i1 === 0) {
            for(var i = 0; i < L; i++) {
                this.des_mat[i][i] = 1;
            }     
                       
            if (lag < 1) {
                this.w_eight[0] = 1;
            }
            else {
                this.w_eight[lagi] = 1;
            }  
        }
        else {
            var il = 0;
            for(var i=0; i < L-1; i++) {
              if (lag<1) {
                this.des_mat[i][i+1] = 1;
                this.des_mat[i][0]    = -1;
              }
              else {
                  il = (i >= lag) ? i+1 : i;
                  this.des_mat[i][il] = 1; 
	          this.des_mat[i][lagi] = -1;                  
              }
          }
        }
    };
 
    global.Regularization = Regularization;
 
 })(mdfajs);
(function(lib) {
  "use strict";
  if (typeof module === "undefined" || typeof module.exports === "undefined") {
    window.mdfajs = lib; // in ordinary browser attach library to window
  } else {
    module.exports = lib; // in nodejs
  }
})(mdfajs);
