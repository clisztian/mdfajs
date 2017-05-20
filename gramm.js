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

"use strict";

var Complex = require('complex.js');
var Weights = require('./weights');

var numeric = require('numeric');
const PI = Math.PI;
const SQRT1_2 = Math.SQRT1_2;


class Gramm {
    
    constructor(L, N, cutoff, lambda, alpha, Lag, eweight) {
      
      var K = N/2;
      const mdfaWeights = new Weights(N, cutoff, alpha); 
            
      this.gamma = mdfaWeights.gamma;
      this.rex = numeric.rep([K+1,L],0);
      this.imx = numeric.rep([K+1,L],0);
      this.rows = K+1;
      this.cols = L;
        
      for(var l = 0; l < L; l++) {
       for(var j = 0; j <= K; j++) {       
           
            var arg = ((l-Lag)*PI*j/K); 
            var weight = new Complex(eweight.real[j]*mdfaWeights.expweight[j], eweight.imag[j]*mdfaWeights.expweight[j]);
            var basis = new Complex(0,arg).exp();
            var w = (basis).mul(weight);
            
            this.rex[j][l] = w.re;
            this.imx[j][l] = Math.sqrt(1.0 + mdfaWeights.gamma[j]*lambda)*w.im;
        }
      }  
    };
    
}

module.exports = Gramm;
