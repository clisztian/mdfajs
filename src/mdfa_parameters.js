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
    
    
    var Parameters = function(data, L, cutoff, i1, lag, lambda, expweight, decay, decay2, smooth) {
        
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
        
        this.N = this.data[0].length;
        this.fdata = new global.Periodogram(this.data);
        this.gram = new global.Gramm(L, this.N, cutoff, lambda, expweight, lag, this.fdata);
        
         
    };
    
    Parameters.changeRegularization = function(decay, decay2, smooth) {
        
        this.decay = decay; 
        this.decay2 = decay2; 
        this.smooth = smooth;
        
        this.reg_mat = new global.Regularization(this.L, this.lag, this.i1, decay, decay2, smooth);
            
    };
    
    Parameters.changeCustomMSE = function(lambda, expweight) {
       
        this.lambda = lambda; 
        this.expweight = expweight; 
        
        this.gram = new global.Gramm(this.L, this.N, this.cutoff, lambda, expweight, this.lag, this.fdata);
    
    };
    
    Parameters.changeData = function(data) {
        
        this.N = this.data[0].length;
        this.fdata = new global.Periodogram(this.data);
        this.gram = new global.Gramm(this.L, this.N, this.cutoff, this.lambda, this.expweight, this.lag, this.fdata);
        
    };
    
    Parameters.changeFilterLengthLag = function(L, lag) {
        
        this.L = L; 
        this.lag = lag;
        this.reg_mat = new global.Regularization(this.L, this.lag, this.i1, this.decay, this.decay2, this.smooth);
        this.gram = new global.Gramm(this.L, this.N, this.cutoff, this.lambda, this.expweight, this.lag, this.fdata);
        
    };
    
    
    
    global.Parameters = Parameters;
    
})(mdfajs);