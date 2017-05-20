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
