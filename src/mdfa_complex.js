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

