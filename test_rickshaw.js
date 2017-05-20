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
var numeric = require('numeric');
var ComplexArray = require('./complex_array');
var Rickshaw = require('rickshaw');

var mdfa = require('./mdfa');
var ComplexArray = require('./complex_array');
var Weights = require('./weights');
var Gramm = require('./gramm');
var Regularization = require('./regularization');
var numeric = require('numeric');


const N = 150;
const cutoff = Math.PI/5;
const lambda = 0.0;
const alpha = 1.0;
const L = 20;
const lag = 0;
const i1 = 1;
const decay = 0.0;
const decay2 = 0.0;
const smooth = 0.0;

var seriesData = [ [], [] ];
var random = new Rickshaw.Fixtures.RandomData(N);

for (var i = 0; i < N; i++) {
	random.addData(seriesData);
}


var data = new ComplexArray(N);

for(var i = 0; i < seriesData[0].length; i++) {
    console.log(seriesData[0][i].y);
    data.real[i] = seriesData[0][i].y;
}

var xt = data.real.slice(0);

//-- Fourier transform of data--
var fdata = data.FFT();

//-- Compute Grammian matrix--
var gram = new Gramm(L, N, cutoff, lambda, alpha, lag, fdata);

//-- Compute regularization matrix--
var reg_mat = new Regularization(L, lag, i1, decay, decay2, smooth);

//-- Assemble matrices and solve--
var coeffs = mdfa.solve(gram, reg_mat);
console.log(numeric.sum(coeffs));

var mysig = mdfa.applyFilter(xt, coeffs);
console.log(mysig);
