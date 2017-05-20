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

var Parameters = global.Parameters;


function refreshSwatch() {
    
  var lr = $("#sliderL").slider("value");
  Parameters.L = lr;
  $("#Lfilter").html('Filter length: ' + Parameters.L);

}

function refreshCutoff() {
    
  var lr = $("#sliderCutoff").slider("value");
  Parameters.cutoff = lr;
  $("#Cutoff").html('Filter length: ' + Parameters.cutoff);

}

function refreshLambda() {
    
  var lr = $("#sliderLambda").slider("value");
  Parameters.lambda = lr;
  $("#Lambda").html('Filter Anticipation: ' + Parameters.lambda);

}

function refreshExpweight() {
    
  var lr = $("#sliderExpweight").slider("value");
  Parameters.expweight = lr;
  $("#Expweight").html('Filter Smoothness: ' + Parameters.expweight);   

}

function refreshLag() {
    
  var lr = $("#sliderLag").slider("value");
  Parameters.lag = lr;
  $("#Lag").html('Forecast/Smooth: ' + Parameters.lag);   

}

function refreshSmooth() {
    
  var lr = $("#sliderSmooth").slider("value");
  Parameters.smooth = lr;
  $("#regSmooth").html('Smoothness: ' + Parameters.smooth);

}

function refreshDecay() {
    
  var lr = $("#sliderDecay").slider("value");
  Parameters.decay = lr;
  $("#regDecay").html('Decay speed: ' + Parameters.decay);

}

function refreshDecay2() {
    
  var lr = $("#sliderDecay2").slider("value");
  Parameters.decay2 = lr;
  $("#regDecay2").html('Decay start: ' + Parameters.decay2);

}

        
        

$(function() {
    
    // set up slider for learning rate
    $("#sliderL").slider({
      orientation: "horizontal",
      min: 3,
      max: 60,
      step: 1,
      value: 20,
      slide: refreshSwatch,
      change: refreshSwatch
    });
    $("#Lfilter").html('Filter length: ' + Parameters.L);

    // set up slider for learning rate
    $("#sliderCutoff").slider({
      orientation: "horizontal",
      min: Math.PI/30,
      max: Math.PI,
      step: .05,
      value: Math.PI/5,
      slide: refreshCutoff,
      change: refreshCutoff
    });
    $("#Cutoff").html('Filter cutoff: ' + Parameters.cutoff);


    // set up slider for learning rate
    $("#sliderLambda").slider({
      orientation: "horizontal",
      min: 0,
      max: 10,
      step: 0.1,
      value: 0,
      slide: refreshLambda,
      change: refreshLambda
    });
    $("#Lambda").html('Filter Anticipation: ' + Parameters.lambda);
   

    // set up slider for learning rate
    $("#sliderExpweight").slider({
      orientation: "horizontal",
      min: 0,
      max: 10,
      step: 0.1,
      value: 0,
      slide: refreshExpweight,
      change: refreshExpweight
    });
    $("#Expweight").html('Filter Smoothness: ' + Parameters.expweight);    
    
    // set up slider for learning rate
    $("#sliderLag").slider({
      orientation: "horizontal",
      min: -4,
      max: 4,
      step: 1,
      value: 0,
      slide: refreshLag,
      change: refreshLag
    });
    $("#Lag").html('Forecast/Smooth: ' + Parameters.lag);
    
    
    $("#reg").html('Filter Regularization');
    
    $("#sliderSmooth").slider({
      orientation: "horizontal",
      min: 0,
      max: .99,
      step: 0.01,
      value: 0,
      slide: refreshSmooth,
      change: refreshSmooth
    });
    $("#regSmooth").html('Smoothness: ' + Parameters.smooth);
    
    $("#sliderDecay").slider({
      orientation: "horizontal",
      min: 0,
      max: .99,
      step: 0.01,
      value: 0,
      slide: refreshDecay,
      change: refreshDecay
    });
    $("#regDecay").html('Decay speed: ' + Parameters.decay);
    
    $("#sliderDecay2").slider({
      orientation: "horizontal",
      min: 0,
      max: .99,
      step: 0.01,
      value: 0,
      slide: refreshDecay2,
      change: refreshDecay2
    });
    $("#regDecay2").html('Decay start: ' + Parameters.decay2);
});
