/*
Copyright (c) 2011, Chris Umbel

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

var FFI = require('ffi');

function fortranArrayToJSMatrix(fortranArray, m, n, elementSize) {
    var op = elementSize == 8 ? 'readDoubleLE' : 'readFloatLE';
    var array = [];
    var rowWidth = elementSize * n;
    var columnOffset = m * elementSize;

    for(var i = 0; i < m; i++) {
        var row = [];
        var rowStart = i * elementSize;

        for(var j = 0; j < n; j++) {
            row.push(fortranArray[op](columnOffset * j + rowStart));
        }

        array.push(row);
    }

    return array;
}

function jsMatrixToFortranArray(array, elementSize) {
    var op = elementSize == 8 ? 'writeDoubleLE' : 'writeFloatLE';
    var m = array.length;
    var n = array[0].length;
    var fortranArrayStart = fortranArray = new Buffer(m * n * elementSize);
    for(var j = 0; j < n; j++) {
        for(var i = 0; i < m; i++) {
            fortranArray[op](array[i][j], elementSize * (j * m + i));
        }
    }

    return fortranArrayStart;
}

function fortranArrayToJSArray(fortranArray, n, op, elementSize) {
    var array = [];
    
    for(var i = 0; i < n; i++) {
        array.push(fortranArray[op](i * elementSize));
    }

    return array;
}

module.exports.fortranArrayToJSMatrix = fortranArrayToJSMatrix;
module.exports.jsMatrixToFortranArray = jsMatrixToFortranArray;
module.exports.fortranArrayToJSArray = fortranArrayToJSArray;
