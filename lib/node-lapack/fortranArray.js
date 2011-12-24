
var FFI = require('node-ffi');

var elementSize = 4;

function fortranArrayToJSMatrix(fortranArray, m, n) {
    var array = [];
    var a = fortranArray;
    var rowWidth = elementSize * n;

    for(var i = 0; i < m; i++) {
        var row = [];
        var rowStart = i * elementSize;

        for(var j = 0; j < n; j++) {
            a = fortranArray.seek(rowStart + (j * rowWidth));
            row.push(a.getFloat());
        }

        array.push(row);
    }

    return array;
}

function jsMatrixToFortranArray(array) {
    var m = array.length;
    var n = array[0].length;
    var fortranArrayStart = fortranArray = new FFI.Pointer(m * n * elementSize);
    for(var j = 0; j < n; j++) {
	for(var i = 0; i < m; i++) {
            fortranArray.putFloat(array[i][j]);
            fortranArray = fortranArray.seek(elementSize);
        }
    }

    return fortranArrayStart;
}

function fortranArrayToJSArray(fortranArray, n) {
    var array = [];
    
    for(var i = 0; i < n; i++) {	
	array.push(fortranArray.getFloat());
	fortranArray = fortranArray.seek(elementSize);
    }

    return array;
}

function fortranIntArrayToJSArray(fortranArray, n) {
    var array = [];
    
    for(var i = 0; i < n; i++) {
	array.push(fortranArray.getInt32());
	fortranArray = fortranArray.seek(4);
    }

    return array;
}

module.exports.fortranArrayToJSMatrix = fortranArrayToJSMatrix;
module.exports.jsMatrixToFortranArray = jsMatrixToFortranArray;
module.exports.fortranArrayToJSArray = fortranArrayToJSArray;
module.exports.fortranIntArrayToJSArray = fortranIntArrayToJSArray;
