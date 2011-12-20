
var FFI = require('node-ffi');

function fortranArrayToJS(fortranArray, m, n) {
    var array = [];
    var a = fortranArray;
    var elementSize = 4;
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

function jsToFortranArray(array) {
    var m = array.length;
    var n = array[0].length;
    var elementSize = 4;
    var fortranArrayStart = fortranArray = new FFI.Pointer(m * n * elementSize);
    for(var i = 0; i < m; i++) {
        for(var j = 0; j < n; j++) {
            fortranArray.putFloat(array[i][j]);
            fortranArray = fortranArray.seek(elementSize);
        }
    }

    return fortranArrayStart;
}

module.exports.fortranArrayToJS = fortranArrayToJS;
module.exports.jsToFortranArray = jsToFortranArray;
