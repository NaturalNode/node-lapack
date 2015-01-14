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

var fortranArray = require('./fortranArray');
var FFI = require('node-ffi');

var LAPACK;

try {
    LAPACK = new FFI.Library('liblapack', {
	"sgeqrf_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", 
			    "pointer", "pointer", "pointer"]],
	"dgeqrf_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", 
			    "pointer", "pointer", "pointer"]],
	"sorgqr_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", 
			    "pointer", "pointer", "pointer"]],
	"sgesvd_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", 
			    "pointer", "pointer", "pointer", "pointer", "pointer", 
			    "pointer", "pointer", "pointer", "pointer", ]],
	"sgetrf_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
	"spotrf_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer"]],
	"dgetrf_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
	"sgesv_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]]
	
    });
} catch(e) {
    console.log("!!! node-lapack requires the native lapack to be built as a shared lib.");
    console.log(e);
}

var FORTRAN_INT = 4;
var FORTRAN_CHAR = 1;
var FORTRAN_FLOAT = 4;
var FORTRAN_DOUBLE = 8;

function eye(m) {
    var matrix = [];

    for(var i = 0; i < m; i++) {
	var row = [];
	matrix.push(row);
	
	for(var j = 0; j < m; j++) {
	    if(i == j)
		row.push(1);
	    else
		row.push(0);
	}
    }

    return matrix;
}

function matrixOp(matrix, elementSize, callback) {
    var m = matrix.length;
    var n = matrix[0].length;    
    var f_m = new FFI.Pointer(FORTRAN_INT);
    var f_n = new FFI.Pointer(FORTRAN_INT);
    var f_a = fortranArray.jsMatrixToFortranArray(matrix, elementSize);
    var f_lda = new FFI.Pointer(FORTRAN_INT);
    var f_letter_u = new FFI.Pointer(FORTRAN_CHAR);
    var f_letter_l = new FFI.Pointer(FORTRAN_CHAR);
    
    f_letter_u.putChar(85);
    f_letter_l.putChar(76);
    f_m.putInt32(m);
    f_n.putInt32(n);
    f_lda.putInt32(Math.max(1, m));
   console.log(f_letter_u); 
    callback(m, n, f_m, f_n, f_a, f_lda, f_letter_u, f_letter_l);
}

function zeroBottomLeft(matrix) {
    // zero out bottom left forming an upper right triangle matrix                
    for(var i = 1; i < matrix.length; i++) {
        for(var j = 0; j < i && j < matrix[0].length; j++)
            matrix[i][j] = 0;
    }

    return matrix
}

function sgesv(a, b) {
    var f_info = new FFI.Pointer(FORTRAN_INT);
    var result = {};

    matrixOp(a, FORTRAN_FLOAT, function(am, an, af_m, af_n, f_a) {
	var f_ipiv = new FFI.Pointer(am  * FORTRAN_INT);

	matrixOp(b, FORTRAN_FLOAT, function(bm, bn, bf_m, bf_n, f_b) {
	    LAPACK.sgesv_(af_m, bf_n, f_a, af_n, f_ipiv, f_b, bf_m, f_info);
	    result.X = fortranArray.fortranArrayToJSMatrix(f_b, bm, bn, FORTRAN_FLOAT);
	    result.P = ipivToP(bm, fortranArray.fortranArrayToJSArray(f_ipiv, bm, 'getInt', FORTRAN_FLOAT));
	});
    });

    return result;
}

function qr(matrix) {
    var result;

    sgeqrf(matrix, function(qr, m, n, f_m, f_n, f_a, f_lda, f_tau, f_work, f_lwork, f_info) {
	var f_k = new FFI.Pointer(FORTRAN_INT);
	f_k.putInt32(Math.min(m, n));
	LAPACK.sorgqr_(f_m, f_n, f_k, f_a, f_lda, f_tau, f_work, f_lwork, f_info);
	qr.Q = fortranArray.fortranArrayToJSMatrix(f_a, m, n, FORTRAN_FLOAT);
	qr.R = zeroBottomLeft(qr.R);
	result = qr;
    });
    
    return result;
}

function sgeqrf(matrix, callback) {
    return geqrf(matrix, 'getFloat', LAPACK.sgeqrf_, FORTRAN_FLOAT, callback);
}

function dgeqrf(matrix, callback) {
    return geqrf(matrix, 'getDouble', LAPACK.dgeqrf_, FORTRAN_DOUBLE, callback);
}

function geqrf(matrix, op, lapackFunc, elementSize, callback) {
    var qr;
    
    matrixOp(matrix, elementSize, function(m, n, f_m, f_n, f_a, f_lda) {
	var f_tau = new FFI.Pointer(m * n * elementSize);
	var f_info = new FFI.Pointer(FORTRAN_INT);
	var f_lwork = new FFI.Pointer(FORTRAN_INT);
	var f_work;
	f_lwork.putInt32(-1);
	
	// get optimal size of workspace
	f_work = new FFI.Pointer(FORTRAN_INT);
	lapackFunc(f_m, f_n, f_a, f_lda, f_tau, f_work, f_lwork, f_info);
	lwork = f_work[op]();
	
	// allocate workspace
	f_work = new FFI.Pointer(lwork * elementSize);    
	f_lwork.putInt32(lwork);
	
	// perform QR decomp
	lapackFunc(f_m, f_n, f_a, f_lda, f_tau, f_work, f_lwork, f_info); 
	
	qr = {
	    R: fortranArray.fortranArrayToJSMatrix(f_a, m, n, elementSize),
	    tau: fortranArray.fortranArrayToJSArray(f_tau, Math.min(m, n), op, elementSize)
	};
	
	if(callback)
	    qr = callback(qr, m, n, f_m, f_n, f_a, f_lda, f_tau, f_work, f_lwork, f_info);
    });
    
    return qr;
}

function cloneMatrix(matrix, height, width) {
    var clone = [];

    height = height || matrix.length;
    width = width || matrix[0].length;

    for(var i = 0; i < height; i++) {
	var row = [];
	clone.push(row);

	for(var j = 0; j < width; j++) {
	    row.push(matrix[i][j]);
	}
    }

    return clone;
}

function swapRows(matrix, i, j) {
    var tmp = matrix[j];
    matrix[j] = matrix[i];
    matrix[i] = tmp;
    return matrix;
}


function cholesky(matrix) {
    var result = spotrf(matrix);
	return result;
}

function lu(matrix) {
    var result = sgetrf(matrix);
    var P = ipivToP(matrix.length, result.IPIV);
    var L = cloneMatrix(result.LU);
    var m = n = Math.min(matrix.length, matrix[0].length);

    for(var i = 0; i < L.length; i++) {
	for(var j = i; j < L[i].length; j++) {
	    if(i == j)
		L[i][j] = 1;
	    else
		L[i][j] = 0;
	}
    }

    return {
	L: L,
	U: zeroBottomLeft(cloneMatrix(result.LU, n, n)),
	P: P,
	IPIV: result.IPIV
    };
}

function ipivToP(m, ipiv){
    var P = eye(m);

    for(var i = 0; i < ipiv.length; i++) {
	if(i != ipiv[i] - 1)
            swapRows(P, i, ipiv[i] - 1);
    }

    return P;
}

function sgetrf(matrix) {
    return getrf(matrix, LAPACK.sgetrf_, FORTRAN_FLOAT);
}
function spotrf(matrix) {
    return potrf(matrix, LAPACK.spotrf_, FORTRAN_FLOAT);
}

function dgetrf(matrix) {
    return getrf(matrix, LAPACK.dgetrf_, FORTRAN_DOUBLE);
}

function potrf(matrix, lapackFunc, elementSize) {
    var result = {};
    matrixOp(matrix, elementSize, function(m, n, f_m, f_n, f_a, f_lda, f_letter_u, f_letter_l) {
        //var f_ipiv = new FFI.Pointer(Math.min(m, n) * FORTRAN_INT);
        var f_info = new FFI.Pointer(FORTRAN_INT);
        lapackFunc(f_letter_u, f_n, f_a, f_lda, f_info);
        result.cholesky = fortranArray.fortranArrayToJSMatrix(f_a, m, n, elementSize);
        result.LDA = f_lda.getInt32(); //fortranArray.fortranArrayToJSMatrix(f_a, m, n, elementSize);
	result.INFO = f_info.getInt32();
        //result.IPIV = fortranArray.fortranArrayToJSArray(f_ipiv, Math.min(m, n), 'getInt', elementSize);
    });

    return result;
}
function getrf(matrix, lapackFunc, elementSize) {
    var result = {};

    matrixOp(matrix, elementSize, function(m, n, f_m, f_n, f_a, f_lda) {
	var f_ipiv = new FFI.Pointer(Math.min(m, n) * FORTRAN_INT);
	var f_info = new FFI.Pointer(FORTRAN_INT);
	lapackFunc(f_m, f_n, f_a, f_m, f_ipiv, f_info);
	result.LU = fortranArray.fortranArrayToJSMatrix(f_a, m, n, elementSize);
	result.IPIV = fortranArray.fortranArrayToJSArray(f_ipiv, Math.min(m, n), 'getInt', elementSize);
    });

    return result;
}

function sgesvd(jobu, jobvt, matrix) {
    var f_jobu = new FFI.Pointer(FORTRAN_CHAR);
    var f_jobvt = new FFI.Pointer(FORTRAN_CHAR);
    f_jobu.putChar(jobu.charCodeAt(0));
    f_jobvt.putChar(jobvt.charCodeAt(0));
    var svd;

    matrixOp(matrix, FORTRAN_FLOAT, function(m, n, f_m, f_n, f_a, f_lda) {
	var f_s = new FFI.Pointer(Math.pow(Math.min(m, n), 2) * FORTRAN_FLOAT);
	var f_u = new FFI.Pointer(Math.pow(m, 2) * FORTRAN_FLOAT);
	var f_ldu = new FFI.Pointer(FORTRAN_INT);
	f_ldu.putInt32(m);

	// TODO: punting on dims for now. revisit with http://www.netlib.org/lapack/single/sgesvd.f
	var f_vt = new FFI.Pointer(Math.pow(n, 2) * FORTRAN_FLOAT);
	var f_ldvt = new FFI.Pointer(FORTRAN_INT);
	f_ldvt.putInt32(n);
	
	var lwork = -1;
	var f_work = new FFI.Pointer(FORTRAN_FLOAT);
	var f_lwork = new FFI.Pointer(FORTRAN_INT);
	f_lwork.putInt32(lwork);
	var f_info = new FFI.Pointer(FORTRAN_INT);

	LAPACK.sgesvd_(f_jobu, f_jobvt, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, 
		      f_work, f_lwork, f_info);

	lwork = f_work.getFloat();
	f_work = new FFI.Pointer(lwork * FORTRAN_FLOAT);
	f_lwork.putInt32(lwork);

	LAPACK.sgesvd_(f_jobu, f_jobvt, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, 
		      f_work, f_lwork, f_info);

	svd = {
	    U: fortranArray.fortranArrayToJSMatrix(f_u, m, m, FORTRAN_FLOAT),
	    S: fortranArray.fortranArrayToJSMatrix(f_s, n, n, FORTRAN_FLOAT),
	    VT: fortranArray.fortranArrayToJSMatrix(f_vt, n, n, FORTRAN_FLOAT)
	};
    });

    return svd;
}

exports.sgeqrf = sgeqrf;
exports.dgeqrf = dgeqrf;
exports.sgesvd = sgesvd;
exports.sgetrf = sgetrf;
exports.dgetrf = sgetrf;
exports.sgesv = sgesv;
exports.qr = qr;
exports.lu = lu;
exports.cholesky = cholesky;
