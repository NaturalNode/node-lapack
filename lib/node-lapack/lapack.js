
var fortranArray = require('./fortranArray');
var FFI = require('node-ffi');

var LAPACK;

try {
    LAPACK = new FFI.Library('liblapack', {
	"SGEQRF": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", 
			    "pointer", "pointer", "pointer"]],
	"SGESVD": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", 
			    "pointer", "pointer", "pointer", "pointer", "pointer", 
			    "pointer", "pointer", "pointer", "pointer", ]]
	
    });
} catch(e) {
    console.log("!!! node-lapack requires the native lapack to be built as a shared lib.");
    console.log(e);
}

var FORTRAN_INT = 4;
var FORTRAN_CHAR = 1;
var FORTRAN_FLOAT = 4;

function matrixOp(matrix, callback) {
    var m = matrix.length;
    var n = matrix[0].length;    
    var f_m = new FFI.Pointer(FORTRAN_INT);
    var f_n = new FFI.Pointer(FORTRAN_INT);
    var f_a = fortranArray.jsToFortranArray(matrix);
    var f_lda = new FFI.Pointer(FORTRAN_INT);

    f_m.putInt32(m);
    f_n.putInt32(n);
    f_lda.putInt32(Math.max(1, m));
    
    callback(m, n, f_m, f_n, f_a, f_lda);
}

function sgeqrf(matrix) {
    var qr;

    matrixOp(matrix, function(m, n, f_m, f_n, f_a, f_lda) {
	var f_tau = new FFI.Pointer(m * n * FORTRAN_FLOAT);
	var f_info = new FFI.Pointer(FORTRAN_INT);
	var f_lwork = new FFI.Pointer(FORTRAN_INT);
	var f_work;
	f_lwork.putInt32(-1);
	
	// get optimal size of workspace
	f_work = new FFI.Pointer(FORTRAN_INT);
	LAPACK.SGEQRF(f_m, f_n, f_a, f_lda, f_tau, f_work, f_lwork, f_info);
	lwork = f_work.getFloat();
	
	// allocate workspace
	f_work = new FFI.Pointer(lwork * FORTRAN_FLOAT);    
	f_lwork.putInt32(lwork);
	
	// perform QR decomp
	LAPACK.SGEQRF(f_m, f_n, f_a, f_lda, f_tau, f_work, f_lwork, f_info);    
	
	qr = {
	    R: fortranArray.fortranArrayToJS(f_a, m, n), 
	    Q: fortranArray.fortranArrayToJS(f_tau, m, n)
	};	
    });

    return qr;
}

function sgesvd(jobu, jobvt, matrix) {
    var f_jobu = new FFI.Pointer(FORTRAN_CHAR);
    var f_jobvt = new FFI.Pointer(FORTRAN_CHAR);
    f_jobu.putChar(jobu.charCodeAt(0));
    f_jobvt.putChar(jobvt.charCodeAt(0));
    var svd;

    matrixOp(matrix, function(m, n, f_m, f_n, f_a, f_lda) {
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

	LAPACK.SGESVD(f_jobu, f_jobvt, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, 
		      f_work, f_lwork, f_info);

	lwork = f_work.getFloat();
	f_work = new FFI.Pointer(lwork * FORTRAN_FLOAT);
	f_lwork.putInt32(lwork);

	LAPACK.SGESVD(f_jobu, f_jobvt, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, 
		      f_work, f_lwork, f_info);

	svd = {
	    U: fortranArray.fortranArrayToJS(f_u, m, m),
	    S: fortranArray.fortranArrayToJS(f_s, n, n),
	    VT: fortranArray.fortranArrayToJS(f_vt, n, n)
	};
    });

    return svd;
}

exports.sgeqrf = sgeqrf;
exports.sgesvd = sgesvd;
