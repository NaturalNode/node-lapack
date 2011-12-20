
var fortranArray = require('./fortranArray');
var FFI = require('node-ffi');

var LAPACK = new FFI.Library('liblapack', {
    "SGEQRF": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]]
});

var FORTRAN_INT = 4;

function sgeqrf(matrix) {
    var m = matrix.length;
    var n = matrix[0].length;
    var f_m = new FFI.Pointer(FORTRAN_INT);
    var f_n = new FFI.Pointer(FORTRAN_INT);
    var f_lda = new FFI.Pointer(FORTRAN_INT);
    var f_tau = new FFI.Pointer(m * n);
    var f_info = new FFI.Pointer(FORTRAN_INT);
    var f_lwork = new FFI.Pointer(FORTRAN_INT);
    var f_work;
    f_m.putInt32(m);
    f_n.putInt32(n);
    f_lda.putInt32(Math.max(1, m));
    f_lwork.putInt32(-1);
    var f_a = fortranArray.jsToFortranArray(matrix);

    // get optimal size of workspace
    f_work = new FFI.Pointer(FORTRAN_INT);
    LAPACK.SGEQRF(f_m, f_n, f_a, f_lda, f_tau, f_work, f_lwork, f_info);
    lwork = f_work.getFloat();    

    // allocate workspace
    f_work = new FFI.Pointer(lwork);    
    f_lwork.putInt32(lwork);

    // perform QR decomp
    LAPACK.SGEQRF(f_m, f_n, f_a, f_lda, f_tau, f_work, f_lwork, f_info);    

    var r = fortranArray.fortranArrayToJS(f_a, m, n);
    var q = fortranArray.fortranArrayToJS(f_tau, m, n);

    return {Q: q, R: r};
}

exports.sgeqrf = sgeqrf;
