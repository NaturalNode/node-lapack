
node-lapack
===========

A node.js wrapper for the high-performance LAPACK linear algebra library.

Prerequisites
=============

This library require LAPACK to be built and installed as a shared library.
In time the entire build process may be unified into this project, but that's
some time away.

Usage
=====

    var lapack = require('lapack');

    /*
    LAPACK functions
    */

    var result = lapack.sgeqrf([
        [1, 2, 3],
        [3, 4, 5],
        [5, 6, 7]
    ]);

    console.log(result.R);
    console.log(result.tau);

    result = sgesvd('A', 'A', [
        [1, 2, 3],
        [3, 4, 5],
        [5, 6, 7]
    ]);

    console.log(result.U);
    console.log(result.S);
    console.log(result.VT);

    result = lapack.sgetrf([
        [1, 2, 3],
        [3, 4, 5],
        [5, 6, 7]
    ]);

    console.log(result.LU);
    console.log(result.IPIV);

    /*
    conveniently wrapped processes
    */

    // perform a complete LU factorization
    var lu = lapack.lu([
        [1, 2, 3],
        [3, 4, 5],
        [5, 6, 7]
    ]);

    console.log(lu.L);
    console.log(lu.U);
    console.log(lu.P);


    // perform a complete QR factorization
    var qr = lapack.qr([
        [1, 2, 3],
        [3, 4, 5],
        [5, 6, 7]
    ]);

    console.log(qr.Q);
    console.log(qr.R);
