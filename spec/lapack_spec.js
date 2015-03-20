
var lapack = require('../lib/node-lapack');
var approxEql = require('./approxEql');
var expect = require('expect');

describe('lapack', function() {
    var M = [
        [2, 1, 1],
        [1, 1, 1],
        [1, 1, 3]
    ];

    var luIn = [
	[11, 3, 11, 3],
	[ 4, 2,  1, 4],
	[-4, 5,  3, 1],
	[-9, 4,  3, 9]
    ];

    var lu2 = [
        [3, 6],
        [2, 3],
        [4, 3],
        [2, 120],
    ];
    
    it('should solve', function() {
	var A = [
	    [2, 4],
	    [2, 8]
	];

	var B = [[2], [4]];
	expect(lapack.sgesv(A, B).X).toEqual([[0], [.5]]);	
	expect(lapack.sgesv(A, B).P).toEqual([[1, 0], [0, 1]]);
    });

    it('should lu rectangular matrix', function() {
	var lu = lapack.lu(lu2);
	expect(lu.IPIV).toEqual([3, 4]);
	expect(lu.L).toEqual([ [ 1, 0 ],
			       [ 0.5, 1 ],
			       [ 0.75, 0.03164556622505188 ],
			       [ 0.5, 0.012658227235078812 ] ]);
	expect(lu.U).toEqual([[4, 3], [0, 118.5]]);
    });

    it('should lu', function() {
	var lu = lapack.lu(luIn);
	expect(approxEql(lu.L, ([
	    [ 1, 0, 0, 0 ],
	    [ -0.8181818723678589, 1, 0, 0 ],
	    [ 0.3636363744735718, 0.14084506034851074, 1, 0 ],
	    [ -0.3636363744735718, 0.9436619877815247,
	      0.9219222664833069, 1 ]	    
	]))).toBe(true);
	
	expect(approxEql(lu.U, [
	    [ 11, 3, 11, 3 ],
	    [ 0, 6.454545497894287,
	      12.000000953674316, 11.454545974731445 ],
	    [ 0, 0, -4.690140724182129, 1.2957748174667358 ],
	    [ 0, 0, 0, -9.912914276123047 ]	    
	])).toBe(true);
	
	expect(approxEql(lu.P, [
	    [ 1, 0, 0, 0 ],
	    [ 0, 0, 0, 1 ],
	    [ 0, 1, 0, 0 ],
	    [ 0, 0, 1, 0 ]	    
	])).toBe(true);
    });
    
    it('should sgetrf', function() {
	var result = lapack.sgetrf(luIn);
	expect(approxEql(result.LU, [ [ 11, 3, 11, 3 ],
				      [ -0.8181818723678589,
					6.454545497894287,
					12.000000953674316,
					11.454545974731445 ],
				      [ 0.3636363744735718,
					0.14084506034851074,
					-4.690140724182129,
					1.2957748174667358 ],
				      [ -0.3636363744735718,
					0.9436619877815247,
					0.9219222664833069,
					-9.912914276123047 ] ])).toBe(true);
	expect(approxEql(result.IPIV, [ 3, 2, 3, 4 ])).toBe(true);
    });

    if('should sgetrf and dgetrf approx eql', function() {
	expect(approxEql(lapack.sgetrf(luIn).LU, lapack.dgetrf(luIn).LU)).toBe(true);
    });

    it('shoud dgeqrf and sgeqrf approximately equal', function() {
        expect(approxEql(lapack.dgeqrf(M).R,
                         lapack.sgeqrf(M).R)).toBe(true);
    });
    
    it('shoud sgeqrf', function() {
	var qr = lapack.sgeqrf(M);

	expect(approxEql(qr.R, 
			 [ [ -2.4494898,
			     -1.6329932,
			     -2.4494896 ],
			   [  0.2247450,
			     -0.57735044,
			     -1.7320511 ],
			   [ 0.2247450,
			     0.41421360,
			     1.4142138 ] ])).toBe(true);
    });

    
    it('shoud dgeqrf', function() {
	var qr = lapack.dgeqrf(M);
	expect(approxEql(qr.R, [[-2.4494898319244385,
				 -1.632993221282959,
				 -2.4494895935058594 ],
			 [ 0.224745035171509,
			  -0.5773504376411438,
			  -1.732051134109497 ],
			 [ 0.224745035171509,
			   0.4142135977745056,
			   1.41421377658844 ]])).toBe(true);
    });

    it('should qr', function() {
	var qr = lapack.qr(M);
	expect(approxEql(qr.R, [ [ -2.4494898,
			     -1.6329932,
			     -2.4494896 ],
			   [ 0.0,
			     -0.57735044,
			     -1.73205113 ],
			   [ 0.0,
			     0.0,
			     1.41421377 ] ]
	)).toBe(true);
	expect(approxEql(qr.Q, [
	    [-0.81649661, 0.57735032, 0.0],
	    [-0.40824828, -0.57735050, -0.7071067],
	    [-0.40824828, -0.57735044, 0.70710683]
	])).toBe(true);
    });

    it('should sgesvd', function() {
        var svd = lapack.sgesvd('A', 'A', M);
	
	expect(approxEql(svd.S, [ [ 4.214320182800293, 0, 2.462372871899283e-38 ],
				  [ 1.4608111381530762,
				    2.4612226861197653e-38,
				    1.401298464324817e-45 ],
				  [ 0.3248690962791443,
				    1.401298464324817e-45,
				    2.462375113976826e-38 ] 
				]));

	expect(approxEql(svd.U, [ [ -0.5206573009490967,
				    0.7392387390136719,
				    0.4271322190761566 ],
				  [ -0.39711257815361023,
				    0.23319198191165924,
				    -0.8876503705978394 ],
				  [ -0.7557893991470337,
				    -0.631781280040741,
				    0.1721479445695877 ] ]));

	expect(approxEql(svd.VT, [ [ -0.5206573605537415,
				     -0.3971126079559326,
				     -0.7557893991470337 ],
				   [ 0.7392386794090271,
				     0.23319187760353088,
				     -0.631781280040741 ],
				   [ 0.42713218927383423,
				     -0.8876502513885498,
				     0.1721479296684265 ] ]));
    });
});

