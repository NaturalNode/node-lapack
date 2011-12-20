
var lapack = require('lib/node-lapack');
var approxEql = require('./approxEql');

describe('lapack', function() {
    it('shoud qr', function() {
	var qr = lapack.sgeqrf([
	    [2, 1, 1],
	    [1, 1, 1],
	    [1, 1, 3]
	]);
	
	expect(approxEql(qr.R, 
			 [ [ -2.4494898319244385,
			     -1.632993221282959,
			     -2.4494895935058594 ],
			   [ 0.22474488615989685,
			     -0.5773501992225647,
			     -1.732050895690918 ],
			   [ 0.22474488615989685,
			     0.41421353816986084,
			     1.4142135381698608 ] ])).toBeTruthy();
    });

    it('should svd', function() {
        var svd = lapack.sgesvd('A', 'A', [
            [2, 1, 1],
            [1, 1, 1],
            [1, 1, 3]
        ]);
    });
});

