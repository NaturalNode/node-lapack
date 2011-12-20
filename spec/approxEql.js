
var precision =  1e-6;

function approxEql(matrixA, matrixB) {
    var i = matrixA.length;
 
    while(i--) {
	var j = matrixA[0].length;
	
	while(j--) {
            if((Math.abs(matrixA[i][j] - matrixB[i][j]) > precision)) {
		console.log(matrixA[i][j]);
		console.log(matrixB[i][j]);
		console.log(Math.abs(matrixA[i][j] - matrixB[i][j]));

		return false;
	    }
        }
    }

    return true;
}

module.exports = approxEql;