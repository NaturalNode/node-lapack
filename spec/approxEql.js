
var precision =  1e-6;

function approxEql(matrixA, matrixB) {
    var i = matrixA.length;
 
    while(i--) {
	var j = matrixA[0].length;
	
	while(j--) {
            if((Math.abs(matrixA[i][j] - matrixB[i][j]) > precision)) {
		return false;
	    }
        }
    }

    return true;
}

module.exports = approxEql;