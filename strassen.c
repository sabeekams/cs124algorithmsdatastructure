#include <stdio.h>
#include <time.h>
#include "strassen.h"
#include <stdlib.h>
#define THRESHOLD 4
#define STRASSEN 8 

/* FLAGS:
    0 : Without padding 
	1: Padding Test 
*/


int main(int argc, char* argv[]) {


	//get the dimension
	int dimension = (int) strtol(argv[2], NULL, 10);
	
	int thresh = THRESHOLD;
    int flag = strtol(argv[1], NULL, 10);
	//get the padding 
	
	int flagcheck = matrixFlagCheck(dimension, thresh, flag);
	

	
	matrix* buildmtx[2]; 

	buildmtx[0] = buildMatrix(flagcheck);
	buildmtx[1] = buildMatrix(flagcheck);

	//fill up the file into the matrixes
	FILE* file = fopen(argv[3], "r");
	
	char buffer[20];
	int n = dimension * dimension;
	
	for (int i = 0; i < 2 * n; i++) {
		fscanf(file, "%s\n", buffer);
		buildmtx[i / n]->array[buildmtx[i / n]->firstrow + (i % n) / dimension][buildmtx[i / n]->firstcol + i % dimension] = (int) strtol(buffer, NULL, 10);
	}

	fclose(file);


	matrix* product = buildMatrix(flagcheck);
	time_t start = clock();
	multiplyMatrix(product, buildmtx[0], buildmtx[1]);
	time_t end = clock();


		for (int i = 0; i < dimension; i++) {
			printf("%i\n", getNumberfromMatrix(product, i, i));
		}
		
		printf("Threshold %i dimension is %i: run time is  %ld\n", thresh, dimension, end - start);
}


int matrixFlagCheck(int dimension, int threshold, int flag) {
	if (flag == 0 || flag == 1) {
		return dimension;
	}
}


// Gets the element from matrix mtx at row i and column j
int getNumberfromMatrix(matrix* mtx, int i, int j) {
	return mtx->array[mtx->firstrow + i][mtx->firstcol + j];
}


// Creates a matrix with side length d
matrix* buildMatrix(int d) {
	matrix* newmatrix = malloc(sizeof(matrix));
	newmatrix->firstrow = 0;
	newmatrix->firstcol = 0;

	newmatrix->dm = d;
	newmatrix->array = malloc(sizeof(int*) * d);
	for (int i = 0; i < d; i++) {
		newmatrix->array[i] = calloc(d, sizeof(int));
	}
	return newmatrix;
}


void multiplyMatrix(matrix* p, matrix* matrix1, matrix* matrix2) {

	if (matrix1->dm <= THRESHOLD){
		/*Standard Multiplication*/
		
		for (int i = 0; i < matrix1->dm; i++) {
		for (int j = 0; j < matrix1->dm; j++) {
			for (int k = 0; k < matrix2->dm; k++) {
			    int a = getNumberfromMatrix(p, i, k) + getNumberfromMatrix(matrix1, i, j) * getNumberfromMatrix(matrix2, j, k);
			    
			    p->array[p->firstrow + i][p->firstcol + k] = a;  

			}
		}
	}	
	    return;
}

    matrix* matrices2[4];
	matrix* matrices1[4]; 

	matrixDivision(matrix2, matrices2);
	matrixDivision(matrix1, matrices1);
	

	// Pieces for Strassen's (index 0 left blank) and temporary matrices from addition
	matrix* product[STRASSEN]; //array of matrices
	for (int i = 1; i < STRASSEN; i++) {
		product[i] = buildMatrix(p->dm / 2); 
	}
	
	int halvedmatrix = p->dm / 2;
	matrix* tmpmatrix1 = buildMatrix(halvedmatrix);
	matrix* tmpmatrix2 = buildMatrix(halvedmatrix);

	// Calculating pieces	
	int f = matrixSubstract(tmpmatrix1, matrices1[0], matrices1[2]);
	int g = matrixAdd(tmpmatrix2, matrices2[0], matrices2[1]);
	
	multiplyMatrix(product[7], f, g);
	multiplyMatrix(product[6], matrixSubstract(tmpmatrix1, matrices1[1], matrices1[3]), matrixAdd(tmpmatrix2, matrices2[2], matrices2[3]));
	
	int a = matrixAdd(tmpmatrix2, matrices2[0], matrices2[3]);
	
	multiplyMatrix(product[5], matrixAdd(tmpmatrix1, matrices1[0], matrices1[3]), a);
	
	int b = matrixSubstract(tmpmatrix1, matrices2[2], matrices2[0]);
	
	multiplyMatrix(product[4], matrices1[3], b); 
	
	int c = matrixAdd(tmpmatrix1, matrices1[2], matrices1[3]);
	
	multiplyMatrix(product[3], c, matrices2[0]);
	
	int d = matrixAdd(tmpmatrix1, matrices1[0], matrices1[1]);
	
	multiplyMatrix(product[2], d , matrices2[3]);
	
	int e = matrixSubstract(tmpmatrix1, matrices2[1], matrices2[3]);
	
	multiplyMatrix(product[1], matrices1[0], e);
	

	int newdim = p->dm / 2;
	
	p->dm = newdim;
	
	matrixAdd(p, product[5], matrixAdd(tmpmatrix1, product[4], matrixSubstract(tmpmatrix2, product[6], product[2])));
	
	p->firstcol = p->dm;
	matrixAdd(p, product[1], product[2]);
	p->firstrow = p->dm;
	
	matrixAdd(p, product[5], matrixSubstract(tmpmatrix1, product[1], matrixAdd(tmpmatrix2, product[3], product[7])));
	p->firstcol = 0;
	
	matrixAdd(p, product[3], product[4]);
	p->firstrow = 0;
	
	p->dm = p->dm * 2;
}



void matrixDivision(matrix* matrixdivide, matrix** matrices) {  //Slicing up the matrixes for strassen calculation 
	
	for (int i = 0; i < 4; i++){
		matrices[i] = malloc(sizeof(matrix));
		matrices[i]->dm = matrixdivide->dm / 2;
		matrices[i]->array = matrixdivide->array;


		if (i == 1 || i == 3){
			
			matrices[i]->firstcol = matrixdivide->firstcol + matrixdivide->dm / 2;
		} else {
            matrices[i]->firstcol = matrixdivide->firstcol;
		}
		

		if (i == 2 || i == 4){
		   
            matrices[i]->firstrow = matrixdivide->firstrow + matrixdivide->dm / 2;
		} else {
		    matrices[i]->firstrow = matrixdivide->firstrow;
		}
	}
}

/*Adds the two matrices together*/
matrix* matrixAdd(matrix* s, matrix* matrixone, matrix* matrixtwo) {
	for (int i = 0; i < matrixtwo->dm; i++) {
		for (int j = 0; j < matrixtwo->dm; j++) {
		
		    int a = matrixone->array[matrixone->firstrow + i][matrixone->firstcol + j] + matrixtwo->array[matrixtwo->firstrow + i][matrixtwo->firstcol + j];

			s->array[s->firstrow + i][s->firstcol + j] = a;
		}
	}
	return s;
}



/*Substracts the two matrices*/

matrix* matrixSubstract(matrix* prod, matrix* matrixone, matrix* matrixtwo) {
	
	for (int i = 0; i < matrixtwo->dm; i++) {
		for (int j = 0; j < matrixtwo->dm; j++) {
		
		int b = matrixone->array[matrixone->firstrow + i][matrixone->firstcol + j] - matrixtwo->array[matrixtwo->firstrow + i][matrixtwo->firstcol + j];
		prod->array[prod->firstrow + i][prod->firstcol + j] = b;
		}
	}

	return prod;
}



