
typedef struct matrix {
 	int firstrow;		
 	int firstcol;	
    int** array;
 	int dm;			
} matrix; 


int getNumberfromMatrix(matrix* mtx, int i, int j);
int main(int argc, char* argv[]);
matrix* buildMatrix(int d);
void matrixDivision(matrix* mtx, matrix** matrices);

matrix* matrixAdd(matrix* s, matrix* m1, matrix* m2);
matrix* matrixSubstract(matrix* s, matrix* m1, matrix* m2);

void printMatrix(matrix* mtx);

void multiplyMatrix(matrix* p, matrix* m1, matrix* m2);
int checkFlag (int flag, int dimension);
