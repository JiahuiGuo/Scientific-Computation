#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

int main()
{
    int nRow, nCol, nnZero;
    int *aijRow, *aijCol;
    double *value;
    int i;

    // Open the file
    FILE* pf;
    pf = fopen("matrix.output", "r");
    if (pf == NULL) {
        fprintf(stderr,  "Can't open input file!\n");
        return EXIT_FAILURE;
    }
    // Read in the matrix in AIJ Format
    fscanf(pf, "%d %d %d", &nRow, &nCol, &nnZero);
    
    // Allocate memory for AIJ format
    aijRow = (int*)malloc(sizeof(int)*nnZero);
    aijCol = (int*)malloc(sizeof(int)*nnZero);
    value = (double*)malloc(sizeof(double)*nnZero);
    
    // Read in the AIJ format matrix
    for(i = 0; i < nnZero; ++i){
        // Note: The index stored is not C style!
        fscanf(pf, "%d %d %lf", &aijRow[i], &aijCol[i], &value[i]);   
        // Change to C style
        aijRow[i] -= 1;
        aijCol[i] -= 1;
    }
    fclose(pf);

    // Reorder the old index
    // i_new = (i_old - 1) / 8660 + 3 * ((i_old - 1) % 8660) + 1
    for(i = 0; i < nnZero; ++i)
    {
        // i_old - 1 is performed previously
        // Reserver C style here
        aijRow[i] = (aijRow[i] / 8660) + 3 * (aijRow[i] % 8660); 
        aijCol[i] = (aijCol[i] / 8660) + 3 * (aijCol[i] % 8660); 
    }

    pf = fopen("matrix.reorder", "w");
    fprintf(pf, "%d %d %d\n", nRow, nCol, nnZero);
    for(i = 0; i < nnZero; ++i)
        fprintf(pf, "%d %d %lf\n", aijRow[i] + 1, aijCol[i] + 1, value[i]);
    fclose(pf);

	char *sort[] = {"sort.sh", NULL};
    execv("./sort.sh", sort);
	
	/*
    // Execute a bash script to order the matrix.reorder file
    char *sort_ccs[] = { "/bin/bash", "-c", "head -1 matrix.reorder > matrix.reorder.ccs | tail -n +2 matrix.reorder | sort -n -k 2,2 -k 1,1 >> matrix.reorder.ccs", NULL};
    execvp(sort_ccs[0], sort_ccs);

    char *sort_crs[] = { "/bin/bash", "-c", "head -1 matrix.reorder > matrix.reorder.crs | tail -n +2 matrix.reorder | sort -n -k 1,1 -k 2,2 >> matrix.reorder.crs", NULL};
    execvp(sort_crs[0], sort_crs);
	*/
	return EXIT_SUCCESS;
}
