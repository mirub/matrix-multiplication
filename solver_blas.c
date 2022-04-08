/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include <string.h>

#include "utils.h"
#include "cblas.h"

/*
 * C = A x B x B_t + A_t x A
 */

/*
 * Function to calculate A x B
 */
double* AB(int N, double *A, double *B) {
	double *AB;
	
	AB = calloc(N * N, sizeof(double));
	if (!AB) {
		exit(EXIT_FAILURE);
	}

	/*
	 * Needed to create a copy so it can be used
	 * and overwritten by cblas_dtrmm 
	 */
	memcpy(AB, B, N * N * sizeof(double));

	cblas_dtrmm(
		CblasRowMajor,
		CblasLeft,
		CblasUpper,
		CblasNoTrans,
		CblasNonUnit,
		N, N, 1.0, A, N,
		AB, N);

	return AB;
}

/*
 * Function to calculate A_t x A
 */
double* A_tA(int N, double *A_t, double *A) {
	double *A_tA;
	
	A_tA = calloc(N * N, sizeof(double));
	if (!A_tA) {
		exit(EXIT_FAILURE);
	}

	/*
	 * Needed to create a copy so it can be used
	 * and overwritten by cblas_dtrmm 
	 */
	memcpy(A_tA, A, N * N * sizeof(double));

	cblas_dtrmm(
		CblasRowMajor,
		CblasLeft,
		CblasUpper,
		CblasTrans,
		CblasNonUnit,
		N, N, 1.0, A, N,
		A_tA, N);

	return A_tA;
}

/*
 * Function to calculate (A x B) x B_t + A_t * A
 */
double* ABB_t(int N, double *AB, double *B, double *A_tA) {
	double *sum;
	
	sum = calloc(N * N, sizeof(double));
	if (!sum) {
		exit(EXIT_FAILURE);
	}

	/*
	 * Needed to create a copy so it can be used
	 * and overwritten by cblas_dgemm
	 */
	memcpy(sum, A_tA, N * N * sizeof(double));

	cblas_dgemm(
		CblasRowMajor,
		CblasNoTrans,
		CblasTrans,
		N, N, N, 1.0, AB,
		N, B, N, 1.0, sum, N);

	return sum;
}

double* my_solver(int N, double *A, double *B) {
	printf("BLAS SOLVER\n");
	double *ptr_AB = AB(N, A, B);
	double *ptr_A_tA = A_tA(N, A, A);
	double *ptr_ABB_t = ABB_t(N, ptr_AB, B, ptr_A_tA);

	free(ptr_AB);
	free(ptr_A_tA);

	return ptr_ABB_t;
}
