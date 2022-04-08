/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"

/*
 * Operation: C = A x B x B_t + A_t x A
 */

/*
 * Function to calculate the transpose
 */
double *transpose_matrix (int N, double *A) {
	int i, j;
	double *A_t;

	A_t = calloc(N * N, sizeof(double));
	if (!A_t) {
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < N; ++i) {
		/* column i from A_t */
		register double *A_t_ptr = A_t + i;
		/* line i from A */
		register double *A_ptr = A + i * N;

		for (j = 0; j < N; ++j) {
			*A_t_ptr = *A_ptr;
			A_t_ptr += N;
			A_ptr++;
		}
	}

	return A_t;
}

/*
 * Function to calculate A x B
 */
double* AB(int N, double *A, double *B) {
	int i, j, k;
	double *AB;
	
	AB = calloc(N * N, sizeof(double));
	if (!AB) {
		exit(EXIT_FAILURE);
	}

	/*
	 * A is upper triangular so it is enough to consider
	 * the columns from i to N
	 * 
	 * switch fors to i-k-j
	 * use registers
	 * 2x2 block calculation
	 * 
	 * operation: AB[i * N + j] += A[i * N + k] * B[k * N + j];
	 * 
	 * let c = AB, a = A, b = B for the following comments
	 */

	for (i = 0; i + 1 < N; i+=2) {
		/* ci */
		register double *ptr_AB_i = AB + i * N;
		/* ci+1 */
		register double *ptr_AB_i1 = AB + (i + 1) * N;
		/* ai */
	 	register double *ptr_A = A + i * (N + 1);
		for (k = i; k < N; k++) {
			/* ci */
			register double *ptr_AB_0 = ptr_AB_i;
			/* ci+1 */
			register double *ptr_AB_1 = ptr_AB_i1;
			/* bk */
			register double *ptr_B = B + k * N;
			/* aik */
			register double cons_value_i = *ptr_A;
			/* ai+1k */
			register double cons_value_i1 = *(ptr_A + N);
			for (j = 0; j + 1 < N; j += 2) {
				/* cij = aik * bkj */
				*ptr_AB_0 += cons_value_i * *ptr_B;
				/* cij+1 = aik + bkj+1 */
				*(ptr_AB_0 + 1) += cons_value_i * *(ptr_B + 1);
				/* ci+1j = ai+1k * bkj */
				*ptr_AB_1 += cons_value_i1 * *ptr_B;
				/* ci+1j+1 = ai+1k * bkj+1 */
				*(ptr_AB_1 + 1) += cons_value_i1 * *(ptr_B + 1);
				ptr_B += 2;
				ptr_AB_0 += 2;
				ptr_AB_1 += 2;
			}
			ptr_A++;
		}
	}

	return AB;
}

/*
 * Function to calculate (A x B) x B_t
 */
double* ABB_t(int N, double *AB, double *B) {
	int i, j, k;
	double *ABB_t;
	
	ABB_t = calloc(N * N, sizeof(double));
	if (!ABB_t) {
		exit(EXIT_FAILURE);
	}

	double *B_t = transpose_matrix(N, B);

	/*
	 * We do not know the form of either AB or B_t
	 * so all the fors go from 0 to N
	 * 
	 * compute the transpose
	 * switch fors to k-i-j
	 * use registers
	 * 2x2 block calculation
	 * 
	 * operation: ABB_t[i * N + j] += AB[i * N + k] * B_t[k * N + j];
	 * 
	 * let c = ABB_t, a = AB, b = B_t for the following comments
	 */

	for (i = 0; i < N; i += 2) {
		/* ai */
		register double *ptr_AB = AB + i * N;
		for (k = 0; k < N; k++) {
			/* ci */
			register double *ptr_ABBt_0 = ABB_t + i * N;
			/* ci+1 */
			register double *ptr_ABBt_1 = ABB_t + (i + 1) * N;
			/* bk0 */
			register double *ptr_B = B_t + k * N;
			/* aik */
			register double cons_value_i = *ptr_AB;
			/* ai+1k */
			register double cons_value_i1 = *(ptr_AB + N);
			for (j = 0; j + 1 < N; j += 2) {
				/* cij = aik * bkj */
				*ptr_ABBt_0 += cons_value_i * *ptr_B;
				/* ci+1j = ai+1k * bkj */
				*(ptr_ABBt_0 + 1) += cons_value_i * *(ptr_B + 1);
				/* cij+1 = aik * bkj+1 */
				*ptr_ABBt_1 += cons_value_i1 * *ptr_B;
				/* ci+1j+1 = ai+1k * bkj+1 */
				*(ptr_ABBt_1 + 1) += cons_value_i1 * *(ptr_B + 1);
				ptr_B += 2;
				ptr_ABBt_0 += 2;
				ptr_ABBt_1 += 2;
			}
			ptr_AB++;
		}
	}

	free(B_t);
	return ABB_t;
}

/*
 * Function to calculate A_t x A
 */
double* A_tA(int N, double *A_t, double *A) {
	int i, j, k;
	double *A_tA;
	
	A_tA = calloc(N * N, sizeof(double));
	if (!A_tA) {
		exit(EXIT_FAILURE);
	}

	/*
	 * A_t is lower triangular so it is enough to consider
	 * the columns from 0 to i + 1
	 *
	 * compute the transpose
	 * switch fors to i-k-j
	 * use registers
	 * 2x2 block calculation
	 * 
	 * operation: A_tA[i * N + j] += A_tt[i * N + k] * A[k * N + j];
	 * 
	 * let c = A_tA, a = A_tt, b = A for the following comments
	 */

	double *A_tt = transpose_matrix(N, A);

	for (i = 0; i < N; i += 2) {
		/* ai */
		register double *ptr_AtA = A_tt + i * N;
		for (k = 0; k < N; k++) {
			/* ci */
			register double *ptr_C_0 = A_tA + i * N;
			/* ci+1 */
			register double *ptr_C_1 = A_tA + (i + 1) * N;
			/* bk0 */
			register double *ptr_A = A + k * N;
			/* aik */
			register double cons_value_i = *ptr_AtA;
			/* ai+1k */
			register double cons_value_i1 = *(ptr_AtA + N);
			for (j = 0; j + 1 < N; j += 2) {
				/* cij = aik * bkj */
				*ptr_C_0 += cons_value_i * *ptr_A;
				/* ci+1j = ai+1k * bkj */
				*(ptr_C_0 + 1) += cons_value_i * *(ptr_A + 1);
				/* cij+1 = aik * bj+1 */
				*ptr_C_1 += cons_value_i1 * *ptr_A;
				/* ci+1j+1 = ai+1k * bkj+1 */
				*(ptr_C_1 + 1) += cons_value_i1 * *(ptr_A + 1);
				ptr_A += 2;
				ptr_C_0 += 2;
				ptr_C_1 += 2;
			}
			ptr_AtA++;
		}
	}

	free(A_tt);
	return A_tA;
}

/*
 * Function to compute the final sum
 */
double* sum(int N, double *A, double *B) {
	int i;
	double *C;
	
	C = calloc(N * N, sizeof(double));
	if (!C) {
		exit(EXIT_FAILURE);
	}

	/*
	 * Add the matrices
	 */
	register double *ptr_A = A;
	register double *ptr_B = B;
	register double *ptr_C = C;
	for (i = 0; i < N * N; i++) {
		*ptr_C = *ptr_A + *ptr_B;
		ptr_A++;
		ptr_B++;
		ptr_C++;
	}

	return C;
}


double* my_solver(int N, double *A, double* B) {
	printf("OPT SOLVER\n");
	register double *ptr_AB = AB(N, A, B);
	register double *ptr_ABB_t = ABB_t(N, ptr_AB, B);
	register double *ptr_A_tA = A_tA(N, A, A);
	double *ptr_sum = sum(N, ptr_ABB_t, ptr_A_tA);

	free(ptr_AB);
	free(ptr_ABB_t);
	free(ptr_A_tA);
	
	return ptr_sum;
}
