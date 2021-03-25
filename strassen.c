#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

size_t n0 = 43;

struct matrix {
	int* array;
	size_t rows;
	size_t cols;

	size_t dimension;
	size_t start_row;
	size_t start_col;
};
const struct matrix NULL_MATRIX = { NULL };

inline size_t min(size_t a, size_t b) {
	return (a < b) ? a : b;
}

inline size_t max(size_t a, size_t b) {
	return (a > b) ? a : b;
}

inline struct matrix quadrant(struct matrix origin, bool i, bool j) {
	size_t newdim = origin.dimension>>1;
	size_t
		start_row = origin.start_row + i*newdim,
		start_col = origin.start_col + j*newdim;

	if (start_row >= origin.rows || start_col >= origin.cols) {
		return (struct matrix) {
			NULL,
			.dimension = newdim
		};
	}

	return (struct matrix) {
		origin.array,
		origin.rows,
		origin.cols,
		newdim,
		start_row,
		start_col
	};
}

inline int get(struct matrix matrix, size_t row, size_t col) {
	assert(row < matrix.dimension && col < matrix.dimension);
	if (!matrix.array) {
		return 0;
	}

	row += matrix.start_row;
	col += matrix.start_col;
	if (row >= matrix.rows || col >= matrix.cols) {
		return 0;
	}

	return matrix.array[row*matrix.cols + col];
}

// N.B. Caller is responsible for freeing iff C's array is null and returned matrix's array is not A.array, B.array, or null.
struct matrix add(struct matrix A,
                  int multiplier,
                  struct matrix B,
                  struct matrix C) {
	assert(A.dimension == B.dimension);

	if (!C.array) {
		if (!B.array) {
			return A;
		} else if (!A.array && multiplier == 1) {
			return B;
		}

		size_t
			rows = min(max(A.rows - A.start_row, B.rows - B.start_row), A.dimension),
			cols = min(max(A.cols - A.start_col, B.cols - B.start_col), A.dimension);
		C = (struct matrix) {
			malloc(rows*cols*sizeof(int)),
			rows,
			cols,
			A.dimension,
			0,
			0
		};
	}

	for (size_t i = 0; i < C.rows - C.start_row && i < C.dimension; i++) {
		for (size_t j = 0; j < C.cols - C.start_col && j < C.dimension; j++) {
			C.array[(C.start_row + i)*C.cols + C.start_col + j] = get(A, i, j) + multiplier*get(B, i, j);
		}
	}

	return C;
}

// N.B. Caller is responsible for freeing iff returned matrix's array is not null.
struct matrix naive(struct matrix A,
                    struct matrix B) {
	assert(A.dimension == B.dimension);

	if (!A.array) {
		return A;
	} else if (!B.array) {
		return B;
	}

	size_t
		rows = min(A.rows, A.dimension),
		cols = min(B.cols, B.dimension);
	struct matrix C = {
		calloc(rows*cols, sizeof(int)),
		rows,
		cols,
		A.dimension,
		0,
		0
	};

	for (size_t i = 0; i < C.rows; i++) {
		for (size_t j = 0; j < C.cols; j++) {
			for (size_t k = 0; k < min(min(A.cols, B.rows), A.dimension); k++) {
				C.array[i*C.cols + j] += get(A, i, k)*get(B, k, j);
			}
		}
	}

	return C;
}

void free_matrix(struct matrix matrix, struct matrix origin) {
	if (matrix.array && matrix.array != origin.array) {
		free(matrix.array);
	}
}

// N.B. Caller is responsible for freeing iff returned matrix's array is not null.
struct matrix strassen(struct matrix A,
                       struct matrix B) {
	assert(A.dimension == B.dimension);

	if (!A.array) {
		return A;
	} else if (!B.array) {
		return B;
	}

	if (A.dimension <= n0) {
		return naive(A, B);
	}

	struct matrix
		a = quadrant(A, 0, 0),
		b = quadrant(A, 0, 1),
		c = quadrant(A, 1, 0),
		d = quadrant(A, 1, 1),
		e = quadrant(B, 0, 0),
		f = quadrant(B, 0, 1),
		g = quadrant(B, 1, 0),
		h = quadrant(B, 1, 1);

	struct matrix
		f_minus_h = add(f, -1, h, NULL_MATRIX),
		a_plus_b  = add(a,  1, b, NULL_MATRIX),
		c_plus_d  = add(c,  1, d, NULL_MATRIX),
		g_minus_e = add(g, -1, e, NULL_MATRIX),
		a_plus_d  = add(a,  1, d, NULL_MATRIX),
		e_plus_h  = add(e,  1, h, NULL_MATRIX),
		b_minus_d = add(b, -1, d, NULL_MATRIX),
		g_plus_h  = add(g,  1, h, NULL_MATRIX),
		a_minus_c = add(a, -1, c, NULL_MATRIX),
		e_plus_f  = add(e,  1, f, NULL_MATRIX);

	struct matrix P[] = {
		[1] = strassen(a, f_minus_h),
		[2] = strassen(a_plus_b, h),
		[3] = strassen(c_plus_d, e),
		[4] = strassen(d, g_minus_e),
		[5] = strassen(a_plus_d, e_plus_h),
		[6] = strassen(b_minus_d, g_plus_h),
		[7] = strassen(a_minus_c, e_plus_f)
	};

	free_matrix(a_minus_c, A);
	free_matrix(a_plus_b,  A);
	free_matrix(a_plus_d,  A);
	free_matrix(b_minus_d, A);
	free_matrix(c_plus_d,  A);
	free_matrix(e_plus_f,  B);
	free_matrix(e_plus_h,  B);
	free_matrix(f_minus_h, B);
	free_matrix(g_minus_e, B);
	free_matrix(g_plus_h,  B);

	size_t
		rows = min(A.rows, A.dimension),
		cols = min(B.cols, B.dimension);
	struct matrix C = {
		malloc(rows*cols*sizeof(int)),
		rows,
		cols,
		A.dimension>>1,
		0,
		0
	};

	add(P[5], 1, P[4], C);
	add(C, -1, P[2], C);
	add(C, 1, P[6], C);

	C.start_col = C.dimension;
	if (C.start_col < C.cols) {
		add(P[1], 1, P[2], C);
	}

	C.start_row = C.dimension;
	C.start_col = 0;
	if (C.start_row < C.rows) {
		add(P[3], 1, P[4], C);

		C.start_col = C.dimension;
		if (C.start_col < C.cols) {
			add(P[5], 1, P[1], C);
			add(C, -1, P[3], C);
			add(C, -1, P[7], C);
		}
	}

	for (int i = 1; i <= 7; i++) {
		free_matrix(P[i], NULL_MATRIX);
	}

	C.dimension = A.dimension;
	C.start_row = 0;
	C.start_col = 0;
	return C;
}

int main(int argc, char** argv) {
	// Make sure every line is flushed immediately:
	setvbuf(stdout, NULL, _IOLBF, 0);

	if (argc < 4) {
		fprintf(stderr, "usage: %s option dimension inputfile\n", (argc > 0) ? argv[0] : "./strassen");
		fprintf(stderr, "\tif option == 0, read 2 dimension*dimension matrices A and B from inputfile and print the diagonal of AB\n");
		fprintf(stderr, "\tif option == -1, find the crossover point n0\n");
		fprintf(stderr, "\tif option > 0, count triangles in a graph of 1024 vertices with each edge included with probability option/100.0\n");
		return 1;
	}

	int option = atoi(argv[1]);
	int dimension = atoi(argv[2]);
	char* inputfile = argv[3];

	// Round up so Strassen's always works on an even dimension:
	int padded = n0;
	if (padded > dimension) {
		padded = dimension;
	} else {
		while (padded < dimension) {
			padded <<= 1;
		}
	}

	if (!option) {
		// Default behavior: multiply matrices in inputfile and print the diagonal.
		int* in[2];
		for (int i = 0; i < 2; i++) {
			in[i] = malloc(dimension*dimension*sizeof(int));
		}

		FILE* input = fopen(inputfile, "r");
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < dimension; j++) {
				for (int k = 0; k < dimension; k++) {
					int cell;
					if (fscanf(input, "%d", &cell) != 1) {
						fprintf(stderr, "failed to read inputfile %s\n", inputfile);
						return -1;
					}
					in[i][j*dimension + k] = cell;
				}
			}
		}
		fclose(input);

		struct matrix A = {
			in[0],
			dimension,
			dimension,
			padded,
			0,
			0
		}, B = {
			in[1],
			dimension,
			dimension,
			padded,
			0,
			0
		};
		struct matrix C = strassen(A, B);

		for (int i = 0; i < dimension; i++) {
			printf("%d\n", C.array[i*dimension + i]);
		}

		for (int i = 0; i < 2; i++) {
			free(in[i]);
		}
		free_matrix(C, NULL_MATRIX);
	} else if (option == -1) {
		// Optimize n0.
		srand(time(NULL));

		printf("  n0:  Time for naive:  Time for Strassen's:\n");
		for (int d = 2; d <= dimension; d++) {
			int* in[2];
			for (int i = 0; i < 2; i++) {
				in[i] = malloc(d*d*sizeof(int));				
			}

			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < d; j++) {
					for (int k = 0; k < d; k++) {
						in[i][j*d + k] = rand() - (RAND_MAX>>1);
					}
				}
			}

			size_t even = d + (d%2);
			struct matrix A = {
				in[0],
				d,
				d,
				even,
				0,
				0
			}, B = {
				in[1],
				d,
				d,
				even,
				0,
				0
			};

			int* trials[10];

			clock_t naive_ticks = clock();
			for (size_t i = 0; i < sizeof(trials)/sizeof(trials[0]); i++) {
				trials[i] = naive(A, B).array;
			}
			naive_ticks = clock() - naive_ticks;
			for (size_t i = 0; i < sizeof(trials)/sizeof(trials[0]); i++) {
				free(trials[i]);
			}

			// Only one level of Strassen's to compare if its addition is an improvement over the naive algorithm at the size d*d:
			n0 = d - 1;
			clock_t strassen_ticks = clock();
			for (size_t i = 0; i < sizeof(trials)/sizeof(trials[0]); i++) {
				trials[i] = strassen(A, B).array;
			}
			strassen_ticks = clock() - strassen_ticks;
			for (size_t i = 0; i < sizeof(trials)/sizeof(trials[0]); i++) {
				free(trials[i]);
			}

			printf("%4d:  %15f  %20f",
				   d,
				   ((double)    naive_ticks)/(sizeof(trials)/sizeof(trials[0]))/CLOCKS_PER_SEC,
				   ((double) strassen_ticks)/(sizeof(trials)/sizeof(trials[0]))/CLOCKS_PER_SEC);
			if (strassen_ticks < naive_ticks) {
				printf(" Won");
			} else if (naive_ticks < strassen_ticks) {
				printf(" Lost");
			}
			printf("\n");

			for (int i = 0; i < 2; i++) {
				free(in[i]);
			}
		}
	} else if (option > 0) {
		// Triangles.
		srand(time(NULL));

		float p = ((float) option)/100.0f;

		int* edges = malloc(1024*1024*sizeof(int));
		for (int i = 0; i < 1024; i++) {
			for (int j = 0; j < i; j++) {
				edges[i*1024 + j] = (rand() <= p*RAND_MAX) ? 1 : 0;
				// Undirected graph:
				edges[j*1024 + i] = edges[i*1024 + j];
			}
			edges[i*1024 + i] = 0;
		}

		struct matrix A = {
			edges,
			1024,
			1024,
			1024,
			0,
			0
		};
		struct matrix Asquared = strassen(A, A);
		struct matrix Acubed = strassen(Asquared, A);

		int triangles = 0;
		for (int i = 0; i < 1024; i++) {
			triangles += Acubed.array[i*1024 + i];
		}
		triangles /= 6;
		printf("%d triangles\n", triangles);

		free(edges);
		free_matrix(Asquared, NULL_MATRIX);
		free_matrix(Acubed, NULL_MATRIX);
	} else {
		fprintf(stderr, "unrecognized option %d\n", option);
		return -1;
	}

	return 0;
}
