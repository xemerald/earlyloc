#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
//#include <immintrin.h>
/* */
#include <matrix.h>

inline static int matrix_same_size( const MATRIX *, const MATRIX * );
inline static int matrix_square( const MATRIX * );
static int matrix_copy( MATRIX *, const MATRIX * );

/**/
MATRIX *matrix_new( const int row, const int column ) {
	MATRIX *res;

	res          = (MATRIX *)calloc(1, sizeof(MATRIX));
	res->i       = row;
	res->j       = column;
	res->total   = row * column;
	res->element = (double *)calloc(res->total, sizeof(double));

	return res;
}

/**/
MATRIX *matrix_identity( const int rank ) {
	MATRIX *res = matrix_new( rank, rank );
	int     i;
	int     step = res->j + 1;

	for ( i=0; i<res->total; i+=step ) res->element[i] = 1.0;

	return res;
}

/**/
MATRIX *matrix_add( const MATRIX *a, const MATRIX *b ) {
	int     i, j;
	MATRIX *res = NULL;

	if ( matrix_same_size(a, b) ) {
		res = matrix_new( a->i, a->j );

		for ( i=0; i<a->i; i++ ) {
			for ( j=0; j<a->j; j++ )
				res->element[i * a->j + j] = a->element[i * a->j + j] + b->element[i * a->j + j];
		}
	}

	return res;
}

/**/
MATRIX *matrix_sub( const MATRIX *a, const MATRIX *b ) {
	int     i, j;
	MATRIX *res = NULL;

	if ( matrix_same_size(a, b) ) {
		res = matrix_new( a->i, a->j );

		for ( i=0; i<a->i; i++ ) {
			for ( j=0; j<a->j; j++ )
				res->element[i * a->j + j] = a->element[i * a->j + j] - b->element[i * a->j + j];
		}
	}

	return res;
}

/**/
MATRIX *matrix_mul( const MATRIX *a, const MATRIX *b ) {
	int     i, j, k;
	double  tmp;
	MATRIX *res = NULL;

	if ( a->j == b->i ) {
		res = matrix_new( a->i, b->j );
		for ( i=0; i<a->i; i++ ) {
			for ( k=0; k<a->j; k++ ) {
				tmp = a->element[i * a->j + k];
				if ( fabs(tmp) > DBL_EPSILON ) {
					for ( j=0; j<b->j; j++ ) {
						if ( fabs(b->element[k * b->j + j]) > DBL_EPSILON )
							res->element[i * b->j + j] += tmp * b->element[k * b->j + j];
					}
				}
			}
		}
	}

	return res;
}

/*
MATRIX *matrix_mul_simd( const MATRIX *a, const MATRIX *b ) {
	int     i, j, k;
	double  tmp;
	MATRIX *res = NULL;
	__m256d I[8], R[8], S[8], Sum[8];


	if ( a->j == b->i ) {
		res = matrix_new( a->i, b->j );
		for ( i=0; i<8; i++ ) Sum[i] = _mm256_setzero_pd();
		for ( i=0; i<a->i; i++ ) {
			for ( k=0; k<8; k++ )
				R[k] = _mm256_loadu_pd()
			for ( k=0; k<a->j; k++ ) {
				tmp = a->element[i * a->j + k];
				if ( fabs(tmp) > DBL_EPSILON ) {
					for ( j=0; j<b->j; j++ ) {
						if ( fabs(b->element[k * b->j + j]) > DBL_EPSILON )
							res->element[i * b->j + j] += tmp * b->element[k * b->j + j];
					}
				}
			}
		}
	}

	return res;
}
*/
/**/
MATRIX *matrix_transpose( const MATRIX *a ) {
	int     i, j;
	MATRIX *res = NULL;

	res = matrix_new( a->j, a->i );

	for ( i=0; i<a->i; i++ ) {
		for ( j=0; j<a->j; j++ ) {
			res->element[j * a->i + i] = a->element[i * a->j + j];
		}
	}

	return res;
}

/*
MATRIX *matrix_adjugate( const MATRIX *a ) {
	int     i, j;
	MATRIX *res = NULL;

	if ( matrix_square( a ) ) {
		res = matrix_init( a->i, a->j );
		for ( i=0; i<a->i; i++ ) {
			for ( j=0; j<a->j; j++ ) {
				res->element[j*(a->i) + i] = a->element[i*(a->j) + j];
			}
		}
	}

	return res;
}
*/
/**/
MATRIX *matrix_inverse( const MATRIX *a ) {
	int     i, j, k;
	int     prow = 0;
	double  pivot;
	MATRIX *res = NULL;
	MATRIX *tmp = NULL;

	if ( matrix_square( a ) ) {
		tmp = matrix_new( a->i, a->j );
		res = matrix_identity( a->i );
		matrix_copy( tmp, a );

		for ( i=0; i<a->i; i++ ) {
		/* */
			pivot = tmp->element[i * a->j + i];
			for ( j=i+1, prow=i; j<a->i; j++ ) {
				if ( tmp->element[j * a->j + i] > pivot ) {
					pivot = tmp->element[j * a->j + i];
					prow  = j;
				}
			}
		/* */
			if ( prow != i ) {
				for ( j=0; j<a->j; j++ ) {
				/* */
					pivot                         = tmp->element[i * a->j + j];
					tmp->element[i * a->j + j]    = tmp->element[prow * a->j + j];
					tmp->element[prow * a->j + j] = pivot;
				/* */
					pivot                         = res->element[i * a->j + j];
					res->element[i * a->j + j]    = res->element[prow * a->j + j];
					res->element[prow * a->j + j] = pivot;
				}
			}
		/* */
			pivot = tmp->element[i * a->j + i];
			for ( j=0; j<a->j; j++ ) {
			/* */
				tmp->element[i * a->j + j] /= pivot;
			/* */
				res->element[i * a->j + j] /= pivot;
			}
		/* */
			for ( j=i+1; j<a->i; j++ ) {
				pivot = tmp->element[j * a->j + i];
				if ( fabs(pivot) > DBL_EPSILON ) {
					for ( k=0; k<a->j; k++ ) {
					/* */
						if ( fabs(tmp->element[i * a->j + k]) > DBL_EPSILON )
							tmp->element[j * a->j + k] -= pivot * tmp->element[i * a->j + k];
					/* */
						if ( fabs(res->element[i * a->j + k]) > DBL_EPSILON )
							res->element[j * a->j + k] -= pivot * res->element[i * a->j + k];
					}
				}
			}
		}

		prow = a->i - 1;
		for ( i=0; i<prow; i++ ) {
			for ( j=i+1; j<a->i; j++ ) {
				pivot = tmp->element[i * a->j + j];
				if ( fabs(pivot) > DBL_EPSILON ) {
					for ( k=0; k<a->j; k++ ) {
					/* */
						if ( fabs(tmp->element[j * a->j + k]) > DBL_EPSILON )
							tmp->element[i * a->j + k] -= pivot * tmp->element[j * a->j + k];
					/* */
						if ( fabs(res->element[j * a->j + k]) > DBL_EPSILON )
							res->element[i * a->j + k] -= pivot * res->element[j* a->j + k];
					}
				}
			}
		}

		matrix_free( tmp );
	}

	return res;
}

/* Least square */
MATRIX *matrix_div( const MATRIX *a, const MATRIX *b ) {
	MATRIX *gt   = matrix_transpose( b );
	MATRIX *gtg  = matrix_mul( gt, b );
	MATRIX *gtd  = matrix_mul( gt, a );
	MATRIX *igtg = matrix_inverse( gtg );
	MATRIX *res  = matrix_mul( igtg, gtd );

	matrix_free( gt );
	matrix_free( gtg );
	matrix_free( gtd );
	matrix_free( igtg );

	return res;
}

/*
 *
 */
MATRIX *matrix_div_weighted( const MATRIX *a, const MATRIX *b, const MATRIX *w ) {
	MATRIX *gt    = matrix_transpose( b );
	MATRIX *wg    = matrix_mul( w, b );
	MATRIX *gtwg  = matrix_mul( gt, wg );
	MATRIX *wd    = matrix_mul( w, a );
	MATRIX *gtwd  = matrix_mul( gt, wd );
	MATRIX *igtwg = matrix_inverse( gtwg );
	MATRIX *res   = matrix_mul( gtwg, gtwd );

	matrix_free( gt );
	matrix_free( wg );
	matrix_free( gtwg );
	matrix_free( wd );
	matrix_free( gtwd );
	matrix_free( igtwg );

	return res;
}

/* Assignment functions */
/*
 *
 */
MATRIX *matrix_assign( MATRIX *dest, const double src, int row, int col ) {
	row--;
	col--;

	if ( row < dest->i && row >= 0 && col < dest->j && col >= 0 ) {
		dest->element[row * dest->j + col] = src;
		return dest;
	}

	return NULL;
}

MATRIX *matrix_assign_seq( MATRIX *dest, const double *src, const int data_size ) {
	if ( dest->total >= data_size ) {
		memcpy(dest->element, src, data_size * sizeof(double));
		return dest;
	}

	return NULL;
}

MATRIX *matrix_assign_row( MATRIX *dest, const double *src, int row_index, const int data_size ) {
	int j;

	row_index--;

	if ( dest->j == data_size ) {
		memcpy(dest->element + row_index * dest->j, src, dest->j * sizeof(double));
		return dest;
	}
	else if ( dest->j > data_size ) {
		memcpy(dest->element + row_index * dest->j, src, data_size * sizeof(double));
		for ( j=data_size; j<dest->j; j++ )
			dest->element[row_index * dest->j + j] = 0.0;
		return dest;
	}

	return NULL;
}

MATRIX *matrix_assign_col( MATRIX *dest, const double *src, int col_index, const int data_size ) {
	int i;

	col_index--;

	if ( dest->i == data_size ) {
		for ( i=0; i<dest->i; i++ )
			dest->element[i * dest->j + col_index] = src[i];
		return dest;
	}
	else if ( dest->i > data_size ) {
		for ( i=0; i<dest->i; i++ ) {
			if ( i < data_size )
				dest->element[i * dest->j + col_index] = src[i];
			else
				dest->element[i * dest->j + col_index] = 0.0;
		}
		return dest;
	}

	return NULL;
}

MATRIX *matrix_assign_diag( MATRIX *dest, const double *src, const int data_size ) {
	int i;

	if ( matrix_square( dest ) ) {
		if ( dest->i == data_size ) {
			for ( i=0; i<dest->i; i++ )
				dest->element[i * dest->j + i] = src[i];
			return dest;
		}
		else if ( dest->i > data_size ) {
			for ( i=0; i<dest->i; i++ ) {
				if ( i < data_size )
					dest->element[i * dest->j + i] = src[i];
				else
					dest->element[i * dest->j + i] = 0.0;
			}
			return dest;
		}
	}

	return NULL;
}

MATRIX *matrix_apply_all( MATRIX *dest, double (*func)( const double ) ) {
	int i;

	for ( i=0; i<dest->total; i++ )
		dest->element[i] = func(dest->element[i]);

	return dest;
}

MATRIX *matrix_apply_row( MATRIX *dest, double (*func)( const double ), int row_index ) {
	int j;

	row_index--;
	for ( j=0; j<dest->j; j++ ) dest->element[row_index * dest->j + j] = func(dest->element[row_index * dest->j + j]);

	return dest;
}

MATRIX *matrix_apply_col( MATRIX *dest, double (*func)( const double ), int col_index ) {
	int i;

	col_index--;
	for ( i=0; i<dest->i; i++ ) dest->element[i * dest->j + col_index] = func(dest->element[i * dest->j + col_index]);

	return dest;
}

MATRIX *matrix_apply_diag( MATRIX *dest, double (*func)( const double ) ) {
	int i;
	int step = dest->j + 1;

	if ( matrix_square( dest ) )
		for ( i=0; i<dest->total; i+=step )
			dest->element[i] = func(dest->element[i]);
	else
		return NULL;

	return dest;
}

double *matrix_prefill_array( double *dest, int data_size, ... ) {
	va_list ap;
	double  dval;
	double *_dest = dest;

	va_start(ap, data_size);
	for ( ; data_size>0; data_size-- ) {
		dval     = va_arg(ap, double);
		*_dest++ = dval;
	}
	va_end(ap);

	return dest;
}


double *matrix_extract_seq( const MATRIX *src, double *dest, const int dest_size ) {
	if ( src->total <= dest_size ) {
		memcpy(dest, src->element, dest_size * sizeof(double));
		return dest;
	}

	return NULL;
}

double matrix_determinant( const MATRIX *a ) {
	if ( matrix_square( a ) ) {
		return 0;
	}
	return -1;
}

void matrix_free( MATRIX *a ) {
	free(a->element);
	free(a);
	return;
}

inline static int matrix_same_size( const MATRIX *a, const MATRIX *b ) {
	return ( a->i == b->i && a->j == b->j );
}

inline static int matrix_square( const MATRIX *a ) {
	return ( a->i == a->j );
}

static int matrix_copy( MATRIX *dest, const MATRIX *src ) {
	if ( matrix_same_size( dest, src ) ) {
		memcpy(dest->element, src->element, src->total * sizeof(double));
		return 0;
	}
	return -1;
}
