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

/* */
#define MATRIX_DIAG_EPS  1.0e-5

#define ARE_SAME_SIZE( _MATRIX_A, _MATRIX_B ) \
		((_MATRIX_A)->i == (_MATRIX_B)->i && (_MATRIX_A)->j == (_MATRIX_B)->j)

#define IS_SQUARE( _MATRIX ) \
		((_MATRIX)->i == (_MATRIX)->j)

#define ADD_EPS_TO_DIAG( _MATRIX ) \
		__extension__({ \
			if ( IS_SQUARE( (_MATRIX) ) ) \
				for ( int i = 0; i < (_MATRIX)->i; i++ ) \
					if ( fabs((_MATRIX)->element[i * (_MATRIX)->j + i]) < MATRIX_DIAG_EPS ) \
						(_MATRIX)->element[i * (_MATRIX)->j + i] += (_MATRIX)->element[i * (_MATRIX)->j + i] > 0.0 ? MATRIX_DIAG_EPS : -MATRIX_DIAG_EPS; \
		})

/***/
static MATRIX *duplicate_matrix( const MATRIX * );

/***/
MATRIX *matrix_new( const int row, const int column ) {
	MATRIX *res = (MATRIX *)calloc(1, sizeof(MATRIX));

	if ( res ) {
		res->i       = row;
		res->j       = column;
		res->total   = row * column;
		if ( !(res->element = (double *)calloc(res->total, sizeof(double))) ) {
			free(res);
			res = NULL;
		}
	}

	return res;
}

/***/
MATRIX *matrix_identity( const int rank ) {
	MATRIX   *res = matrix_new( rank, rank );
	const int step = res->j + 1;

	if ( res )
		for ( long i = 0; i < res->total; i += step ) 
			res->element[i] = 1.0;

	return res;
}

/***/
MATRIX *matrix_add( const MATRIX *a, const MATRIX *b ) {
	MATRIX *res = NULL;

	if ( ARE_SAME_SIZE( a, b ) && (res = matrix_new( a->i, a->j )) ) {
		for ( int i = 0; i < a->i; i++ )
			for ( int j = 0; j < a->j; j++ )
				res->element[i * a->j + j] = a->element[i * a->j + j] + b->element[i * a->j + j];
	}

	return res;
}

/***/
MATRIX *matrix_sub( const MATRIX *a, const MATRIX *b ) {
	MATRIX *res = NULL;

	if ( ARE_SAME_SIZE( a, b ) && (res = matrix_new( a->i, a->j )) ) {
		for ( int i = 0; i < a->i; i++ )
			for ( int j = 0; j < a->j; j++ )
				res->element[i * a->j + j] = a->element[i * a->j + j] - b->element[i * a->j + j];
	}

	return res;
}

/***/
MATRIX *matrix_mul( const MATRIX *a, const MATRIX *b ) {
	double  tmp;
	MATRIX *res = NULL;

	if ( a->j == b->i && (res = matrix_new( a->i, b->j )) ) {
		for ( int i = 0; i < a->i; i++ ) {
			for ( int k = 0; k < a->j; k++ ) {
				if ( fabs(tmp = a->element[i * a->j + k]) > DBL_EPSILON ) {
					for ( int j = 0; j < b->j; j++ ) {
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

/***/
MATRIX *matrix_transpose( const MATRIX *a ) {
	MATRIX *res = matrix_new( a->j, a->i );

	if ( res ) {
		for ( int i = 0; i < a->i; i++ )
			for ( int j = 0; j < a->j; j++ )
				res->element[j * a->i + i] = a->element[i * a->j + j];
	}

	return res;
}

/*
MATRIX *matrix_adjugate( const MATRIX *a ) {
	MATRIX *res = NULL;
	if ( IS_SQUARE( a ) ) {
		res = matrix_init( a->i, a->j );
		for ( int i = 0; i < a->i; i++ ) {
			for ( int j = 0; j < a->j; j++ ) {
				res->element[j*(a->i) + i] = a->element[i*(a->j) + j];
			}
		}
	}

	return res;
}
*/

/***/
MATRIX *matrix_inverse( const MATRIX *a ) 
{
	int     prow = 0;
	double  pivot;
	MATRIX *res = NULL;
	MATRIX *tmp = NULL;

	if ( IS_SQUARE( a ) ) {
		tmp = duplicate_matrix( a );
		res = matrix_identity( a->i );

		for ( int i = 0; i < a->i; i++ ) {
		/* */
			pivot = tmp->element[i * a->j + i];
			prow = i;
			for ( int j = i + 1; j < a->i; j++ ) {
				if ( tmp->element[j * a->j + i] > pivot ) {
					pivot = tmp->element[j * a->j + i];
					prow  = j;
				}
			}
		/* */
			if ( prow != i ) {
				for ( int j = 0; j < a->j; j++ ) {
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
			for ( int j = 0; j < a->j; j++ ) {
			/* */
				tmp->element[i * a->j + j] /= pivot;
			/* */
				res->element[i * a->j + j] /= pivot;
			}
		/* */
			for ( int j = i + 1; j < a->i; j++ ) {
				pivot = tmp->element[j * a->j + i];
				if ( fabs(pivot) > DBL_EPSILON ) {
					for ( int k = 0; k < a->j; k++ ) {
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
	/* */
		prow = a->i - 1;
		for ( int i = 0; i < prow; i++ ) {
			for ( int j = i + 1; j < a->i; j++ ) {
				pivot = tmp->element[i * a->j + j];
				if ( fabs(pivot) > DBL_EPSILON ) {
					for ( int k = 0; k < a->j; k++ ) {
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
	/* */
		matrix_free( tmp );
	}

	return res;
}

/* Least square */
MATRIX *matrix_div( const MATRIX *a, const MATRIX *b ) {
	MATRIX *gt;
	MATRIX *gtg;
	MATRIX *gtd;
	MATRIX *igtg;
	MATRIX *res;

/* */
	gt  = matrix_transpose( b );
	gtg = matrix_mul( gt, b );
	gtd = matrix_mul( gt, a );
/* */
	ADD_EPS_TO_DIAG( gtg );
	igtg = matrix_inverse( gtg );
	res  = matrix_mul( igtg, gtd );

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
	MATRIX *wg;
	MATRIX *gtw;
	MATRIX *gtwg;
	MATRIX *gtwd;
	MATRIX *igtwg;
	MATRIX *res;

/* */
	wg = duplicate_matrix( b );
	for ( int i = 0; i < wg->i; i++ ) {
		double tmp = w->element[i * w->j + i];
		for ( int j = 0; j < wg->j; j++ ) 
			wg->element[i * wg->j + j] *= tmp;
	}
/* */
	gtw  = matrix_transpose( wg );
	gtwg = matrix_mul( gtw, b );
	gtwd = matrix_mul( gtw, a );
/* */
	ADD_EPS_TO_DIAG( gtwg );
	igtwg = matrix_inverse( gtwg );
	res   = matrix_mul( igtwg, gtwd );
/* */
	matrix_free( wg );
	matrix_free( gtw );
	matrix_free( gtwg );
	matrix_free( gtwd );
	matrix_free( igtwg );

	return res;
}

/* Assignment functions */

/***/
MATRIX *matrix_assign( MATRIX *dest, const double src, int row, int col ) 
{
	row--;
	col--;

	if ( row < dest->i && row >= 0 && col < dest->j && col >= 0 ) {
		dest->element[row * dest->j + col] = src;
		return dest;
	}

	return NULL;
}

/***/
MATRIX *matrix_assign_seq( MATRIX *dest, const double *src, const long data_size ) 
{
	if ( dest->total >= data_size ) {
		memcpy(dest->element, src, data_size * sizeof(double));
		return dest;
	}

	return NULL;
}

/***/
MATRIX *matrix_assign_row( MATRIX *dest, const double *src, int row_index, const int data_size ) 
{
	row_index--;

	if ( dest->j == data_size ) {
		memcpy(dest->element + row_index * dest->j, src, dest->j * sizeof(double));
		return dest;
	}
	else if ( dest->j > data_size ) {
		memcpy(dest->element + row_index * dest->j, src, data_size * sizeof(double));
		for ( int j = data_size; j < dest->j; j++ )
			dest->element[row_index * dest->j + j] = 0.0;
		return dest;
	}

	return NULL;
}

/***/
MATRIX *matrix_assign_col( MATRIX *dest, const double *src, int col_index, const int data_size ) 
{
	col_index--;

	if ( dest->i == data_size ) {
		for ( int i = 0; i < dest->i; i++ )
			dest->element[i * dest->j + col_index] = src[i];
		return dest;
	}
	else if ( dest->i > data_size ) {
		for ( int i = 0; i < dest->i; i++ ) {
			if ( i < data_size )
				dest->element[i * dest->j + col_index] = src[i];
			else
				dest->element[i * dest->j + col_index] = 0.0;
		}
		return dest;
	}

	return NULL;
}

/***/
MATRIX *matrix_assign_diag( MATRIX *dest, const double *src, const int data_size ) 
{
	if ( IS_SQUARE( dest ) ) {
		if ( dest->i == data_size ) {
			for ( int i = 0; i < dest->i; i++ )
				dest->element[i * dest->j + i] = src[i];
			return dest;
		}
		else if ( dest->i > data_size ) {
			for ( int i = 0; i < dest->i; i++ ) {
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

/***/
MATRIX *matrix_apply_all( MATRIX *dest, double (*func)( const double ) ) 
{
	for ( long i = 0; i < dest->total; i++ )
		dest->element[i] = func(dest->element[i]);

	return dest;
}

/***/
MATRIX *matrix_apply_row( MATRIX *dest, double (*func)( const double ), int row_index ) 
{
	row_index--;
	for ( int j = 0; j < dest->j; j++ ) 
		dest->element[row_index * dest->j + j] = func(dest->element[row_index * dest->j + j]);

	return dest;
}

/***/
MATRIX *matrix_apply_col( MATRIX *dest, double (*func)( const double ), int col_index ) 
{
	col_index--;
	for ( int i = 0; i < dest->i; i++ ) 
		dest->element[i * dest->j + col_index] = func(dest->element[i * dest->j + col_index]);

	return dest;
}

/***/
MATRIX *matrix_apply_diag( MATRIX *dest, double (*func)( const double ) ) 
{
	const int step = dest->j + 1;

	if ( IS_SQUARE( dest ) )
		for ( long i = 0; i < dest->total; i += step )
			dest->element[i] = func(dest->element[i]);
	else
		return NULL;

	return dest;
}

/***/
double *matrix_prefill_array( double *dest, int data_size, ... ) 
{
	va_list ap;
	double  dval;
	double *_dest = dest;

	va_start(ap, data_size);
	for ( ; data_size > 0; data_size-- ) {
		dval     = va_arg(ap, double);
		*_dest++ = dval;
	}
	va_end(ap);

	return dest;
}

/***/
double *matrix_extract_seq( const MATRIX *src, double *dest, const long dest_size ) 
{
	if ( src->total <= dest_size ) {
		memcpy(dest, src->element, dest_size * sizeof(double));
		return dest;
	}

	return NULL;
}

/***/
double matrix_determinant( const MATRIX *a ) 
{
	if ( IS_SQUARE( a ) ) {
		return 0;
	}
	return -1;
}

/***/
void matrix_free( MATRIX *a ) 
{
	free(a->element);
	free(a);
	return;
}

/***/
static MATRIX *duplicate_matrix( const MATRIX *src ) 
{
	MATRIX *result = matrix_new( src->i, src->j );
	
	if ( result )
		memcpy(result->element, src->element, src->total * sizeof(double));
		
	return result;
}
