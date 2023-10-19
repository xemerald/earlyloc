/**
 * @file matrix.h
 * @author Benjamin Yang in Department of Geology, National Taiwan University
 * @brief
 * @version 0.1
 * @date 2018-03-05
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

/**
 * @brief
 *
 */
typedef struct {
	int     i;
	int     j;
	int     total;
	double *element;
} MATRIX;

/* */
MATRIX *matrix_new( const int, const int );
MATRIX *matrix_identity( const int );

MATRIX *matrix_add( const MATRIX *, const MATRIX * );
MATRIX *matrix_sub( const MATRIX *, const MATRIX * );
MATRIX *matrix_mul( const MATRIX *, const MATRIX * );
MATRIX *matrix_div( const MATRIX *, const MATRIX * );
MATRIX *matrix_div_weighted( const MATRIX *, const MATRIX *, const MATRIX * );
MATRIX *matrix_transpose( const MATRIX * );
MATRIX *matrix_inverse( const MATRIX * );

MATRIX *matrix_assign( MATRIX *, const double, int, int );
MATRIX *matrix_assign_seq( MATRIX *, const double *, const int );
MATRIX *matrix_assign_row( MATRIX *, const double *, int, const int );
MATRIX *matrix_assign_col( MATRIX *, const double *, int, const int );
MATRIX *matrix_assign_diag( MATRIX *, const double *, const int );

MATRIX *matrix_apply_all( MATRIX *, double (*)( const double ) );
MATRIX *matrix_apply_row( MATRIX *, double (*)( const double ), int );
MATRIX *matrix_apply_col( MATRIX *, double (*)( const double ), int );
MATRIX *matrix_apply_diag( MATRIX *, double (*)( const double ) );

double *matrix_prefill_array( double *, int, ... );
double *matrix_extract_seq( const MATRIX *, double *, const int );

double  matrix_determinant( const MATRIX * );

void matrix_free( MATRIX * );
