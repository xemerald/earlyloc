/**
 * @file earlyloc_misc.h
 * @author Benjamin Yang in Department of Geology, National Taiwan University
 * @brief
 * @version 0.1
 * @date 2023-10-05
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#define EL_SIMPLE_TIMESTAMP_BUF_SIZE    16

/* */
double el_misc_timenow( void );
char  *el_misc_simple_timestamp_gen( char *, const int, const double );
double el_misc_geog2distf( const double, const double, const double, const double );
