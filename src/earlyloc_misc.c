/**
 * @file earlyloc_misc.c
 * @author Benjamin Yang in Department of Geology, National Taiwan University
 * @brief
 * @version 0.1
 * @date 2023-10-05
 *
 * @copyright Copyright (c) 2023
 *
 */
/* Standard C header include */
#include <stdio.h>
#include <time.h>
#include <math.h>
/* */
#include <earlyloc_misc.h>

#define TIMESTAMP_SIMPLE_FORMAT  "%04d%02d%02d%02d%02d%02d"

/**
 * @brief
 *
 * @return double
 */
double el_misc_timenow( void )
{
	struct timespec time_sp;

/* */
	clock_gettime(CLOCK_REALTIME_COARSE, &time_sp);

	return time_sp.tv_sec + time_sp.tv_nsec * 1.0e-9;
}

/*
 *
 */
char *el_misc_simple_timestamp_gen( char *buffer, const int buffer_size, const double timestamp )
{
	struct tm    sptime;
	const time_t _timestamp = (time_t)timestamp;

/* Checking for buffer size */
	if ( buffer == NULL || buffer_size < EL_SIMPLE_TIMESTAMP_BUF_SIZE )
		return buffer;
/* */
	gmtime_r(&_timestamp, &sptime);
	sprintf(
		buffer, TIMESTAMP_SIMPLE_FORMAT,
		sptime.tm_year + 1900, sptime.tm_mon + 1, sptime.tm_mday,
		sptime.tm_hour, sptime.tm_min, sptime.tm_sec
	);

	return buffer;
}

/**
 * @brief Transforms the geographic coordinate(latitude & longitude) into distance(unit: km)
 *
 * @param elon
 * @param elat
 * @param slon
 * @param slat
 * @return double
 */
double el_misc_geog2distf( const double elon, const double elat, const double slon, const double slat )
{
	const double avlat = (elat + slat)*0.5;

	double a = 1.840708 + avlat*(.0015269 + avlat*(-.00034 + avlat*(1.02337e-6)));
	double b = 1.843404 + avlat*(-6.93799e-5 + avlat*(8.79993e-6 + avlat*(-6.47527e-8)));

	a *= (slon - elon) * 60.0;
	b *= (slat - elat) * 60.0;

	return sqrt(a * a + b * b);
}
