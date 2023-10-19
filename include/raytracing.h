/**
 * @file raytracing.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-10-13
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
/* */
#define RT_MAX_NODE  16384
/* */
#define RT_P_WAVE_VELOCITY  0
#define RT_S_WAVE_VELOCITY  1

/**
 * @brief
 *
 */
typedef struct {
	double r;
	double a;
	double b;
	double v;
} RAY_INFO;

/* */
int rt_main( RAY_INFO *, int *, double *, double, double, double, double, double, double, const int );
int rt_velmod_load( const char * );
void rt_velmod_free( void );
void rt_tko_azi_cal( const RAY_INFO *, const int, double *, double * );
void rt_drvt_cal( const RAY_INFO *, const int, double *, double *, double * );
/* */
RAY_INFO *rt_ray_rad2deg( RAY_INFO *, const int, double, double );
