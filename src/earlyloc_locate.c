/**
 * @file earlyloc_locate.c
 * @author Benjamin Yang in Department of Geology, National Taiwan University
 * @brief
 * @version 0.1
 * @date 2023-10-05
 *
 * @copyright Copyright (c) 2023
 *
 */
/* */
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* */
#include <constants.h>
#include <raytracing.h>
#include <dl_chain_list.h>
#include <matrix.h>
#include <earlyloc.h>
#include <earlyloc_misc.h>

/* */
#define GEIGER_ERROR_RETURN  -1.0f
/* */
#define REFINE_DAMPING_FACTOR  4.0f
#define REFINE_MIN_DAMPING     0.5f
#define NEAR_HYPO_DISTANCE     80.0f
#define FAR_HYPO_DISTANCE      600.0f
/*
 *
 */
typedef struct {
/* Hypo & station point */
	double hyp_depth;
	double epc_dist;
	double epc_dist_x;
	double epc_dist_y;
/* Velocity parameters */
	double vel_init;
	double vel_grad;
/* Center of the circle */
	double center_x;
	double center_z;
/* Angles */
	double angle_a;
	double angle_b;
} LINEAR_RAY_PATH_PARAMS;

/* */
static double step_geiger_method(
	double *, double *, double *, double *, PICKS_POOL *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL *, const double [HYPO_PARAMS_NUMBER]
);
static double step_geiger_method_3D( double *, double *, double *, double *, PICKS_POOL *, const double [HYPO_PARAMS_NUMBER] );
static void update_picks_state(
	const double, const double, const double, const double, PICKS_POOL *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL *
);
static void update_picks_state_3D( const double, const double, const double, const double, PICKS_POOL * );
static void update_hypo_state( const double, const double, const double, const double, HYPO_STATE * );
static LINEAR_RAY_PATH_PARAMS get_linear_ray_path(
	const double, const double, const double, const double, const double, const double, const double, const double, const double
);
static double  get_linear_travel_time( const LINEAR_RAY_PATH_PARAMS * );
static double *get_travel_time_derivatives( const LINEAR_RAY_PATH_PARAMS *, double [HYPO_PARAMS_NUMBER] );
static double *get_travel_time_derivatives_3D( const RAY_INFO *, const int, double [HYPO_PARAMS_NUMBER] );
static double  get_r_weight( const double, const double, const double, const int, const int );
static double *get_damping_matrix( double [HYPO_PARAMS_NUMBER], const HYPO_STATE * );
static double  get_gap_degree( const double, const double, const PICKS_POOL * );
static int get_hypo_quality( const int, const double, const double, const double );
static int compare_gap( const void *, const void * );

/* */
#define SELECT_VEL_DEPTH(VEL_INIT, VEL_GRAD, DEPTH, MODEL) \
		__extension__({ \
			(VEL_INIT) = (DEPTH) < (MODEL)->boundary ? (MODEL)->shallow_init : (MODEL)->deep_init; \
			(VEL_GRAD) = (DEPTH) < (MODEL)->boundary ? (MODEL)->shallow_grad : (MODEL)->deep_grad; \
		})

#define SELECT_VEL_PHASE_DEPTH(VEL_INIT, VEL_GRAD, DEPTH, PHASE, PMODEL, SMODEL) \
		__extension__({ \
			const LAYER_VEL_MODEL *model_ptr_in_macro = NULL; \
			if ( !strcmp((PHASE), "P") ) \
				model_ptr_in_macro = (PMODEL); \
			else if ( !strcmp((PHASE), "S") ) \
				model_ptr_in_macro = (SMODEL); \
			if ( model_ptr_in_macro ) { \
				SELECT_VEL_DEPTH((VEL_INIT), (VEL_GRAD), (DEPTH), model_ptr_in_macro); \
			} \
		})

#define INIT_RAY_PATH_PARAMS(__RAY_PATH) \
		((__RAY_PATH) = (LINEAR_RAY_PATH_PARAMS){ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 })

/* */
static uint8_t VelocityModel3DReady = 0;

/*
 *
 */
int el_loc_primary_locate( HYPO_STATE *hyp, const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model )
{
	double lon0, lat0, depth0, time0;
	double tmp;

/* */
	lon0   = hyp->longitude;
	lat0   = hyp->latitude;
	depth0 = hyp->depth;
	time0  = hyp->origin_time;

/* First, Geiger's method step for 10 times */
	for ( int i = 0; i < 20; i++ ) {
		if ( VelocityModel3DReady )
			tmp = step_geiger_method_3D( &lon0, &lat0, &depth0, &time0, &hyp->pool, NULL );
		else
			tmp = step_geiger_method( &lon0, &lat0, &depth0, &time0, &hyp->pool, p_model, s_model, NULL );
	/* */
		if ( tmp <= 1.0 )
			break;
	/* Go back to the initial if the result is too far */
		if ( fabs(time0 - hyp->origin_time) > 60.0 ) {
			lon0   = hyp->longitude;
			lat0   = hyp->latitude;
			depth0 = hyp->depth;
			time0  = hyp->origin_time;
		}
	}
/* Something error when doing geiger method */
	if ( tmp < 0.0 )
		return -1;
/* Calculate relative parameters one more by new hyp depth */
	if ( VelocityModel3DReady )
		update_picks_state_3D( lon0, lat0, depth0, time0, &hyp->pool );
	else
		update_picks_state( lon0, lat0, depth0, time0, &hyp->pool, p_model, s_model );
/* Update all the parameters to the hypo state */
	update_hypo_state( lon0, lat0, depth0, time0, hyp );

	return 0;
}

/*
 *
 */
void el_loc_all_states_update( HYPO_STATE *hyp, const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model )
{
/* Calculate relative parameters one more by new hyp depth */
	if ( VelocityModel3DReady )
		update_picks_state_3D( hyp->longitude, hyp->latitude, hyp->depth, hyp->origin_time, &hyp->pool );
	else
		update_picks_state( hyp->longitude, hyp->latitude, hyp->depth, hyp->origin_time, &hyp->pool, p_model, s_model );
/* Update all the parameters to the hypo state */
	update_hypo_state( hyp->longitude, hyp->latitude, hyp->depth, hyp->origin_time, hyp );

	return;
}

/*
 *
 */
void el_loc_location_guess( HYPO_STATE *hyp, const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model )
{
	DL_NODE    *node;
	PICK_STATE *pick, *first = NULL;
	int         wfactor   = 0;
	int         num_pick  = 0;
	double      avg_lon   = 0.0;
	double      avg_lat   = 0.0;
	double      avg_ptime = 0.0;

/* Find the first pick & the average arrival time */
	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) && !strcmp(pick->observe.phase_name, "P") ) {
			avg_ptime += pick->observe.picktime;
			num_pick++;
			if ( !first || pick->observe.picktime < first->observe.picktime )
				first = pick;
		}
	}
	avg_ptime /= num_pick;
/* */
	num_pick = 0;
	hyp->origin_time = -1.0;
	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) && !strcmp(pick->observe.phase_name, "P") && pick->observe.picktime < avg_ptime ) {
		/* Those primary picks whould have more weighting */
			wfactor   = pick->flag & PICK_FLAG_PRIMARY ? 3 : 1;
			wfactor  += pick == first ? 1 : 0;
			avg_lon  += pick->observe.longitude * wfactor;
			avg_lat  += pick->observe.latitude * wfactor;
			num_pick += wfactor;
			if ( hyp->origin_time < 0.0 || pick->observe.picktime < hyp->origin_time )
				hyp->origin_time = pick->observe.picktime;
		}
	}

/* */
	if ( num_pick ) {
		hyp->origin_time -= 2.0;
		hyp->longitude    = avg_lon / num_pick;
		hyp->latitude     = avg_lat / num_pick;
	}
	else {
		hyp->origin_time = first->observe.picktime - 2.0;
		hyp->longitude   = first->observe.longitude + 0.01;
		hyp->latitude    = first->observe.latitude + 0.01;
	}
	hyp->depth = 10.0;

	return;
}

/*
 *
 */
void el_loc_location_refine( HYPO_STATE *hyp, const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model )
{
	double lon0, lat0, depth0, time0;
	double dmatrix[HYPO_PARAMS_NUMBER];

/* Find damping factors */
	get_damping_matrix( dmatrix, hyp );
/* */
	lon0   = hyp->longitude;
	lat0   = hyp->latitude;
	depth0 = 10.0;
	time0  = hyp->origin_time;
/* */
	update_picks_state( lon0, lat0, depth0, time0, &hyp->pool, p_model, s_model );
	update_hypo_state( lon0, lat0, depth0, time0, hyp );
/* */
	for ( int i = 0; i < 30; i++ ) {
	/* */
		step_geiger_method( &lon0, &lat0, &depth0, &time0, &hyp->pool, p_model, s_model, dmatrix );
	/* Update all the parameters to the hypo state */
		update_picks_state( lon0, lat0, depth0, time0, &hyp->pool, p_model, s_model );
		update_hypo_state( lon0, lat0, depth0, time0, hyp );
	#ifdef _DEBUG
		printf("earlyloc: Refine: %lf %lf %lf %lf, average error: %lf\n", lon0, lat0, depth0, time0, hyp->avg_error );
	#endif
	/* */
		if ( hyp->avg_error < 0.9 )
			break;
	}

	return;
}

/**
 * @brief
 *
 * @param hyp
 * @param p_model
 * @param s_model
 * @return double
 */
double el_loc_origintime_adjust( HYPO_STATE *hyp, const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model )
{
	DL_NODE    *node;
	PICK_STATE *pick;

	double result   = 0.0;
	double residual = 0.0;
	double r_weight = 0.0;
	double sum_wei  = 0.0;
	double veli     = 0.0;
	double velg     = 0.0;
	LINEAR_RAY_PATH_PARAMS ray_path;
/* */
	const double delta_x = el_misc_geog2distf( hyp->longitude - 0.5, hyp->latitude, hyp->longitude + 0.5, hyp->latitude );
	const double delta_y = el_misc_geog2distf( hyp->longitude, hyp->latitude - 0.5, hyp->longitude, hyp->latitude + 0.5 );

	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
		/* */
			SELECT_VEL_PHASE_DEPTH( veli, velg, hyp->depth, pick->observe.phase_name, p_model, s_model );
		/* */
			ray_path = get_linear_ray_path(
				hyp->longitude, hyp->latitude, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, hyp->depth, veli, velg
			);
		/* Assign derived values to picking */
			residual = pick->observe.picktime - (hyp->origin_time + get_linear_travel_time( &ray_path ));
			r_weight = get_r_weight( ray_path.epc_dist, hyp->depth, residual, pick->observe.weight, pick->flag );
		/* */
			sum_wei += r_weight;
			result  += residual * r_weight;
		}
	}
/* Recalculate the travel time residual & weight */
	result /= sum_wei;
	hyp->origin_time += result;

	return result;
}

/**
 * @brief
 *
 * @param hyp
 * @param pick
 * @param p_model
 * @param s_model
 * @return double
 */
double el_loc_residual_estimate( const HYPO_STATE *hyp, const PICK_STATE *pick, const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model )
{
	LINEAR_RAY_PATH_PARAMS ray_path;
	double veli = 0.0;
	double velg = 0.0;
/* */
	const double delta_x = el_misc_geog2distf( hyp->longitude - 0.5, hyp->latitude, hyp->longitude + 0.5, hyp->latitude );
	const double delta_y = el_misc_geog2distf( hyp->longitude, hyp->latitude - 0.5, hyp->longitude, hyp->latitude + 0.5 );

/* */
	INIT_RAY_PATH_PARAMS( ray_path );
	SELECT_VEL_PHASE_DEPTH( veli, velg, hyp->depth, pick->observe.phase_name, p_model, s_model );
/* */
	ray_path = get_linear_ray_path(
		hyp->longitude, hyp->latitude, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, hyp->depth, veli, velg
	);
/* */

	return pick->observe.picktime - (hyp->origin_time + get_linear_travel_time( &ray_path ));
}

/**
 * @brief
 *
 * @param model_path
 * @return int
 */
int el_loc_3dvelmod_load( const char *model_path )
{
	if ( rt_velmod_load( model_path ) )
		return -1;

	VelocityModel3DReady = 1;

	return 0;
}

/**
 * @brief
 *
 */
void el_loc_3dvelmod_free( void )
{
	if ( VelocityModel3DReady )
		rt_velmod_free();

	return;
}

/*
 *
 */
static double step_geiger_method(
	double *lon0, double *lat0, double *depth0, double *time0, PICKS_POOL *pool,
	const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model, const double damping_matrix[HYPO_PARAMS_NUMBER]
) {
	int         i;
	int         valids = 0;
	DL_NODE    *node;
	PICK_STATE *pick;
/* */
	double result  = 0.0;
	double sum_wei = 0.0;
	double veli = 0.0;
	double velg = 0.0;
	double residual;
	double trv_time[pool->totals];
	double distance[pool->totals];
	double r_weight[pool->totals];
	double g_params[HYPO_PARAMS_NUMBER] = { 0.0 };
/* */
	MATRIX *matrix_g = NULL, *matrix_d = NULL, *matrix_w = NULL, *matrix_m = NULL;
/* */
	LINEAR_RAY_PATH_PARAMS ray_path;
/* */
	const double delta_x = el_misc_geog2distf( *lon0 - 0.5, *lat0, *lon0 + 0.5, *lat0 );
	const double delta_y = el_misc_geog2distf( *lon0, *lat0 - 0.5, *lon0, *lat0 + 0.5 );

/* Get the valids pick in the pool */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) )
			valids++;
	}
/* */
	INIT_RAY_PATH_PARAMS( ray_path );
	matrix_g = matrix_new( valids, HYPO_PARAMS_NUMBER );
	matrix_d = matrix_new( valids, 1 );
	matrix_w = matrix_new( valids, valids );
/* */
	i = 0;
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
		/* */
			SELECT_VEL_PHASE_DEPTH( veli, velg, *depth0, pick->observe.phase_name, p_model, s_model );
		/* */
			ray_path = get_linear_ray_path(
				*lon0, *lat0, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, *depth0, veli, velg
			);
		/* Assign derived values to picking */
			distance[i] = ray_path.epc_dist;
			trv_time[i] = get_linear_travel_time( &ray_path );
			residual    = pick->observe.picktime - (*time0 + trv_time[i]);
			r_weight[i] = get_r_weight( distance[i], *depth0, residual, pick->observe.weight, pick->flag );
		/* Get the derivatives of T */
			get_travel_time_derivatives( &ray_path, g_params );
		/* */
			if ( !damping_matrix ) {
				sum_wei += r_weight[i];
				result  += residual * r_weight[i];
			}
			else {
			/* Apply the damping values */
				for ( int j = 0; j < HYPO_PARAMS_NUMBER; j++ )
					g_params[j] *= damping_matrix[j];
			}
		/* Assign derived values to matrix */
			matrix_assign_row( matrix_g, g_params, i + 1, HYPO_PARAMS_NUMBER );
			i++;
		}
	}
/* Recalculate the travel time residual & weight */
	if ( !damping_matrix ) {
		result /= sum_wei;
		*time0 += result;
	}
/* */
	result = 0.0;
	i = 0;
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
		/* */
			residual = pick->observe.picktime - (*time0 + trv_time[i]);
			r_weight[i] = get_r_weight( distance[i], *depth0, residual, pick->observe.weight, pick->flag );
		/* */
			matrix_assign( matrix_d, residual, i + 1, 1 );
			result += r_weight[i];
			i++;
		}
	}
/* Construct the weighting matrix */
	result /= (double)valids;
	i = 0;
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
			matrix_assign( matrix_w, r_weight[i] / result, i + 1, i + 1);
			i++;
		}
	}

/* Go through the least square procedure & get the adjustments */
	matrix_m = matrix_div_weighted( matrix_d, matrix_g, matrix_w );
	matrix_extract_seq( matrix_m, g_params, HYPO_PARAMS_NUMBER );
	*lon0   += g_params[0] / delta_x;
	*lat0   += g_params[1] / delta_y;
	*depth0 += g_params[2];
	*time0  += g_params[3];
/* */
	for ( i = 0, result = 0.0; i < HYPO_PARAMS_NUMBER; i++ )
		result += fabs(g_params[i]);
/* Avoid the air quake & keep it always deeper than 5 km & less than 100 km */
	if ( *depth0 < 0.0 )
		*depth0 = fabs(*depth0);
	if ( *depth0 < MIN_HYPO_DEPTH )
		*depth0 = MIN_HYPO_DEPTH + EARLYLOC_EPSILON;
	else if ( *depth0 > MAX_HYPO_DEPTH )
		*depth0 = MAX_HYPO_DEPTH - EARLYLOC_EPSILON;
/* */
	matrix_free( matrix_g );
	matrix_free( matrix_d );
	matrix_free( matrix_w );
	matrix_free( matrix_m );

	return result;
}

/*
 *
 */
static double step_geiger_method_3D( double *lon0, double *lat0, double *depth0, double *time0, PICKS_POOL *pool, const double damping_matrix[HYPO_PARAMS_NUMBER] )
{
	int         i;
	int         valids = 0;
	DL_NODE    *node;
	PICK_STATE *pick;
/* */
	double result  = 0.0;
	double sum_wei = 0.0;
	double residual;
	double _x, _y;
	double trv_time[pool->totals];
	double distance[pool->totals];
	double r_weight[pool->totals];
	double g_params[HYPO_PARAMS_NUMBER] = { 0.0 };
/* */
	MATRIX *matrix_g = NULL, *matrix_d = NULL, *matrix_w = NULL, *matrix_m = NULL;
/* */
	RAY_INFO ray_path[RT_MAX_NODE + 1];
	int      np;
/* */
	const double delta_x = el_misc_geog2distf( *lon0 - 0.5, *lat0, *lon0 + 0.5, *lat0 );
	const double delta_y = el_misc_geog2distf( *lon0, *lat0 - 0.5, *lon0, *lat0 + 0.5 );

/* Get the valids pick in the pool */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) )
			valids++;
	}
/* */
	matrix_g = matrix_new( valids, HYPO_PARAMS_NUMBER );
	matrix_d = matrix_new( valids, 1 );
	matrix_w = matrix_new( valids, valids );
/* */
	i = 0;
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
		/* Do the ray tracing & get the travel time */
			if (
				rt_main(
					ray_path, &np, &trv_time[i], *lat0, *lon0, *depth0, pick->observe.latitude, pick->observe.longitude, pick->observe.elevation,
					!strcmp(pick->observe.phase_name, "S") ? RT_S_WAVE_VELOCITY : RT_P_WAVE_VELOCITY
				)
			) {
				return GEIGER_ERROR_RETURN;
			}
		/* */
			_x = (pick->observe.longitude - *lon0) * delta_x;
			_y = (pick->observe.latitude - *lat0) * delta_y;
		/* Assign derived values to picking */
			residual    = pick->observe.picktime - (*time0 + trv_time[i]);
			distance[i] = sqrt(_x * _x + _y * _y + EARLYLOC_EPSILON);
			r_weight[i] = get_r_weight( distance[i], *depth0, residual, pick->observe.weight, pick->flag );
			sum_wei    += r_weight[i];
			result     += residual * r_weight[i];
		/* Get the derivatives of T */
			get_travel_time_derivatives_3D( ray_path, np, g_params );
		/* Apply the damping values */
			if ( damping_matrix )
				for ( int j = 0; j < HYPO_PARAMS_NUMBER; j++ )
					g_params[j] *= damping_matrix[j];
		/* Assign derived values to matrix */
			matrix_assign_row( matrix_g, g_params, i + 1, HYPO_PARAMS_NUMBER );
			i++;
		}
	}
/* Recalculate the travel time residual & weight */
	result /= sum_wei;
	*time0 += result;
/* */
	result = 0.0;
	i = 0;
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
		/* */
			residual = pick->observe.picktime - (*time0 + trv_time[i]);
			r_weight[i] = get_r_weight( distance[i], *depth0, residual, pick->observe.weight, pick->flag );
		/* */
			matrix_assign( matrix_d, residual, i + 1, 1 );
			result += r_weight[i];
			i++;
		}
	}
/* Construct the weighting matrix */
	result /= (double)valids;
	i = 0;
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
			matrix_assign( matrix_w, r_weight[i] / result, i + 1, i + 1);
			i++;
		}
	}

/* Go through the least square procedure & get the adjustments */
	matrix_m = matrix_div_weighted( matrix_d, matrix_g, matrix_w );
	matrix_extract_seq( matrix_m, g_params, HYPO_PARAMS_NUMBER );
	*lon0   += g_params[0] / delta_x;
	*lat0   += g_params[1] / delta_y;
	*depth0 += g_params[2];
	*time0  += g_params[3];
/* */
	for ( i = 0, result = 0.0; i < HYPO_PARAMS_NUMBER; i++ )
		result += fabs(g_params[i]);
//#ifdef _DEBUG
	printf("%lf %lf %lf %lf\n", g_params[0], g_params[1], g_params[2], g_params[3]);
//#endif
/* Avoid the air quake & keep it always deeper than 5 km & less than 100 km */
	if ( *depth0 < 0.0 )
		*depth0 = fabs(*depth0);
	if ( *depth0 < MIN_HYPO_DEPTH )
		*depth0 = MIN_HYPO_DEPTH + EARLYLOC_EPSILON;
	else if ( *depth0 > MAX_HYPO_DEPTH )
		*depth0 = MAX_HYPO_DEPTH - EARLYLOC_EPSILON;
/* */
	matrix_free( matrix_g );
	matrix_free( matrix_d );
	matrix_free( matrix_w );
	matrix_free( matrix_m );

	return result;
}

/*
 *
 */
static void update_picks_state(
	const double lon0, const double lat0, const double depth0, const double time0, PICKS_POOL *pool,
	const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model
) {
	LINEAR_RAY_PATH_PARAMS ray_path;
	double      veli = 0.0;
	double      velg = 0.0;
	DL_NODE    *node;
	PICK_STATE *pick;
/* */
	const double delta_x = el_misc_geog2distf( lon0 - 0.5, lat0, lon0 + 0.5, lat0 );
	const double delta_y = el_misc_geog2distf( lon0, lat0 - 0.5, lon0, lat0 + 0.5 );

/* */
	INIT_RAY_PATH_PARAMS( ray_path );
/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		SELECT_VEL_PHASE_DEPTH( veli, velg, depth0, pick->observe.phase_name, p_model, s_model );
	/* */
		ray_path = get_linear_ray_path(
			lon0, lat0, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, depth0, veli, velg
		);
	/* */
		pick->distance = ray_path.epc_dist;
		pick->trv_time = get_linear_travel_time( &ray_path );
		pick->residual = pick->observe.picktime - (time0 + pick->trv_time);
		pick->r_weight = get_r_weight( pick->distance, depth0, pick->residual, pick->observe.weight, pick->flag );
	}

	return;
}

/*
 *
 */
static void update_picks_state_3D( const double lon0, const double lat0, const double depth0, const double time0, PICKS_POOL *pool )
{
	RAY_INFO    ray_path[RT_MAX_NODE + 1];
	int         np;
	double      _x;
	double      _y;
	DL_NODE    *node;
	PICK_STATE *pick;
/* */
	const double delta_x = el_misc_geog2distf( lon0 - 0.5, lat0, lon0 + 0.5, lat0 );
	const double delta_y = el_misc_geog2distf( lon0, lat0 - 0.5, lon0, lat0 + 0.5 );

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		rt_main(
			ray_path, &np, &pick->trv_time, lat0, lon0, depth0, pick->observe.latitude, pick->observe.longitude, pick->observe.elevation,
			!strcmp(pick->observe.phase_name, "S") ? RT_S_WAVE_VELOCITY : RT_P_WAVE_VELOCITY
		);
	/* */
		_x = (pick->observe.longitude - lon0) * delta_x;
		_y = (pick->observe.latitude - lat0) * delta_y;
	/* Assign derived values to picking */
		pick->distance = sqrt(_x * _x + _y * _y + EARLYLOC_EPSILON);
		pick->residual = pick->observe.picktime - (time0 + pick->trv_time);
		pick->r_weight = get_r_weight( pick->distance, depth0, pick->residual, pick->observe.weight, pick->flag );
	}

	return;
}

/*
 *
 */
static void update_hypo_state(
	const double lon0, const double lat0, const double depth0, const double time0, HYPO_STATE *hyp
) {
	DL_NODE    *node;
	PICK_STATE *pick;
/* */
	int    valids = 0;
	double error  = 0.0;
	double weight = 0.0;
	double min_epcdist = -1.0;

/* */
	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
			error  += pick->residual * pick->residual;
			weight += pick->r_weight * pick->r_weight;
		/* */
			if ( min_epcdist < 0.0 || pick->distance < min_epcdist )
				min_epcdist = pick->distance;
		/* */
			valids++;
		}
	}
/* */
	error  = sqrt(error / valids);
	weight = sqrt(weight / valids);

/* */
	hyp->latitude    = lat0;
	hyp->longitude   = lon0;
	hyp->depth       = depth0;
	hyp->origin_time = time0;
	hyp->gap         = get_gap_degree( lon0, lat0, &hyp->pool );
	hyp->avg_error   = error;
	hyp->avg_weight  = weight;
	hyp->q           = get_hypo_quality( valids, hyp->gap, hyp->depth, min_epcdist );

	return;
}

/*
 *
 */
static LINEAR_RAY_PATH_PARAMS get_linear_ray_path(
	const double epc_lon, const double epc_lat, const double sta_lon, const double sta_lat,
	const double delta_x, const double delta_y, const double hyp_depth, const double veli, const double velg
) {
	LINEAR_RAY_PATH_PARAMS result;

	result.epc_dist_x = (sta_lon - epc_lon) * delta_x;
	result.epc_dist_y = (sta_lat - epc_lat) * delta_y;
	result.epc_dist   = sqrt(result.epc_dist_x * result.epc_dist_x + result.epc_dist_y * result.epc_dist_y + EARLYLOC_EPSILON);
	result.hyp_depth  = hyp_depth;
	result.vel_init   = veli;
	result.vel_grad   = velg;

	result.center_z = -veli / velg;
	result.center_x = (result.epc_dist * result.epc_dist + 2.0 * result.center_z * hyp_depth - hyp_depth * hyp_depth) / (2. * result.epc_dist);
	result.angle_a  = atan((hyp_depth - result.center_z) / result.center_x);
	result.angle_b  = atan(-result.center_z / (result.epc_dist - result.center_x));
/* */
	if ( result.angle_a < 0.0 )
		result.angle_a += EARLYLOC_PI;
	result.angle_a = EARLYLOC_PI - result.angle_a;

	return result;
}

/*
 *
 */
static double get_linear_travel_time( const LINEAR_RAY_PATH_PARAMS *ray_path )
{
	return (-1.0 / ray_path->vel_grad) * log(fabs(tan(ray_path->angle_b * 0.5) / tan(ray_path->angle_a * 0.5)));
}

/*
 *
 */
static double *get_travel_time_derivatives(
	const LINEAR_RAY_PATH_PARAMS *ray_path, double derivatives[HYPO_PARAMS_NUMBER]
) {
	double tmp1 = ray_path->vel_init + ray_path->vel_grad * ray_path->hyp_depth;
	double tmp2 = -sin(ray_path->angle_a) / (tmp1 * ray_path->epc_dist);

/* Spatial derivative of T */
	derivatives[0] = tmp2 * ray_path->epc_dist_x;
	derivatives[1] = tmp2 * ray_path->epc_dist_y;
	derivatives[2] = -cos(ray_path->angle_a) / tmp1;
	derivatives[3] = 1.0;

	return derivatives;
}

/**
 * @brief Get the travel time derivatives in 3D velocity model
 *
 * @param ray_path
 * @param np
 * @param derivatives
 * @return double*
 */
static double *get_travel_time_derivatives_3D(
	const RAY_INFO *ray_path, const int np, double derivatives[HYPO_PARAMS_NUMBER]
) {
/* Spatial derivative of T */
	rt_drvt_cal( ray_path, np, &derivatives[0], &derivatives[1], &derivatives[2] );
/* The derivative for origin time always be 1 */
	derivatives[3] = 1.0;

	return derivatives;
}

/*
 *
 */
static double get_r_weight( const double epc_dist, const double hyp_depth, const double residual, const int pick_weight, const int pick_flag )
{
	double tmp;
	double result = 1.0;
/* */
	const double tres     = 1.0 + (pick_flag & PICK_FLAG_PRIMARY ? 1.0: 0.0) + (pick_weight < 2 ? 1.0 : 0.0);
	const double hyp_dist = sqrt(hyp_depth * hyp_depth + epc_dist * epc_dist + EARLYLOC_EPSILON);

/* */
	if ( hyp_dist > NEAR_HYPO_DISTANCE ) {
		result *=
			(FAR_HYPO_DISTANCE - NEAR_HYPO_DISTANCE) / (9. * hyp_dist + FAR_HYPO_DISTANCE - 10. * NEAR_HYPO_DISTANCE);
	}
/* */
	tmp     = tres / (tres + fabs(residual));
	result *= tmp * tmp;
	//result *= tmp * (tres / (tres + log10(epc_dist))) * (tres / (tres + pick_weight));

	return result;
}

/***/
static double *get_damping_matrix( double dmatrix[HYPO_PARAMS_NUMBER], const HYPO_STATE *hyp )
{
	DL_NODE    *node, *node_i;
	PICK_STATE *pick, *pick_i;
	double      t_diff, x_diff, y_diff;
/* */
	const double delta_x = el_misc_geog2distf( hyp->longitude - 0.5, hyp->latitude, hyp->longitude + 0.5, hyp->latitude );
	const double delta_y = el_misc_geog2distf( hyp->longitude, hyp->latitude - 0.5, hyp->longitude, hyp->latitude + 0.5 );

/* */
	dmatrix[0] = dmatrix[1] = dmatrix[2] = dmatrix[3] = EARLYLOC_EPSILON;
/* Go thru the pool */
	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( !(pick->flag & PICK_FLAG_REJECT) && !strcmp(pick->observe.phase_name, "P") ) {
			DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node_i, pick_i ) {
				if ( pick_i != pick && !(pick_i->flag & PICK_FLAG_REJECT) && !strcmp(pick_i->observe.phase_name, "P") ) {
					t_diff = pick->observe.picktime - pick_i->observe.picktime;
					t_diff *= t_diff;
					x_diff = fabs(pick->observe.longitude - pick_i->observe.longitude) * delta_x;
					y_diff = fabs(pick->observe.latitude - pick_i->observe.latitude) * delta_y;
				/* */
					dmatrix[0] += t_diff / x_diff;
					dmatrix[1] += t_diff / y_diff;
					dmatrix[2] += x_diff;
					dmatrix[3] += y_diff;
				}
			}
		}
	}
/* Get the ratio in t_diff */
	t_diff = (dmatrix[0] * dmatrix[3]) / (dmatrix[1] * dmatrix[2] + EARLYLOC_EPSILON);
	t_diff *= (dmatrix[3] * dmatrix[3]) / (dmatrix[2] * dmatrix[2] + EARLYLOC_EPSILON);
/* According to the ratio define the damping factors */
	dmatrix[0] = t_diff > 1.0 ? REFINE_DAMPING_FACTOR : REFINE_DAMPING_FACTOR * t_diff;
	dmatrix[1] = t_diff > 1.0 ? REFINE_DAMPING_FACTOR / t_diff : REFINE_DAMPING_FACTOR;
	dmatrix[2] = REFINE_MIN_DAMPING;
	dmatrix[3] = 0.0;
	if ( dmatrix[0] < REFINE_MIN_DAMPING )
		dmatrix[0] = REFINE_MIN_DAMPING;
	if ( dmatrix[1] < REFINE_MIN_DAMPING )
		dmatrix[1] = REFINE_MIN_DAMPING;

	return dmatrix;
}

/*
 *
 */
static double get_gap_degree( const double epc_lon, const double epc_lat, const PICKS_POOL *pool )
{
	int    i = 0;
	int    valids = 0;
	double result = 0.0;
	double ngap[pool->totals + 1];
	DL_NODE    *node;
	PICK_STATE *pick;

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) )
		/* */
			ngap[i++] = atan2(pick->observe.longitude - epc_lon, pick->observe.latitude - epc_lat);
	}
	valids = i;
/* */
	qsort(ngap, valids, sizeof(double), compare_gap);
	ngap[valids] = ngap[0] + EARLYLOC_PI2;
/* */
	for( i = 0; i < valids; i++ ) {
		if ( ngap[i + 1] - ngap[i] > result )
			result = ngap[i + 1] - ngap[i];
	}
	result *= EARLYLOC_RAD2DEG;

	return result;
}

/*
 *
 */
static int get_hypo_quality( const int valid_picks, const double gap, const double depth, const double min_epcdist )
{
	int result = HYPO_RESULT_QUALITY_D;

	if ( valid_picks >= 6 ) {
		if ( gap <= 90.0 && (min_epcdist <= depth || min_epcdist <= 5.0) ) {
			result = HYPO_RESULT_QUALITY_A;
		}
		else if ( gap <= 135.0 && (min_epcdist <= (depth * 2.0) || min_epcdist <= 10.0) ) {
			result = HYPO_RESULT_QUALITY_B;
		}
		else if ( gap <= 180.0 && min_epcdist <= 50.0 ) {
			result = HYPO_RESULT_QUALITY_C;
		}
	}

	return result;
}

/*
 *  compare_gap()
 */
static int compare_gap( const void *a, const void *b )
{
	if ( *(double *)a < *(double *)b )
		return -1;
	if ( *(double *)a > *(double *)b )
		return 1;

	return 0;
}
