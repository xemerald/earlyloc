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
#define NEAR_HYPO_DISTANCE    80.0f
#define FAR_HYPO_DISTANCE     600.0f
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
/* Travel time */
	double traveltime;
} LINEAR_RAY_INFO;

/* */
static double step_geiger_method(
	double *, double *, double *, double *, PICKS_POOL *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL *, const double [HYPO_PARAMS_NUMBER]
);
static double step_geiger_method_3D( double *, double *, double *, double *, PICKS_POOL *, const double [HYPO_PARAMS_NUMBER] );
static double step_geiger_method_tdiff( double *, double *, double *, double *, PICKS_POOL *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL *, const int );
static void update_picks_state(
	const double, const double, const double, const double, PICKS_POOL *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL *
);
static void update_picks_state_3D( const double, const double, const double, const double, PICKS_POOL * );
static void update_hypo_state( const double, const double, const double, const double, HYPO_STATE * );
static LINEAR_RAY_INFO *get_linear_ray(
	LINEAR_RAY_INFO *, const double, const double, const double, const double, const double, const double, const double, const double, const double
);
static double *get_travel_time_derivatives( const LINEAR_RAY_INFO *, double [HYPO_PARAMS_NUMBER] );
static double *get_travel_time_derivatives_3D( const RAY_INFO *, const int, double [HYPO_PARAMS_NUMBER] );
static double  get_r_weight( const double, const double, const double, const int, const int );
static double  get_gap_degree( const double, const double, const PICKS_POOL * );
static double  get_pool_residual_avg( const double, const double, const double, const double, PICKS_POOL *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL * );
static int     get_hypo_quality( const int, const double, const double, const double );
static int     compare_gap( const void *, const void * );

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

#define INIT_LINEAR_RAY_INFO(__RAY_PATH) \
		((__RAY_PATH) = (LINEAR_RAY_INFO){ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 })

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

/* First, Geiger's method step for 20 times */
	for ( int i = 0; i < 20; i++ ) {
		if ( VelocityModel3DReady )
			tmp = step_geiger_method_3D( &lon0, &lat0, &depth0, &time0, &hyp->pool, NULL );
		else
			tmp = step_geiger_method( &lon0, &lat0, &depth0, &time0, &hyp->pool, p_model, s_model, NULL );
	/* */
		if ( tmp < 1.0 )
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
/* Calculate relative parameters one more by new hyp location */
	time0 += get_pool_residual_avg( lon0, lat0, depth0, time0, &hyp->pool, p_model, s_model );
	if ( VelocityModel3DReady )
		update_picks_state_3D( lon0, lat0, depth0, time0, &hyp->pool );
	else
		update_picks_state( lon0, lat0, depth0, time0, &hyp->pool, p_model, s_model );
/* Update all the parameters to the hypo state */
	update_hypo_state( lon0, lat0, depth0, time0, hyp );
#ifdef _DEBUG
	printf("earlyloc: Primary locate to: %lf %lf %lf %lf, average error: %lf\n", lon0, lat0, depth0, time0, hyp->avg_error );
#endif

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
void el_loc_location_refine( HYPO_STATE *hyp, const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model, const int use_weight )
{
	double lon0, lat0, depth0, time0;
	double lon1, lat1, depth1, time1;
	double werr0;

/* Update all the parameters to the hypo state */
	update_picks_state( hyp->longitude, hyp->latitude, hyp->depth, hyp->origin_time, &hyp->pool, p_model, s_model );
	update_hypo_state( hyp->longitude, hyp->latitude, hyp->depth, hyp->origin_time, hyp );
/* */
	lon1   = lon0   = hyp->longitude;
	lat1   = lat0   = hyp->latitude;
	depth1 = depth0 = hyp->depth;
	time1  = time0  = hyp->origin_time;
	werr0  = hyp->avg_error / hyp->avg_weight;
/* */
	depth1 = 10.0;
	for ( int i = 0; i < 20; i++ ) {
		if ( step_geiger_method_tdiff( &lon1, &lat1, &depth1, &hyp->origin_time, &hyp->pool, p_model, s_model, use_weight ) < 0.5 ) {
			break;
		}
		if ( fabs(hyp->longitude - lon1) > 2.0 || fabs(hyp->latitude - lat1) > 2.0 ) {
			lon1   = hyp->longitude;
			lat1   = hyp->latitude;
			depth1 = 10.0;
		}
	}
/* Update all the parameters to the hypo state */
	time1 += get_pool_residual_avg( lon1, lat1, depth1, time1, &hyp->pool, p_model, s_model );
	update_picks_state( lon1, lat1, depth1, time1, &hyp->pool, p_model, s_model );
	update_hypo_state( lon1, lat1, depth1, time1, hyp );
	if ( (hyp->avg_error / hyp->avg_weight) > werr0 ) {
		update_picks_state( lon0, lat0, depth0, time0, &hyp->pool, p_model, s_model );
		update_hypo_state( lon0, lat0, depth0, time0, hyp );
#ifdef _DEBUG
		printf("earlyloc: Back to: %lf %lf %lf %lf, average error: %lf\n", lon0, lat0, depth0, time0, hyp->avg_error );
	}
	else {
		printf("earlyloc: Refine: %lf %lf %lf %lf, average error: %lf\n", lon1, lat1, depth1, time1, hyp->avg_error );
#endif
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
	LINEAR_RAY_INFO ray_path;
/* */
	const double delta_x = el_misc_geog2distf( hyp->longitude - 0.5, hyp->latitude, hyp->longitude + 0.5, hyp->latitude );
	const double delta_y = el_misc_geog2distf( hyp->longitude, hyp->latitude - 0.5, hyp->longitude, hyp->latitude + 0.5 );

	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
		/* */
			SELECT_VEL_PHASE_DEPTH( veli, velg, hyp->depth, pick->observe.phase_name, p_model, s_model );
		/* */
			get_linear_ray(
				&ray_path, hyp->longitude, hyp->latitude, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, hyp->depth, veli, velg
			);
		/* Assign derived values to picking */
			residual = pick->observe.picktime - (hyp->origin_time + ray_path.traveltime);
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
	LINEAR_RAY_INFO ray_path;
	double veli = 0.0;
	double velg = 0.0;
/* */
	const double delta_x = el_misc_geog2distf( hyp->longitude - 0.5, hyp->latitude, hyp->longitude + 0.5, hyp->latitude );
	const double delta_y = el_misc_geog2distf( hyp->longitude, hyp->latitude - 0.5, hyp->longitude, hyp->latitude + 0.5 );

/* */
	INIT_LINEAR_RAY_INFO( ray_path );
	SELECT_VEL_PHASE_DEPTH( veli, velg, hyp->depth, pick->observe.phase_name, p_model, s_model );
/* */
	get_linear_ray(
		&ray_path, hyp->longitude, hyp->latitude, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, hyp->depth, veli, velg
	);
/* */

	return pick->observe.picktime - (hyp->origin_time + ray_path.traveltime);
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
	LINEAR_RAY_INFO ray_path;
/* */
	const double delta_x = el_misc_geog2distf( *lon0 - 0.5, *lat0, *lon0 + 0.5, *lat0 );
	const double delta_y = el_misc_geog2distf( *lon0, *lat0 - 0.5, *lon0, *lat0 + 0.5 );

/* Get the valids pick in the pool */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) )
			valids++;
	}
/* */
	INIT_LINEAR_RAY_INFO( ray_path );
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
			get_linear_ray(
				&ray_path, *lon0, *lat0, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, *depth0, veli, velg
			);
		/* Assign derived values to picking */
			distance[i] = ray_path.epc_dist;
			trv_time[i] = ray_path.traveltime;
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
	for ( i = 0; i < valids; i++ )
		matrix_assign( matrix_w, r_weight[i] / result, i + 1, i + 1);
/* Go through the least square procedure & get the adjustments */
	matrix_m = matrix_div_weighted( matrix_d, matrix_g, matrix_w );
	matrix_extract_seq( matrix_m, g_params, HYPO_PARAMS_NUMBER );
	*lon0   += g_params[0] / delta_x;
	*lat0   += g_params[1] / delta_y;
	*depth0 += g_params[2];
	*time0  += g_params[3];
/* */
	result = fabs(g_params[0]) + fabs(g_params[1]);
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
	for ( i = 0; i < valids; i++ )
		matrix_assign( matrix_w, r_weight[i] / result, i + 1, i + 1);
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
static double step_geiger_method_tdiff(
	double *lon0, double *lat0, double *depth0, double *time0, PICKS_POOL *pool, 
	const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model, const int use_weight
) {
	int         i;
	int         valids = 0;
	int         npairs = 0;
	DL_NODE    *node;
	PICK_STATE *pick;
/* */
	double result = 0.0;
	double veli = 0.0;
	double velg = 0.0;
	double residual[pool->totals];
	double r_weight[pool->totals];
	double drvts[pool->totals][HYPO_PARAMS_NUMBER];
	double g_params[HYPO_PARAMS_NUMBER];
/* */
	MATRIX *matrix_g = NULL, *matrix_d = NULL, *matrix_w = NULL, *matrix_m = NULL;
/* */
	LINEAR_RAY_INFO ray_path;
/* */
	const int    params_num = HYPO_PARAMS_NUMBER - 1;
	const double delta_x = el_misc_geog2distf( *lon0 - 0.5, *lat0, *lon0 + 0.5, *lat0 );
	const double delta_y = el_misc_geog2distf( *lon0, *lat0 - 0.5, *lon0, *lat0 + 0.5 );

/* Get the valids pick in the pool */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) )
			valids++;
	}
/* */
	for ( int j = 0; j < valids; j++ )
		for ( int k = j + 1; k < valids; k++ )
			npairs++;
/* */
	INIT_LINEAR_RAY_INFO( ray_path );
	matrix_g = matrix_new( npairs, params_num );
	matrix_d = matrix_new( npairs, 1 );
	if ( use_weight )
		matrix_w = matrix_new( npairs, npairs );
/* */
	i = 0;
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
		/* */
			SELECT_VEL_PHASE_DEPTH( veli, velg, *depth0, pick->observe.phase_name, p_model, s_model );
		/* */
			get_linear_ray(
				&ray_path, *lon0, *lat0, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, *depth0, veli, velg
			);
		/* Assign derived values to picking */
			residual[i] = pick->observe.picktime - (*time0 + ray_path.traveltime);
			if ( use_weight ) {
				r_weight[i] = get_r_weight( ray_path.epc_dist, *depth0, residual[i], pick->observe.weight, pick->flag );
				result     += r_weight[i];
			}
		/* Get the derivatives of T */
			get_travel_time_derivatives( &ray_path, drvts[i] );
			i++;
		}
	}
/* */
	if ( use_weight ) {
		result /= valids;
		result *= result;
	}
	i = 0;
	for ( int j = 0; j < valids; j++ ) {
		for ( int k = j + 1; k < valids; k++ ) {
			g_params[0] = drvts[j][0] - drvts[k][0];
			g_params[1] = drvts[j][1] - drvts[k][1];
			g_params[2] = drvts[j][2] - drvts[k][2];
		/* Assign derived values to matrix */
			matrix_assign_row( matrix_g, g_params, i + 1, params_num );
			matrix_assign( matrix_d, residual[j] - residual[k], i + 1, 1 );
			if ( use_weight )
				matrix_assign( matrix_w, (r_weight[j] * r_weight[k]) / result, i + 1, i + 1 );
			i++;
		}
	}
/* Go through the least square procedure & get the adjustments */
	if ( use_weight )
		matrix_m = matrix_div_weighted( matrix_d, matrix_g, matrix_w );
	else
		matrix_m = matrix_div( matrix_d, matrix_g );
	matrix_extract_seq( matrix_m, g_params, params_num );
	*lon0   += g_params[0] / delta_x;
	*lat0   += g_params[1] / delta_y;
	*depth0 += g_params[2];
#ifdef _DEBUG
	printf("Adjust: %lf, %lf, %lf\n", g_params[0], g_params[1], g_params[2]);
#endif
/* */
	result = fabs(g_params[0]) + fabs(g_params[1]);
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
	matrix_free( matrix_m );
	if ( use_weight )
		matrix_free( matrix_w );

	return result;
}

/**
 * @brief
 *
 * @param lon0
 * @param lat0
 * @param depth0
 * @param time0
 * @param pool
 * @param p_model
 * @param s_model
 */
static void update_picks_state(
	const double lon0, const double lat0, const double depth0, const double time0, PICKS_POOL *pool,
	const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model
) {
	LINEAR_RAY_INFO ray_path;
	double      veli = 0.0;
	double      velg = 0.0;
	DL_NODE    *node;
	PICK_STATE *pick;
/* */
	const double delta_x = el_misc_geog2distf( lon0 - 0.5, lat0, lon0 + 0.5, lat0 );
	const double delta_y = el_misc_geog2distf( lon0, lat0 - 0.5, lon0, lat0 + 0.5 );

/* */
	INIT_LINEAR_RAY_INFO( ray_path );
/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		SELECT_VEL_PHASE_DEPTH( veli, velg, depth0, pick->observe.phase_name, p_model, s_model );
	/* */
		get_linear_ray(
			&ray_path, lon0, lat0, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, depth0, veli, velg
		);
	/* */
		pick->distance = ray_path.epc_dist;
		pick->trv_time = ray_path.traveltime;
		pick->residual = pick->observe.picktime - (time0 + pick->trv_time);
		pick->r_weight = get_r_weight( pick->distance, depth0, pick->residual, pick->observe.weight, pick->flag );
	}

	return;
}

/**
 * @brief
 *
 * @param lon0
 * @param lat0
 * @param depth0
 * @param time0
 * @param pool
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

/**
 * @brief
 *
 * @param lon0
 * @param lat0
 * @param depth0
 * @param time0
 * @param hyp
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

/**
 * @brief Get the linear ray object
 *
 * @param dest
 * @param epc_lon
 * @param epc_lat
 * @param sta_lon
 * @param sta_lat
 * @param delta_x
 * @param delta_y
 * @param hyp_depth
 * @param veli
 * @param velg
 * @return LINEAR_RAY_INFO*
 */
static LINEAR_RAY_INFO *get_linear_ray(
	LINEAR_RAY_INFO *dest, const double epc_lon, const double epc_lat, const double sta_lon, const double sta_lat,
	const double delta_x, const double delta_y, const double hyp_depth, const double veli, const double velg
) {
	dest->epc_dist_x = (sta_lon - epc_lon) * delta_x;
	dest->epc_dist_y = (sta_lat - epc_lat) * delta_y;
	dest->epc_dist   = sqrt(dest->epc_dist_x * dest->epc_dist_x + dest->epc_dist_y * dest->epc_dist_y + EARLYLOC_EPSILON);
	dest->hyp_depth  = hyp_depth;
	dest->vel_init   = veli;
	dest->vel_grad   = velg;

	dest->center_z = -veli / velg;
	dest->center_x = (dest->epc_dist * dest->epc_dist + 2.0 * dest->center_z * hyp_depth - hyp_depth * hyp_depth) / (2. * dest->epc_dist);
	dest->angle_a  = atan((hyp_depth - dest->center_z) / dest->center_x);
	dest->angle_b  = atan(-dest->center_z / (dest->epc_dist - dest->center_x));
/* */
	if ( dest->angle_a < 0.0 )
		dest->angle_a += EARLYLOC_PI;
	dest->angle_a = EARLYLOC_PI - dest->angle_a;
/* */
	dest->traveltime = (-1.0 / dest->vel_grad) * log(fabs(tan(dest->angle_b * 0.5) / tan(dest->angle_a * 0.5)));

	return dest;
}

/**
 * @brief Get the travel time derivatives object
 *
 * @param ray_path
 * @param derivatives
 * @return double*
 */
static double *get_travel_time_derivatives( const LINEAR_RAY_INFO *ray_path, double derivatives[HYPO_PARAMS_NUMBER] )
{
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

/**
 * @brief Get the r weight object
 *
 * @param epc_dist
 * @param hyp_depth
 * @param residual
 * @param pick_weight
 * @param pick_flag
 * @return double
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

/**
 * @brief Get the gap degree object
 *
 * @param epc_lon
 * @param epc_lat
 * @param pool
 * @return double
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

/**
 * @brief Get the hypo quality object
 *
 * @param valid_picks
 * @param gap
 * @param depth
 * @param min_epcdist
 * @return int
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

/***/
static double get_pool_residual_avg( const double lon0, const double lat0, const double depth0, const double time0, PICKS_POOL *pool, const LAYER_VEL_MODEL *p_model, const LAYER_VEL_MODEL *s_model )
{
	DL_NODE    *node;
	PICK_STATE *pick;

	int    valids = 0;
	double result = 0.0;
	double veli   = 0.0;
	double velg   = 0.0;
	LINEAR_RAY_INFO ray_path;
/* */
	const double delta_x = el_misc_geog2distf( lon0 - 0.5, lat0, lon0 + 0.5, lat0 );
	const double delta_y = el_misc_geog2distf( lon0, lat0 - 0.5, lon0, lat0 + 0.5 );

	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( EL_PICK_VALID_LOCATE( pick ) ) {
		/* */
			SELECT_VEL_PHASE_DEPTH( veli, velg, depth0, pick->observe.phase_name, p_model, s_model );
		/* */
			get_linear_ray(
				&ray_path, lon0, lat0, pick->observe.longitude, pick->observe.latitude, delta_x, delta_y, depth0, veli, velg
			);
		/* */
			result += pick->observe.picktime - (time0 + ray_path.traveltime);
			valids++;
		}
	}
/* Recalculate the travel time residual & weight */
	result /= valids;

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
