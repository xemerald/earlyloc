/**
 * @file earlyloc_locate.h
 * @author Benjamin Yang in Department of Geology, National Taiwan University
 * @brief
 * @version 0.1
 * @date 2023-10-05
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
/* */
#include <earlyloc.h>

/* */
int el_loc_primary_locate( HYPO_STATE *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL * );
void el_loc_all_states_update( HYPO_STATE *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL * );
void el_loc_location_guess( HYPO_STATE *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL * );
void el_loc_location_refine( HYPO_STATE *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL *, const int );
double el_loc_origintime_adjust( HYPO_STATE *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL * );
double el_loc_residual_estimate( const HYPO_STATE *, const PICK_STATE *, const LAYER_VEL_MODEL *, const LAYER_VEL_MODEL * );
/* */
int el_loc_3dvelmod_load( const char * );
void el_loc_3dvelmod_free( void );
