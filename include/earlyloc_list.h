/*
 *
 */
#pragma once
/* */
#include <time.h>
/* */
#include <dbinfo.h>
#include <early_event_msg.h>
#include <earlyloc.h>
/* */
#define EARLYLOC_LIST_INITIALIZING  0
#define EARLYLOC_LIST_UPDATING      1
/* */
int      el_list_db_fetch( const char *, const DBINFO *, const int );
int      el_list_sta_line_parse( const char *, const int );
void     el_list_end( void );
USE_SNL *el_list_find( const EARLY_PICK_MSG * );
void     el_list_tree_activate( void );
void     el_list_tree_abandon( void );
int      el_list_total_station_get( void );
time_t   el_list_timestamp_get( void );
