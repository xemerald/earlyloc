/**
 * @file earlyloc_list.c
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-10-06
 *
 * @copyright Copyright (c) 2023
 *
 */
#define _GNU_SOURCE
/* Standard C header include */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <search.h>
#include <time.h>
#include <ctype.h>
/* Earthworm environment header include */
#include <earthworm.h>
#include <trace_buf.h>
/* Local header include */
#include <dblist.h>
#include <dbinfo.h>
#include <dl_chain_list.h>
#include <earlyloc.h>
#include <earlyloc_list.h>
/* */
typedef struct {
	int     count;      /* Number of clients in the list */
	time_t  timestamp;  /* Time of the last time updated */
	void   *entry;      /* Pointer to first client       */
	void   *root;       /* Root of binary searching tree */
	void   *root_t;     /* Temporary root of binary searching tree */
} SNLList;

/* */
static int      fetch_list_sql( SNLList *, const char *, const DBINFO *, const int );
static int      compare_snl( const void *, const void * );
static SNLList *init_snl_list( void );
static void     destroy_snl_list( SNLList * );
static void     dummy_func( void * );
static USE_SNL *update_stainfo( USE_SNL *, const USE_SNL * );
static USE_SNL *append_usesnl_list( SNLList *, USE_SNL *, const int );
static USE_SNL *create_new_usesnl(
	const char *, const char *, const char *, const double, const double, const double
);
static USE_SNL *enrich_stainfo_raw(
	USE_SNL *, const char *, const char *, const char *, const double, const double, const double
);
/* */
#if defined( _USE_SQL )
static void extract_stainfo_mysql(
	char *, char *, char *, double *, double *, double *, const MYSQL_ROW, const unsigned long []
);
#endif
/* */
static SNLList *SList = NULL;

/*
 * el_list_db_fetch() -
 */
int el_list_db_fetch( const char *table_sta, const DBINFO *dbinfo, const int update )
{
	if ( !SList ) {
		SList = init_snl_list();
		if ( !SList ) {
			logit("e", "earlyloc: Fatal! Trace list memory initialized error!\n");
			return -3;
		}
	}

	if ( strlen(dbinfo->host) > 0 && strlen(table_sta) > 0 )
		return fetch_list_sql( SList, table_sta, dbinfo, update );
	else
		return 0;
}

/*
 *
 */
int el_list_sta_line_parse( const char *line, const int update )
{
	int     result = 0;
	char    sta[TRACE2_STA_LEN] = { 0 };
	char    net[TRACE2_NET_LEN] = { 0 };
	char    loc[TRACE2_LOC_LEN] = { 0 };
	double  lat = 0.0;
	double  lon = 0.0;
	double  elv = 0.0;

/* */
	if ( !SList ) {
		SList = init_snl_list();
		if ( !SList ) {
			logit("e", "earlyloc: Fatal! Trace list memory initialized error!\n");
			return -3;
		}
	}
/* */
	if ( sscanf(line, "%s %s %s %lf %lf %lf", sta, net, loc, &lat, &lon, &elv) >= 6 ) {
	/* */
		if ( append_usesnl_list( SList, create_new_usesnl( sta, net, loc, lat, lon, elv ), update ) == NULL )
			result = -2;
	}
	else {
		logit("e", "earlyloc: ERROR, lack of some trace information in local list!\n");
		result = -1;
	}

	return result;
}

/*
 * el_list_end() -
 */
void el_list_end( void )
{
	destroy_snl_list( SList );
	SList = NULL;

	return;
}

/*
 * earlyloc_list_find() -
 */
USE_SNL *el_list_find( const EARLY_PICK_MSG *pick )
{
	USE_SNL *result = NULL;
	USE_SNL  key;

/* */
	memcpy(key.sta, pick->station, TRACE2_STA_LEN);
	memcpy(key.net, pick->network, TRACE2_NET_LEN);
	memcpy(key.loc, pick->location, TRACE2_LOC_LEN);
/* Find which station */
	if ( (result = tfind(&key, &SList->root, compare_snl)) != NULL ) {
	/* Found in the main Palert table */
		result = *(USE_SNL **)result;
	}

	return result;
}

/*
 * el_list_tree_activate() -
 */
void el_list_tree_activate( void )
{
	void *_root = SList->root;

	SList->root      = SList->root_t;
	SList->root_t    = NULL;
	SList->timestamp = time(NULL);

	if ( _root ) {
		sleep_ew(1000);
		tdestroy(_root, dummy_func);
	}

	return;
}

/*
 * el_list_tree_abandon() -
 */
void el_list_tree_abandon( void )
{
	if ( SList->root_t )
		tdestroy(SList->root_t, dummy_func);

	SList->root_t = NULL;

	return;
}

/*
 * el_list_total_station_get() -
 */
int el_list_total_station_get( void )
{
	DL_NODE *node   = NULL;
	int      result = 0;

/* */
	DL_LIST_FOR_EACH( (DL_NODE *)SList->entry, node ) {
		result++;
	}
/* */
	SList->count = result;

	return result;
}

/*
 * el_list_timestamp_get() -
 */
time_t el_list_timestamp_get( void )
{
	return SList->timestamp;
}

/*
 *
 */
#if defined( _USE_SQL )
/*
 * fetch_list_sql() - Get stations list from MySQL server
 */
static int fetch_list_sql( SNLList *list, const char *table_sta, const DBINFO *dbinfo, const int update )
{
	int     result = 0;
	char    sta[TRACE2_STA_LEN] = { 0 };
	char    net[TRACE2_NET_LEN] = { 0 };
	char    loc[TRACE2_LOC_LEN] = { 0 };
	double  lat;
	double  lon;
	double  elv;

	MYSQL_RES *sql_res = NULL;
	MYSQL_ROW  sql_row;

/* Connect to database */
	printf("earlyloc: Querying the stations information from MySQL server %s...\n", dbinfo->host);
	sql_res = dblist_sta_query_sql(
		dbinfo, table_sta, EARLYLOC_INFO_FROM_SQL,
		COL_STA_STATION, COL_STA_NETWORK, COL_STA_LOCATION,
		COL_STA_LATITUDE, COL_STA_LONGITUDE, COL_STA_ELEVATION
	);
	if ( sql_res == NULL )
		return -1;
	printf("earlyloc: Queried the stations information success!\n");

/* Start the SQL server connection for channel */
	dblist_start_persistent_sql( dbinfo );
/* Read station list from query result */
	while ( (sql_row = dblist_fetch_row_sql( sql_res )) != NULL ) {
	/* */
		extract_stainfo_mysql(
			sta, net, loc, &lat, &lon, &elv,
			sql_row, dblist_fetch_lengths_sql( sql_res )
		);
	/* */
		if (
			append_usesnl_list(
				list, create_new_usesnl( sta, net, loc, lat, lon, elv ), update
			) != NULL
		) {
			result++;
		}
		else {
			result = -2;
			break;
		}
	}
/* Close the connection for channel */
	dblist_close_persistent_sql();
	dblist_free_result_sql( sql_res );
	dblist_end_thread_sql();

	if ( result > 0 )
		logit("o", "earlyloc: Read %d stations information from MySQL server success!\n", result);
	else
		logit("e", "earlyloc: Some errors happened when fetching stations information from MySQL server!\n");

	return result;
}

/*
 * extract_stainfo_mysql() -
 */
static void extract_stainfo_mysql(
	char *sta, char *net, char *loc, double *lat, double *lon, double *elv,
	const MYSQL_ROW sql_row, const unsigned long row_lengths[]
) {
	char _str[32] = { 0 };

/* */
	dblist_field_extract_sql( sta, TRACE2_STA_LEN, sql_row[0], row_lengths[0] );
	dblist_field_extract_sql( net, TRACE2_NET_LEN, sql_row[1], row_lengths[1] );
	dblist_field_extract_sql( loc, TRACE2_LOC_LEN, sql_row[2], row_lengths[2] );
	*lat = atof(dblist_field_extract_sql( _str, sizeof(_str), sql_row[3], row_lengths[3] ));
	*lon = atof(dblist_field_extract_sql( _str, sizeof(_str), sql_row[4], row_lengths[4] ));
	*elv = atof(dblist_field_extract_sql( _str, sizeof(_str), sql_row[5], row_lengths[5] ));

	return;
}

#else
/*
 * fetch_list_sql() - Fake function
 */
static int fetch_list_sql( SNLList *list, const char *table_sta, const DBINFO *dbinfo, const int update )
{
	printf(
		"earlyloc: Skip the process of fetching station list from remote database "
		"'cause you did not define the _USE_SQL tag when compiling.\n"
	);
	return 0;
}
#endif

/*
 * init_snl_list() -
 */
static SNLList *init_snl_list( void )
{
	SNLList *result = (SNLList *)calloc(1, sizeof(SNLList));

	if ( result ) {
		result->count     = 0;
		result->timestamp = time(NULL);
		result->entry     = NULL;
		result->root      = NULL;
		result->root_t    = NULL;
	}

	return result;
}

/*
 * destroy_snl_list() -
 */
static void destroy_snl_list( SNLList *list )
{
	if ( list != (SNLList *)NULL ) {
	/* */
		tdestroy(list->root, dummy_func);
		dl_list_destroy( (DL_NODE **)&list->entry, free );
		free(list);
	}

	return;
}

/*
 *  append_usesnl_list() - Appending the new client to the client list.
 */
static USE_SNL *append_usesnl_list( SNLList *list, USE_SNL *usesnl, const int update )
{
	USE_SNL *result = NULL;
	void    **_root  = update == EARLYLOC_LIST_UPDATING ? &list->root : &list->root_t;

/* */
	if ( list && usesnl ) {
		if ( (result = tfind(usesnl, _root, compare_snl)) == NULL ) {
		/* Insert the station information into binary tree */
			if ( dl_node_append( (DL_NODE **)&list->entry, usesnl ) == NULL ) {
				logit("e", "earlyloc: Error insert channel into linked list!\n");
				goto except;
			}
			if ( (result = tsearch(usesnl, &list->root_t, compare_snl)) == NULL ) {
				logit("e", "earlyloc: Error insert channel into binary tree!\n");
				goto except;
			}
		}
		else if ( update == EARLYLOC_LIST_UPDATING ) {
			update_stainfo( *(USE_SNL **)result, usesnl );
			if ( (result = tsearch(usesnl, &list->root_t, compare_snl)) == NULL ) {
				logit("e", "earlyloc: Error insert channel into binary tree!\n");
				goto except;
			}
		}
		else {
			logit(
				"o", "earlyloc: SNL(%s.%s.%s) is already in the list, skip it!\n",
				usesnl->sta, usesnl->net, usesnl->loc
			);
			free(usesnl);
		}
	}

	return result ? *(USE_SNL **)result : NULL;
/* Exception handle */
except:
	free(usesnl);
	return NULL;
}

/*
 *  create_new_usesnl() - Creating new channel info memory space with the input value.
 */
static USE_SNL *create_new_usesnl(
	const char *sta, const char *net, const char *loc, const double lat, const double lon, const double elv
) {
	USE_SNL *result = (USE_SNL *)calloc(1, sizeof(USE_SNL));

/* */
	if ( result )
		enrich_stainfo_raw( result, sta, net, loc, lat, lon, elv );

	return result;
}

/*
 * enrich_stainfo_raw() -
 */
static USE_SNL *enrich_stainfo_raw(
	USE_SNL *usesnl, const char *sta, const char *net, const char *loc,
	const double lat, const double lon, const double elv
) {
/* */
	memcpy(usesnl->sta, sta, TRACE2_STA_LEN);
	memcpy(usesnl->net, net, TRACE2_NET_LEN);
	memcpy(usesnl->loc, loc, TRACE2_LOC_LEN);
	usesnl->latitude  = lat;
	usesnl->longitude = lon;
	usesnl->elevation = elv;

	return usesnl;
}

/*
 * update_stainfo() -
 */
static USE_SNL *update_stainfo( USE_SNL *dest, const USE_SNL *src )
{
/* */
	dest->latitude  = src->latitude;
	dest->longitude = src->longitude;
	dest->elevation = src->elevation;

	return dest;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return int
 */
static int compare_snl( const void *a, const void *b )
{
	int result;
	USE_SNL *tmpa = (USE_SNL *)a;
	USE_SNL *tmpb = (USE_SNL *)b;

/* */
	if ( (result = strcmp( tmpa->sta, tmpb->sta )) )
		return result;
	if ( (result = strcmp( tmpa->net, tmpb->net )) )
		return result;

	return strcmp( tmpa->loc, tmpb->loc );
}

/*
 * dummy_func() -
 */
static void dummy_func( void *node )
{
	return;
}
