/**
 * @file dl_chain_list.h
 * @author Benjamin Yang in Department of Geology, National Taiwan University
 * @brief
 * @version 0.1
 * @date 2019-03-05
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
/* */
typedef struct dl_node {
	void           *data;
	struct dl_node *prev;
	struct dl_node *next;
} DL_NODE;
/* External macros */
#define DL_NODE_GET_DATA(NODE)  ((NODE) ? (NODE)->data : NULL)
#define DL_NODE_GET_NEXT(NODE)  ((NODE) ? (NODE)->next : NULL)
#define DL_NODE_GET_PREV(NODE)  ((NODE) ? (NODE)->prev : NULL)
/* For loop for dl-list */
#define DL_LIST_FOR_EACH(HEAD, NODE) \
		for ( (NODE) = (HEAD); (NODE) != (DL_NODE *)NULL; (NODE) = (NODE)->next )

#define DL_LIST_FOR_EACH_DATA(HEAD, NODE, DATAP) \
		for ( (NODE) = (HEAD), (DATAP) = (__typeof__(DATAP))DL_NODE_GET_DATA(NODE); \
			(NODE) != (DL_NODE *)NULL; \
			(NODE) = (NODE)->next, (DATAP) = (__typeof__(DATAP))DL_NODE_GET_DATA(NODE) )

#define DL_LIST_FOR_EACH_SAFE(HEAD, NODE, SAFE) \
		for ( (NODE) = (HEAD), (SAFE) = DL_NODE_GET_NEXT(NODE); \
			(NODE) != (DL_NODE *)NULL; \
			(NODE) = (SAFE), (SAFE) = DL_NODE_GET_NEXT(NODE) )

#define DL_LIST_FOR_EACH_DATA_SAFE(HEAD, NODE, DATAP, SAFE) \
		for ( (NODE) = (HEAD), (SAFE) = DL_NODE_GET_NEXT(NODE), (DATAP) = (__typeof__(DATAP))DL_NODE_GET_DATA(NODE); \
			(NODE) != (DL_NODE *)NULL; \
			(NODE) = (SAFE), (SAFE) = DL_NODE_GET_NEXT(NODE), (DATAP) = (__typeof__(DATAP))DL_NODE_GET_DATA(NODE) )

/* Export functions' prototypes */
DL_NODE *dl_node_append( DL_NODE **, const void * );
DL_NODE *dl_node_insert( DL_NODE *, const void * );
DL_NODE *dl_node_push( DL_NODE **, const void * );
DL_NODE *dl_node_pop( DL_NODE ** );
DL_NODE *dl_node_delete( DL_NODE *, void (*)( void * ) );
void    *dl_node_data_extract( DL_NODE * );
void     dl_list_walk( DL_NODE *, void (*)( void *, const int, void * ), void * );
DL_NODE *dl_list_filter( DL_NODE **, int (*)( void *, void * ), void *, void (*)( void * ) );
void     dl_list_destroy( DL_NODE **, void (*)( void * ) );
