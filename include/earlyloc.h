/**
 * @file earlyloc.h
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
#include <stdint.h>
#include <threads.h>
/* */
#include <earthworm.h>
#include <trace_buf.h>
/* */
#include <early_event_msg.h>
/* */
#define EARLYLOC_INFO_FROM_SQL  6
#define MAX_ALLOW_UPDATE_SEC    10
#define MAX_PATH_STR            512
#define MAX_POSTFIX_STR         32
/*
 *
 */
#define CONVERGE_CRITERIA      0.8f
#define ROUGH_ASSOC_CRITERIA   (CONVERGE_CRITERIA * 4.0f)
#define REJECT_CRITERIA        (CONVERGE_CRITERIA * 6.0f)
#define MIN_STABLE_SECOND      10
#define MIN_STABLE_REPCOUNT    MIN_STABLE_SECOND
#define MAX_ALLOW_REPCOUNT     512
#define PICK_COSITE_DIST       0.5f
#define MIN_HYPO_DEPTH         5.0f
#define MAX_HYPO_DEPTH         100.0f
#define MIN_DIFF_PICK_TIME     0.1f
#define RESIDUAL_EPSILON       1.0e-6
/*
 * Picking flag in Main & hypo pool
 */
#define PICK_FLAG_INUSE    0x01
#define PICK_FLAG_PRIMARY  0x02
#define PICK_FLAG_MAGMASK  0x04
#define PICK_FLAG_LOCMASK  0x08
#define PICK_FLAG_COSITE   0x10
/* ... */
#define PICK_FLAG_REJECT   0x80
/*
 * Hypo state flag
 */
#define HYPO_IS_UNUSED      0
#define HYPO_IS_WAITING     1
#define HYPO_IS_PROCESSING  2
#define HYPO_IS_FINISHED    3
/*
 * Hypo & Main pool flag
 */
#define POOL_ALREADY_USED  0
#define POOL_HAS_NEW_PICK  1
/*
 *
 */
#define HYPO_PARAMS_NUMBER  4
#define MIN_LOCATE_PICKS    HYPO_PARAMS_NUMBER
#define MIN_VALID_MAGS      MIN_LOCATE_PICKS
/*
 *
 */
#define PICK_FETCH_STRATEGY_TABLE \
		X(PICK_FETCH_STRATEGY_GREEDY, "greedy") \
		X(PICK_FETCH_STRATEGY_STEP,   "step"  ) \
		X(PICK_FETCH_STRATEGY_HYBRID, "hybrid") \
		X(PICK_FETCH_STRATEGY_COUNT,  "null"  )

#define X(a, b) a,
typedef enum {
	PICK_FETCH_STRATEGY_TABLE
} PICK_FETCH_STRATEGIES;
#undef X

/*
 *
 */
#define HYPO_RESULT_QUALITY_TABLE \
		X(HYPO_RESULT_QUALITY_A, 'A') \
		X(HYPO_RESULT_QUALITY_B, 'B') \
		X(HYPO_RESULT_QUALITY_C, 'C') \
		X(HYPO_RESULT_QUALITY_D, 'D')

#define X(a, b) a,
typedef enum {
	HYPO_RESULT_QUALITY_TABLE
} HYPO_RESULT_QUALITIES;
#undef X

/*
 *
 */
typedef struct {
	double boundary;
	double shallow_init;
	double shallow_grad;
	double deep_init;
	double deep_grad;
} LAYER_VEL_MODEL;


/**
 * @brief SNL info related struct
 *
 */
typedef struct {
	char sta[TRACE2_STA_LEN];   /* Site name (NULL-terminated) */
	char net[TRACE2_NET_LEN];   /* Network name (NULL-terminated) */
	char loc[TRACE2_LOC_LEN];   /* Location code (NULL-terminated) */

	double latitude;      /* Latitude of station in degree */
	double longitude;     /* Longitude of station in degree */
	double elevation;     /* Elevation of station in meter */
} USE_SNL;

/**
 * @brief
 *
 */
typedef struct {
/* Calculated & used in the locate process */
	uint8_t flag;
	double  recvtime;
	double  trv_time;
	double  distance;
	double  residual;
	double  r_weight;
/* */
	EARLY_PICK_MSG observe;
} PICK_STATE;

/**
 * @brief
 *
 */
typedef struct {
	DL_NODE *entry;
	DL_NODE *last;
	int      totals;
	int      valids;
	int      rejects;
	int      cosites;
} PICKS_POOL;

/*
 *
 */
typedef struct hypo_state {
/* */
	int      eid;
	thrd_t   tid;            /* Thread ID */
	uint8_t  flag;
	uint16_t rep_count;
/* */
	double latitude;
	double longitude;
	double depth;
	double origin_time;
	double trigger_time;
	double gap;
	double avg_error;
	double avg_weight;
	int    q;
/* Initial guess */
	double ig_latitude;
	double ig_longitude;
	double ig_depth;
	double ig_origin_time;
/* */
	PICKS_POOL pool;
	PICKS_POOL pick_queue;
	mtx_t      queue_mutex;
/* */
	struct hypo_state *best;
} HYPO_STATE;

/*
 *
 */
typedef struct {
	DL_NODE *entry;
	DL_NODE *last;
	int      totals;
} HYPOS_POOL;

/* */
#define EL_PICK_POOL_INIT(__PICK_POOL) \
		((__PICK_POOL) = (PICKS_POOL){ NULL, NULL, 0, 0, 0, 0 })
/* */
#define EL_HYPO_POOL_INIT(__HYPO_POOL) \
		((__HYPO_POOL) = (HYPOS_POOL){ NULL, NULL, 0 })
/* */
#define EL_PICK_VALID_LOCATE(__PICK) \
		(!((__PICK)->flag & PICK_FLAG_REJECT || (__PICK)->flag & PICK_FLAG_LOCMASK || (__PICK)->flag & PICK_FLAG_COSITE))
/* */
#define EL_MARK_PICK_FLAG(__FLAG, __MARK_FLAG) \
		((__FLAG) |= (__MARK_FLAG))
/* */
#define EL_UNMARK_PICK_FLAG(__FLAG, __UNMARK_FLAG) \
		((__FLAG) &= ~(__UNMARK_FLAG))
/* */
#define EL_MARK_UNMARK_PICK_FLAG(__FLAG, __MARK_FLAG, __UNMARK_FLAG) \
		((__FLAG) |= (__MARK_FLAG), (__FLAG) &= ~(__UNMARK_FLAG))
/* */
#define EL_MARK_PICK_REJECT(__PICK) \
		(EL_MARK_UNMARK_PICK_FLAG((__PICK)->flag, PICK_FLAG_REJECT, PICK_FLAG_LOCMASK))
/* */
#define EL_MARK_PICK_COSITE(__PICK) \
		(EL_MARK_UNMARK_PICK_FLAG((__PICK)->flag, PICK_FLAG_COSITE, PICK_FLAG_LOCMASK))
/* */
#define EL_MARK_PICK_LOCMASK(__PICK) \
		(EL_MARK_PICK_FLAG((__PICK)->flag, PICK_FLAG_LOCMASK))
/* */
#define EL_MARK_PICK_MAGMASK(__PICK) \
		(EL_MARK_PICK_FLAG((__PICK)->flag, PICK_FLAG_MAGMASK))
/* */
#define EL_MARK_PICK_PRIMARY(__PICK) \
		(EL_MARK_PICK_FLAG((__PICK)->flag, PICK_FLAG_PRIMARY))
/* */
#define EL_MARK_PICK_INUSE(__PICK) \
		(EL_MARK_PICK_FLAG((__PICK)->flag, PICK_FLAG_INUSE))
/* */
#define EL_PICK_LOCMASK_UNMARKABLE(__PICK) \
		((__PICK)->flag & PICK_FLAG_LOCMASK)
