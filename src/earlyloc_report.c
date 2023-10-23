/**
 * @file earlyloc_report.c
 * @author Benjamin Yang in Department of Geology, National Taiwan University
 * @brief
 * @version 0.1
 * @date 2023-10-05
 *
 * @copyright Copyright (c) 2023
 *
 */
/* */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
/* */
#include <earthworm.h>
#include <transport.h>
/* */
#include <dl_chain_list.h>
#include <early_event_msg.h>
#include <earlyloc.h>
#include <earlyloc_misc.h>

/* */
static char *gen_flag_mark( char *, const uint8_t );
static char *gen_timestamp_str( char *, const int, const double );
static void  get_picks_pool_time_ref( const PICKS_POOL *, double *, double * );
static void  send_result_message_ring( SHM_INFO *, MSG_LOGO *, const void *, const int );

/* */
#define REPORT_QUALITY_DESC_FORMAT  \
		"\nReporting time: %s, Avg. error: %.2f, Q: %c, Gap: %.0f, Avg. weight: %.2f, Total picks: %d(%d), Valid picks: %d, Valid mag.: %d, Padj: %.2f\n"

#define REPORT_HYPO_STATE_COLUMN  \
		" Year  Month  Day  Hour  Min.   Sec.     Lat.      Lon.   Depth   Mall  Mpd_s    Mpv    Mpd    Mtc  Proc_time  After_OT\n"

#define REPORT_HYPO_STATE_FORMAT  \
		" %4d     %02d   %02d    %02d    %02d  %05.2f  %7.4f  %8.4f  %6.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f     %6.2f    %6.2f\n"

#define REPORT_DIVIDER  \
		"==============================================================================================================================================================================\n"

#define REPORT_PICK_RESULT_COLUMN  \
		"  F     S   C  N  L      Lat.       Lon.         Pa        Pv        Pd        Tc   Mpv   Mpd   Mtc  Residual  Dist. Weight(P)       Arrival_time\n"

#define REPORT_PICK_RESULT_FORMAT  \
		" %2s %5s %s %s %s %9.5f %10.5f %10.6f %9.6f %9.6f %9.6f %5.2f %5.2f %5.2f  %8.4f  %5.1f   %.2f(%d)  %s\n"

#define REPORT_FLAG_MARK_CAPTION \
		"Flags: (x)Rejected, (-)Locate-masked, (o)Co-sited, (*)Primary, (+)In-used, (^)Magnitude-masked\n"

#define REPORT_RING_HYPO_OUTPUT_FORMAT  \
		"%d %f %f %d Mpd %.1f Mall %.1f %.2f %.2f %.1f %2d %2d %2d %.1f %.1f %d %d %.1f %s %.1f\n"

#define REPORT_RING_PICK_OUTPUT_FORMAT  \
		"%5s %s %s %s %9.5f %10.5f %9.6f %9.6f %9.6f %9.6f  %7.4f %5.1f %d %d %.1f %.2f %d %d %s\n"

#define REPORT_TIMESTAMP_BUFFER_SIZE  32
#define REPORT_TIMESTAMP_FORMAT      "%04d/%02d/%02d %02d:%02d:%05.2lf"

#define NEGLECT_DELAY_TIME            1.0f
/*
 *
 */
void el_report_file_output( const HYPO_STATE *hyp, const char *path, const char *postfix, const int output_reject )
{
	FILE  *fp = NULL;
	char   filename[MAX_PATH_STR] = { 0 };
	char   buffer[MAX_PATH_STR] = { 0 };
	char   flag_m[4] = { 0 };
	double tmp;
	double first_recv;
	double time_shift;
/* Nanosecond Timer */
	const double time_report = el_misc_timenow();

	time_t    time_tmp;
	struct tm brk_time;

	DL_NODE    *node;
	PICK_STATE *pick;
/* Quality strings */
#define X(a, b) b,
	const char quality[] = {
		HYPO_RESULT_QUALITY_TABLE
	};
#undef X

/* */
	get_picks_pool_time_ref( &hyp->pool, &first_recv, &time_shift );
/* */
	if ( postfix && strlen(postfix) ) {
		sprintf(
			filename, "%s_%d_%s_n%d.rep",
			el_misc_simple_timestamp_gen( buffer, MAX_PATH_STR, hyp->trigger_time ),
			hyp->eid, postfix, (int)hyp->rep_count
		);
	}
	else {
		sprintf(
			filename, "%s_%d_n%d.rep",
			el_misc_simple_timestamp_gen( buffer, MAX_PATH_STR, hyp->trigger_time ),
			hyp->eid, (int)hyp->rep_count
		);
	}
/* */
	sprintf(buffer, "%s%s", path, filename);
	fp = fopen(buffer, "w");
/* The first line: report information */
	fprintf(
		fp, REPORT_QUALITY_DESC_FORMAT,
		gen_timestamp_str( buffer, MAX_PATH_STR, time_report ),
		hyp->avg_error, quality[hyp->q], hyp->gap, hyp->avg_weight,
		hyp->pool.totals, hyp->pool.rejects, hyp->pool.valids,
		0, -1.0
	);
/* The second line: hypo result fields */
	fprintf(fp, REPORT_HYPO_STATE_COLUMN);
/* The third line: hypo result */
	time_tmp = (time_t)hyp->origin_time;
	gmtime_r(&time_tmp, &brk_time);
	tmp  = hyp->origin_time - (double)time_tmp;
	tmp += (double)brk_time.tm_sec;
	fprintf(
		fp, REPORT_HYPO_STATE_FORMAT,
		brk_time.tm_year + 1900, brk_time.tm_mon + 1, brk_time.tm_mday,
		brk_time.tm_hour, brk_time.tm_min, tmp,
		hyp->latitude, hyp->longitude, hyp->depth,
		-1.0, -1.0, -1.0, -1.0, -1.0,
		time_report - first_recv, time_report - (hyp->origin_time + time_shift)
	);
/* Divider only */
	fprintf(fp, REPORT_DIVIDER);
/* Following every pick */
	fprintf(fp, REPORT_PICK_RESULT_COLUMN);
	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( output_reject || !(pick->flag & PICK_FLAG_REJECT) ) {
		/* Get the hypo distance */
			tmp = sqrt(hyp->depth * hyp->depth + pick->distance * pick->distance);
			fprintf(
				fp, REPORT_PICK_RESULT_FORMAT,
				gen_flag_mark( flag_m, pick->flag ),
				pick->observe.station, pick->observe.channel, pick->observe.network, pick->observe.location,
				pick->observe.latitude, pick->observe.longitude,
				-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
				pick->residual, tmp, pick->r_weight, (int)pick->observe.weight,
				gen_timestamp_str( buffer, MAX_PATH_STR, pick->observe.picktime )
			);
		}
	}
/* Caption below */
	fprintf(fp, REPORT_DIVIDER);
	fprintf(fp, REPORT_FLAG_MARK_CAPTION);
/* Finish file */
	fclose(fp);
	logit("ot", "earlyloc: Hypo(#%d) finished report file %s!!\n", hyp->eid, filename);

	return;
}

/**
 * @brief
 *
 * @param hyp
 * @param region
 * @param logo
 * @param postfix
 * @param output_reject
 */
void el_report_ring_output(
	const HYPO_STATE *hyp, SHM_INFO *region, MSG_LOGO *logo, const char *postfix, const int output_reject
) {
#define X(a, b) b,
	const char quality[] = {
		HYPO_RESULT_QUALITY_TABLE
	};
#undef X
	char                    timestamp[REPORT_TIMESTAMP_BUFFER_SIZE] = { 0 };
	uint8_t                 buffer[EARLY_EVENT_SIZE];
	EARLY_EVENT_MSG_HEADER *eevh  = (EARLY_EVENT_MSG_HEADER *)buffer;
	EARLY_PICK_MSG         *epick = (EARLY_PICK_MSG   *)(eevh + 1);
	DL_NODE                *node;
	PICK_STATE             *pick;

/* */
	el_misc_simple_timestamp_gen( timestamp, REPORT_TIMESTAMP_BUFFER_SIZE, hyp->trigger_time );
/* */
	sprintf(eevh->event_id, "%s_%.16s", timestamp, postfix);
	eevh->seq = hyp->rep_count;
	eevh->origin_id = (long)hyp->origin_time;
	eevh->origin_time = hyp->origin_time;
	eevh->evlat = hyp->latitude;
	eevh->evlon = hyp->longitude;
	eevh->evdepth = hyp->depth;
	eevh->mag = -1.0;
	eevh->gap = hyp->gap;
	eevh->nsta = eevh->npha = hyp->pool.totals;
	eevh->suse = eevh->puse = hyp->pool.valids;
	eevh->dmin = 0.0;
	eevh->depth_flag = 'F';
	eevh->oterr = eevh->laterr = eevh->lonerr = eevh->deperr = -1.0;
	eevh->se = hyp->avg_error;
	eevh->errh = eevh->errz = eevh->avh = -1.0;
	eevh->q = quality[hyp->q];

	eevh->npicks = 0;
	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( output_reject || !(pick->flag & PICK_FLAG_REJECT) ) {
			pick->observe.dist = pick->distance;
			pick->observe.residual = pick->residual;
			*epick = pick->observe;
			epick->elevation *= -1000.0;
			eevh->npicks++;
			epick++;
		}
	}
/* */
	send_result_message_ring(region, logo, buffer, sizeof(EARLY_EVENT_MSG_HEADER) + eevh->npicks * sizeof(EARLY_PICK_MSG));

	return;
}

/**
 * @brief
 *
 * @param dest
 * @param flag
 * @return char*
 */
static char *gen_flag_mark( char *dest, const uint8_t flag )
{
/* */
	dest[1] = '\0';
/* */
	if ( flag & PICK_FLAG_REJECT ) {
		dest[0] = 'x';
	}
	else {
		if ( flag & PICK_FLAG_LOCMASK ) {
			dest[0] = '-';
		}
		else if ( flag & PICK_FLAG_COSITE ) {
			dest[0] = 'o';
		}
		else if ( flag & PICK_FLAG_PRIMARY ) {
			dest[0] = '*';
		}
		else if ( flag & PICK_FLAG_INUSE ) {
			dest[0] = '+';
		}
	/* */
		if ( flag & PICK_FLAG_MAGMASK ) {
			dest[1] = '^';
			dest[2] = '\0';
		}
	}

	return dest;
}

/**
 * @brief
 *
 * @param buffer
 * @param buffer_size
 * @param timestamp
 * @return char*
 */
static char *gen_timestamp_str( char *buffer, const int buffer_size, const double timestamp )
{
	struct tm    sptime;
	const time_t _timestamp = (time_t)timestamp;
	double sec_f = timestamp - (double)_timestamp;

/* Checking for buffer size */
	if ( buffer == NULL || buffer_size < REPORT_TIMESTAMP_BUFFER_SIZE )
		return buffer;
/* */
	gmtime_r(&_timestamp, &sptime);
	sec_f += (double)sptime.tm_sec;
	sprintf(
		buffer, REPORT_TIMESTAMP_FORMAT,
		sptime.tm_year + 1900, sptime.tm_mon + 1, sptime.tm_mday,
		sptime.tm_hour, sptime.tm_min, sec_f
	);

	return buffer;
}

/**
 * @brief Get the picks pool time ref object
 *
 * @param pool
 * @param first_recv
 * @param time_shift
 * @return double
 */
static void get_picks_pool_time_ref( const PICKS_POOL *pool, double *first_recv, double *time_shift )
{
	DL_NODE    *node;
	PICK_STATE *pick;
	int         npicks = 0;

/* */
	*first_recv = el_misc_timenow();
	*time_shift = 0.0;
/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( pick->recvtime < *first_recv )
			*first_recv = pick->recvtime;
	/* */
		*time_shift += pick->recvtime - pick->observe.picktime;
		npicks++;
	}
/* */
	*time_shift /= npicks;
	*time_shift -= NEGLECT_DELAY_TIME;
	if ( *time_shift < 0.0 )
		*time_shift = 0.0;

	return;
}

/**
 * @brief
 *
 * @param region
 * @param logo
 * @param outmsg
 * @param msg_size
 */
static void send_result_message_ring(
	SHM_INFO *region, MSG_LOGO *logo, const void *outmsg, const int msg_size
) {
	if ( tport_putmsg(region, logo, msg_size, (char *)outmsg) != PUT_OK )
		logit("et", "earlyloc: Error sending result to the output ring.\n");
	else
		logit("t", "earlyloc: Sending result message success!\n");

	return;
}
