#ifdef _OS2
#define INCL_DOSMEMMGR
#define INCL_DOSSEMAPHORES
#include <os2.h>
#endif
/* Standard C header include */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <threads.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/stat.h>
/* Earthworm environment header include */
#include <earthworm.h>
#include <kom.h>
#include <transport.h>
#include <lockfile.h>
#include <trace_buf.h>
/* Local header include */
#include <dbinfo.h>
#include <early_event_msg.h>
#include <dl_chain_list.h>
#include <earlyloc.h>
#include <earlyloc_list.h>
#include <earlyloc_misc.h>
#include <earlyloc_locate.h>
#include <earlyloc_report.h>

/* Functions prototype in this source file */
static void earlyloc_config( char * );
static void earlyloc_lookup( void );
static void earlyloc_status( uint8_t, short, char * );
static void earlyloc_end( void );                /* Free all the local memory & close socket */


static HYPO_STATE *save_to_best_state( HYPO_STATE * );
static HYPO_STATE *use_init_guess( HYPO_STATE * );
static HYPO_STATE *save_to_init_guess( HYPO_STATE * );
static HYPO_STATE *reset_init_guess( HYPO_STATE * );

static EARLY_PICK_MSG *fill_coor2epick( EARLY_PICK_MSG *, const USE_SNL * );
static PICK_STATE *parse_epick2pickstate( PICK_STATE *, const EARLY_PICK_MSG * );
static int mk_outdir_by_evt( char *, const char *, const double, const int, const char * );

static int thread_proc_trigger( void * );

static HYPO_STATE *create_hypo_state( void );
static void associate_pick_hypos( HYPOS_POOL *, PICK_STATE * );
static HYPO_STATE *insert_hypo_to_pool( HYPOS_POOL *, HYPO_STATE * );
static void bootstrap_hypo_in_pool( HYPOS_POOL * );
static void close_obsolete_hypos( HYPOS_POOL * );
static HYPOS_POOL *destroy_hypo_pool( HYPOS_POOL * );
static PICKS_POOL *destroy_pick_pool( PICKS_POOL * );
static int check_pick_exist_pool( const PICKS_POOL *, const PICK_STATE * );
static int check_scnlp_exist_pool( const PICKS_POOL *, const PICK_STATE * );
static int check_pick_cluster_pool( const PICKS_POOL *, const PICK_STATE * );
static int check_pick_cosite_pool( const PICKS_POOL *, const PICK_STATE * );
static PICK_STATE *insert_pick_to_pool( PICKS_POOL *, const PICK_STATE *, const int );
static PICK_STATE *pop_pick_from_pool( PICKS_POOL * );
static PICK_STATE *mark_locmask_pick_max_residual( PICKS_POOL * );
static void remove_obsolete_picks( PICKS_POOL * );
static void copy_associated_picks( PICKS_POOL *, const PICKS_POOL *, const int );
static void mark_reject_picks( PICKS_POOL *, double );
static void unmark_locmask_picks( PICKS_POOL *, double );
static void mark_primary_picks( PICKS_POOL * );
static int gen_hypo_eid( const HYPO_STATE * );
static int cluster_pick( const PICK_STATE *, const PICK_STATE * );
static int compare_pick( const void *, const void * );
static int compare_scnl( const void *, const void * );

/* Ring messages things */
static  SHM_INFO  InRegion;      /* shared memory region to use for i/o    */
static  SHM_INFO  OutRegion;     /* shared memory region to use for i/o    */

#define MAXLOGO 3

MSG_LOGO  Getlogo[MAXLOGO];       /* array for requesting module, type, instid */
pid_t     myPid;                  /* for restarts by startstop                 */

#define MAXLIST  5
/* Things to read or derive from configuration file */
static char     InRingName[MAX_RING_STR];   /* name of transport ring for i/o    */
static char     OutRingName[MAX_RING_STR];  /* name of transport ring for i/o    */
static char     MyModName[MAX_MOD_STR];     /* speak as this module name/id      */
static uint8_t  LogSwitch;                  /* 0 if no logfile should be written */
static uint64_t HeartBeatInterval;          /* seconds between heartbeats        */
static uint16_t nLogo = 0;
static char     ReportPath[MAX_PATH_STR] = { 0 };        /*  */
static char     OutputPostfix[MAX_POSTFIX_STR] = { 0 };
static uint8_t  OutputRejectPick = 0;         /* 0 if don't want to output rejected pick */
static uint8_t  OutputPickRing = 0;           /* 0 if don't want to output pick to ring  */
static uint8_t  PickFetchStrategy = PICK_FETCH_STRATEGY_GREEDY;
static uint16_t ReportTermNum = 0;
static uint16_t IgnoreWeightP = 5;
static uint16_t IgnoreWeightS = 5;
static double   PickAliveTime = 60.0;    /* */
static uint16_t TriggerPicks  = 5;
static double   HypoAliveTime = 60.0;    /* */
static double   ClusterTimeDiff;         /* */
static double   ClusterDist;             /* */
static uint16_t ClusterPicks = 2;
static LAYER_VEL_MODEL PWaveModel;
static LAYER_VEL_MODEL SWaveModel;
static DBINFO   DBInfo;
static char     SQLStationTable[MAXLIST][MAX_TABLE_LEGTH];
static uint16_t nList = 0;

/* Things to look up in the earthworm.h tables with getutil.c functions */
static int64_t InRingKey;       /* key of transport ring for i/o     */
static int64_t OutRingKey;      /* key of transport ring for i/o     */
static uint8_t InstId;          /* local installation id             */
static uint8_t MyModId;         /* Module Id for this program        */
static uint8_t TypeHeartBeat;
static uint8_t TypeError;
static uint8_t TypePickInput;
static uint8_t TypeHypoOutput;

/* Error messages used by earlyloc */
#define  ERR_MISSMSG       0   /* message missed in transport ring       */
#define  ERR_TOOBIG        1   /* retreived msg too large for buffer     */
#define  ERR_NOTRACK       2   /* msg retreived; tracking limit exceeded */
#define  ERR_QUEUE         3   /* error queueing message for sending      */
static char Text[150];         /* string for log/error messages          */

/* Station list update status flag
 *********************************/
#define  LIST_IS_UPDATED      0
#define  LIST_NEED_UPDATED    1
#define  LIST_UNDER_UPDATE    2

static volatile uint8_t UpdateStatus = LIST_IS_UPDATED;
/* Macros */
#define HYPO_IS_CONVERGED(__HYPO) \
		((__HYPO)->avg_error <= CONVERGE_CRITERIA && (__HYPO)->avg_error > 0.1 && (__HYPO)->avg_weight > 0.1)

#define PICK_IS_ON_SURFACE(__PICK) \
		(!(memcmp((__PICK)->observe.location, "--", 2)) || \
		!(memcmp((__PICK)->observe.location, "10", 2)) || \
		!(memcmp((__PICK)->observe.location, "01", 2)) || \
		!(memcmp((__PICK)->observe.location, "  ", 2)) || \
		(__PICK)->observe.location[0] == '\0')

#define PICK_IS_ASSOCIATED_WITH_POOL(__POOL, __PICK) \
		(check_pick_cluster_pool( (__POOL), (__PICK) ) || \
		check_pick_cosite_pool( (__POOL), (__PICK) ))

#define PICKS_ARE_CLUSTER(PICK_A, PICK_B) \
		(compare_scnl( (PICK_A), (PICK_B) ) && cluster_pick( (PICK_A), (PICK_B) ))

#define PICKS_ARE_COSITE(PICK_A, PICK_B) \
		(compare_scnl( (PICK_A), (PICK_B) ) && \
		el_misc_geog2distf( (PICK_A)->observe.longitude, (PICK_A)->observe.latitude, (PICK_B)->observe.longitude, (PICK_B)->observe.latitude ) <= PICK_COSITE_DIST)

/*
 *
 */
int main ( int argc, char **argv )
{
	int res;
	int i;

	int64_t  recsize = 0;
	MSG_LOGO reclogo;
	time_t   time_now;           /* current time                  */
	time_t   time_last_beat;     /* time last heartbeat was sent  */
	time_t   time_last_scan;     /* time last heartbeat was sent  */
	char    *lockfile;
	int32_t  lockfile_fd;
/* */
	EARLY_PICK_MSG input_pick;
/* */
	USE_SNL    *usesnl = NULL;
	PICK_STATE  pick_state;
	HYPO_STATE *hypo_state;
	PICKS_POOL  main_pool;
	HYPOS_POOL  hypo_pool;

/* Check command line arguments */
	if ( argc != 2 ) {
		fprintf(stderr, "Usage: earlyloc <configfile>\n");
		exit(0);
	}
/* Initialize name of log-file & open it */
	logit_init(argv[1], 0, 256, 1);
/* Read the configuration file(s) */
	earlyloc_config( argv[1] );
	logit("" , "%s: Read command file <%s>\n", argv[0], argv[1]);
	/* Read the channels list from remote database */
	for ( i = 0; i < nList; i++ ) {
		if ( el_list_db_fetch( SQLStationTable[i], &DBInfo, EARLYLOC_LIST_INITIALIZING ) < 0 ) {
			fprintf(stderr, "Something error when fetching stations list %s. Exiting!\n", SQLStationTable[i]);
			exit(-1);
		}
	}
/* Checking total station number again */
	if ( !(i = el_list_total_station_get()) ) {
		fprintf(stderr, "There is not any station in the list after fetching. Exiting!\n");
		exit(-1);
	}
	else {
		logit("o", "earlyloc: There are total %d station(s) in the list.\n", i);
		el_list_tree_activate();
	}
/* Look up important info from earthworm.h tables */
	earlyloc_lookup();
/* Reinitialize logit to desired logging level */
	logit_init(argv[1], 0, 256, LogSwitch);
	lockfile = ew_lockfile_path(argv[1]);
	if ( (lockfile_fd = ew_lockfile(lockfile) ) == -1 ) {
		fprintf(stderr, "one instance of %s is already running, exiting\n", argv[0]);
		exit(-1);
	}
/* Get process ID for heartbeat messages */
	myPid = getpid();
	if ( myPid == -1 ) {
		logit("e","earlyloc: Cannot get pid. Exiting.\n");
		exit(-1);
	}

/* Attach to Input/Output shared memory ring */
	tport_attach(&InRegion, InRingKey);
	logit("", "earlyloc: Attached to public memory region %s: %ld\n", InRingName, InRingKey);
/* Flush the transport ring */
	tport_flush(&InRegion, Getlogo, nLogo, &reclogo);
	tport_attach(&OutRegion, OutRingKey);
	logit("", "earlyloc: Attached to public memory region %s: %ld\n", OutRingName, OutRingKey);

/* */
	EL_HYPO_POOL_INIT( hypo_pool );
/* */
	EL_PICK_POOL_INIT( main_pool );

/* Force a heartbeat to be issued in first pass thru main loop */
	time_last_beat = time(&time_now) - HeartBeatInterval - 1;
	time_last_scan = time_last_beat + 1;
/*----------------------- setup done; start main loop -------------------------*/
	while ( 1 ) {
	/* Send earlyloc's heartbeat */
		if ( time(&time_now) - time_last_beat >= (int64_t)HeartBeatInterval ) {
			time_last_beat = time_now;
			earlyloc_status( TypeHeartBeat, 0, "" );
		}
	/*  */
		if ( time_now - time_last_scan >= 1 ) {
			time_last_scan = time_now;
			remove_obsolete_picks( &main_pool );
			close_obsolete_hypos( &hypo_pool );
		}

	/* Process all new messages */
		do {
		/* See if a termination has been requested */
			i = tport_getflag(&InRegion);
			if ( i == TERMINATE || i == myPid ) {
			/* write a termination msg to log file */
				logit("t", "earlyloc: Termination requested; exiting!\n");
				fflush(stdout);
				goto exit_procedure;
			}

		/* Get msg & check the return code from transport */
			res = tport_getmsg(&InRegion, Getlogo, nLogo, &reclogo, &recsize, (char *)&input_pick, sizeof(EARLY_PICK_MSG));
		/* no more new messages */
			if ( res == GET_NONE ) {
				break;
			}
		/* next message was too big */
			else if ( res == GET_TOOBIG ) {
			/* complain and try again */
				sprintf(
					Text, "Retrieved msg[%ld] (i%u m%u t%u) too big for Buffer[%ld]",
					recsize, reclogo.instid, reclogo.mod, reclogo.type, sizeof(EARLY_PICK_MSG)
				);
				earlyloc_status( TypeError, ERR_TOOBIG, Text );
				continue;
			}
		/* got a msg, but missed some */
			else if ( res == GET_MISS ) {
				sprintf(
					Text, "Missed msg(s)  i%u m%u t%u  %s.",
					reclogo.instid, reclogo.mod, reclogo.type, InRingName
				);
				earlyloc_status( TypeError, ERR_MISSMSG, Text );
			}
		/* got a msg, but can't tell */
			else if ( res == GET_NOTRACK ) {
			/* if any were missed */
				sprintf(
					Text, "Msg received (i%u m%u t%u); transport.h NTRACK_GET exceeded",
					reclogo.instid, reclogo.mod, reclogo.type
				);
				earlyloc_status( TypeError, ERR_NOTRACK, Text );
			}

		/* Process the message */
			if ( reclogo.type == TypePickInput ) {
			#ifdef _DEBUG
				printf("earlyloc: Got a new pick from %s.%s.%s.%s!\n",
				input_pick.station, input_pick.channel, input_pick.network, input_pick.location);
			#endif
				if ( !(usesnl = el_list_find( &input_pick )) ) {
				#ifdef _DEBUG
				/* Not found in trace table */
					printf("earlyloc: Pick with %s.%s.%s.%s not found in SNL table, maybe it's a new SNL.\n",
					input_pick.station, input_pick.channel, input_pick.network, input_pick.location);
				#endif
				/* Force to update the table */
					if ( UpdateStatus == LIST_IS_UPDATED )
						UpdateStatus = LIST_NEED_UPDATED;
					continue;
				}
			/* Make the global pick information to local pick type */
				parse_epick2pickstate( &pick_state, fill_coor2epick( &input_pick, usesnl ) );
			/* Drop the existed pick */
				if ( check_pick_exist_pool( &main_pool, &pick_state ) )
					continue;
			/* Mark the pick with higher weight (w\ low SNR) */
				if (
					(!strcmp(pick_state.observe.phase_name, "P") && pick_state.observe.weight >= IgnoreWeightP) ||
					(!strcmp(pick_state.observe.phase_name, "S") && pick_state.observe.weight >= IgnoreWeightS)
				) {
					EL_MARK_PICK_REJECT( &pick_state );
				}
			/* Insert the new picking to the existed hypos' pick pool */
				associate_pick_hypos( &hypo_pool, &pick_state );
			/*
			 * If this new picking was not been inserted to any hypo pool,
			 * check it could cluster with pickings inside main pool
			 */
				if ( !(pick_state.flag & PICK_FLAG_INUSE) && !(pick_state.flag & PICK_FLAG_REJECT) && main_pool.totals > ClusterPicks ) {
					if ( check_pick_cluster_pool( &main_pool, &pick_state ) ) {
					/* Creating new hypo state to process trigger... */
						if ( (hypo_state = create_hypo_state()) ) {
						/* Mark the pick as in-use & primary */
							EL_MARK_PICK_INUSE( &pick_state );
							insert_pick_to_pool( &hypo_state->pool, &pick_state, 1 );
							copy_associated_picks( &hypo_state->pool, &main_pool, 1 );
							mark_primary_picks( &hypo_state->pool );
						/* Mark the last pick insert time for temporary use, not the real trigger time */
							hypo_state->trigger_time = el_misc_timenow();
							hypo_state->flag = HYPO_IS_WAITING;
							hypo_state->eid = gen_hypo_eid( hypo_state );
						/* */
							insert_hypo_to_pool( &hypo_pool, hypo_state );
						}
					}
				}
			/* Starting the hypo processing thread... */
				bootstrap_hypo_in_pool( &hypo_pool );
			/* All the picks should exist in the main pool */
				insert_pick_to_pool( &main_pool, &pick_state, 0 );
			}
		} while ( 1 ); /* end of message-processing-loop */
	/* no more messages; wait for new ones to arrive */
		sleep_ew(50);
	}
/*-----------------------------end of main loop-------------------------------*/
exit_procedure:
	earlyloc_end();
	destroy_hypo_pool( &hypo_pool );
	destroy_pick_pool( &main_pool );
	ew_unlockfile(lockfile_fd);
	ew_unlink_lockfile(lockfile);

	return 0;
}

/*
 * earlyloc_config() processes command file(s) using kom.c functions;
 *                     exits if any errors are encountered.
 */
static void earlyloc_config( char *configfile )
{
	char  init[23];     /* init flags, one byte for each required command */
	char *com;
	char *str;

	int    ncommand;     /* # of required commands you expect to process   */
	int    nmiss;        /* number of required commands that were missed   */
	int    nfiles;
	int    success;
	int    i;
	double _val = 0.0;
	char   filepath[MAX_PATH_STR];

#define X(a, b) b,
	const char *strategy[] = {
		PICK_FETCH_STRATEGY_TABLE
	};
#undef X

/* Set to zero one init flag for each required command */
	ncommand = 23;
	for ( i = 0; i < ncommand; i++ )
		init[i] = 0;

/* Open the main configuration file */
	nfiles = k_open( configfile );
	if ( nfiles == 0 ) {
		logit("e", "earlyloc: Error opening command file <%s>; exiting!\n", configfile);
		exit(-1);
	}
/* Process all command files */
/* While there are command files open */
	while ( nfiles > 0 ) {
	/* Read next line from active file  */
		while ( k_rd() ) {
		/* Get the first token from line */
			com = k_str();
		/* Ignore blank lines & comments */
			if ( !com )
				continue;
			if ( com[0] == '#' )
				continue;
		/* Open a nested configuration file */
			if ( com[0] == '@' ) {
				success = nfiles + 1;
				nfiles  = k_open(&com[1]);
				if ( nfiles != success ) {
					logit("e", "earlyloc: Error opening command file <%s>; exiting!\n", &com[1]);
					exit(-1);
				}
				continue;
			}

		/* Process anything else as a command */
		/* 0 */
			if( k_its("LogFile") ) {
				LogSwitch = k_int();
				init[0] = 1;
			}
		/* 1 */
			else if( k_its("MyModuleId") ) {
				if ( (str = k_str()) )
					strcpy( MyModName, str );
				init[1] = 1;
			}
		/* 2 */
			else if( k_its("InRingName") ) {
				if ( (str = k_str()) )
					strcpy( InRingName, str );
				init[2] = 1;
			}
		/* 3 */
			else if( k_its("OutRingName") ) {
				if ( (str = k_str()) )
					strcpy( OutRingName, str );
				init[3] = 1;
			}
		/* 4 */
			else if( k_its("HeartBeatInterval") ) {
				HeartBeatInterval = k_long();
				init[4] = 1;
			}
			else if( k_its("ReportPath") ) {
				if ( (str = k_str()) )
					strcpy( ReportPath, str );
				if ( ReportPath[strlen(ReportPath) - 1] != '/' )
					strncat(ReportPath, "/", 1);
				if ( !access(ReportPath, F_OK | W_OK ) ) {
					logit("o", "earlyloc: Output the report file to '%s'.\n", ReportPath);
				}
				else {
					logit("e", "earlyloc: ERROR with the output path: '%s'!\n", ReportPath);
					exit(-1);
				}
			}
			else if( k_its("OutputPostfix") ) {
				if ( (str = k_str()) )
					strcpy( OutputPostfix, str );
				logit("o", "earlyloc: Output report will add '%s' as postfix.\n", OutputPostfix);
			}
			else if( k_its("ReportTermNum") ) {
				ReportTermNum = k_int();
				logit("o", "earlyloc: Report output will terminate at the seq: %d.\n", ReportTermNum ? ReportTermNum : MAX_ALLOW_REPCOUNT);
			}
			else if( k_its("OutputRejectPick") ) {
				OutputRejectPick = k_int();
				logit("o", "earlyloc: Rejected picks %s output to the report file & ring.\n", OutputRejectPick ? "will" : "will not");
			}
			else if( k_its("OutputPickRing") ) {
				OutputPickRing = k_int();
				logit("o", "earlyloc: Pick result %s output to the ring.\n", OutputPickRing ? "will" : "will not");
			}
			else if( k_its("IgnoreWeightP") ) {
				IgnoreWeightP = k_int();
				logit("o", "earlyloc: P-phase picks with weight larger than & equal to %d will be ignored.\n", IgnoreWeightP);
			}
			else if( k_its("IgnoreWeightS") ) {
				IgnoreWeightS = k_int();
				logit("o", "earlyloc: S-phase picks with weight larger than & equal to %d will be ignored.\n", IgnoreWeightS);
			}
			else if( k_its("PickAliveTime") ) {
				_val = k_val();
				logit("o", "earlyloc: Picks in the main pool will alive for %.2lf sec. (default is %.2lf).\n", _val, PickAliveTime);
				PickAliveTime = _val;
			}
			else if( k_its("PickFetchStrategy") ) {
				if ( (str = k_str()) ) {
					for ( i = 0; i < PICK_FETCH_STRATEGY_COUNT; i++ ) {
						if ( !strcmp(str, strategy[i]) )
							break;
					}
					if ( i < PICK_FETCH_STRATEGY_COUNT ) {
						PickFetchStrategy = i;
					}
				}
				logit("o", "earlyloc: Using the '%s' fetching strategy of picks.\n", strategy[PickFetchStrategy]);
			}
			else if( k_its("TriggerPicks") ) {
				i = k_int();
				logit("o", "earlyloc: Triggering picks number is change to %d (default is %d)\n", i, TriggerPicks);
				TriggerPicks = i;
			}
			else if( k_its("HypoAliveTime") ) {
				_val = k_val();
				logit("o", "earlyloc: Hypo process will alive for %.2lf sec. after locating (default is %.2lf).\n", _val, HypoAliveTime);
				HypoAliveTime = _val;
			}
		/* 5 */
			else if( k_its("ClusterTimeDiff") ) {
				ClusterTimeDiff = k_val();
				init[5] = 1;
			}
		/* 6 */
			else if( k_its("ClusterDist") ) {
				ClusterDist = k_val();
				init[6] = 1;
			}
			else if( k_its("ClusterPicks") ) {
				i = k_int();
				logit("o", "earlyloc: Clustering picks number is change to %d (default is %d)\n", i, ClusterPicks);
				ClusterPicks = i;
			}
		/* 7 */
			else if( k_its("PWaveBoundary") ) {
				PWaveModel.boundary = k_val();
				init[7] = 1;
			}
		/* 8 */
			else if( k_its("ShallowPWaveVel") ) {
				PWaveModel.shallow_init = k_val();
				init[8] = 1;
			}
		/* 9 */
			else if( k_its("ShallowPWaveGrad") ) {
				PWaveModel.shallow_grad = k_val();
				init[9] = 1;
			}
		/* 10 */
			else if( k_its("DeepPWaveVel") ) {
				PWaveModel.deep_init = k_val();
				init[10] = 1;
			}
		/* 11 */
			else if( k_its("DeepPWaveGrad") ) {
				PWaveModel.deep_grad = k_val();
				init[11] = 1;
			}
		/* 12 */
			else if( k_its("SWaveBoundary") ) {
				SWaveModel.boundary = k_val();
				init[12] = 1;
			}
		/* 13 */
			else if( k_its("ShallowSWaveVel") ) {
				SWaveModel.shallow_init = k_val();
				init[13] = 1;
			}
		/* 14 */
			else if( k_its("ShallowSWaveGrad") ) {
				SWaveModel.shallow_grad = k_val();
				init[14] = 1;
			}
		/* 15 */
			else if( k_its("DeepSWaveVel") ) {
				SWaveModel.deep_init = k_val();
				init[15] = 1;
			}
		/* 16 */
			else if( k_its("DeepSWaveGrad") ) {
				SWaveModel.deep_grad = k_val();
				init[16] = 1;
			}
			else if ( k_its("3DVelocityModelFile") ) {
				str = k_str();
				if ( str )
					strcpy(filepath, str);
				logit("o", "earlyloc: 3D velocity model file: %s\n", filepath);

				if ( el_loc_3dvelmod_load( filepath ) ) {
					logit("e", "earlyloc: Error reading 3D velocity model file; exiting!\n");
					exit(-1);
				}
				else {
					logit("o", "earlyloc: Reading 3D velocity model file finish!\n");
				}
			}
			else if ( k_its("SQLHost") ) {
				str = k_str();
				if ( str )
					strcpy(DBInfo.host, str);
#if defined( _USE_SQL )
				for ( i = 18; i < ncommand; i++ )
					init[i] = 0;
#endif
			}
		/* 18 */
			else if ( k_its("SQLPort") ) {
				DBInfo.port = k_long();
				init[18] = 1;
			}
		/* 19 */
			else if ( k_its("SQLUser") ) {
				str = k_str();
				if ( str )
					strcpy(DBInfo.user, str);
				init[19] = 1;
			}
		/* 20 */
			else if ( k_its("SQLPassword") ) {
				str = k_str();
				if ( str )
					strcpy(DBInfo.password, str);
				init[20] = 1;
			}
		/* 21 */
			else if ( k_its("SQLDatabase") ) {
				str = k_str();
				if ( str )
					strcpy(DBInfo.database, str);
				init[21] = 1;
			}
		/* 22 */
			else if ( k_its("SQLStationTable") ) {
				if ( nList >= MAXLIST ) {
					logit("e", "earlyloc: Too many <SQLStationTable> commands in <%s>", configfile);
					logit("e", "; max=%d; exiting!\n", (int)MAXLIST);
					exit(-1);
				}
				if ( (str = k_str()) )
					strcpy(SQLStationTable[nList], str);
				nList++;
				init[22] = 1;
			}
			else if ( k_its("UseSNL") ) {
				str = k_get();
				for ( str += strlen(str) + 1; isspace(*str); str++ );
				if ( el_list_sta_line_parse( str, EARLYLOC_LIST_INITIALIZING ) ) {
					logit(
						"e", "earlyloc: ERROR, lack of some stations information for in <%s>. Exiting!\n",
						configfile
					);
					exit(-1);
				}
			}
		/* Enter installation & module to get event messages from */
		/* 17 */
			else if( k_its("GetEventsFrom") ) {
				if ( nLogo + 1 >= MAXLOGO ) {
					logit("e", "earlyloc: Too many <GetEventsFrom> commands in <%s>", configfile);
					logit("e", "; max=%d; exiting!\n", (int) MAXLOGO);
					exit(-1);
				}
				if ( ( str = k_str() ) ) {
					if ( GetInst(str, &Getlogo[nLogo].instid) != 0 ) {
						logit("e", "earlyloc: Invalid installation name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				if ( ( str = k_str() ) ) {
					if ( GetModId(str, &Getlogo[nLogo].mod) != 0 ) {
						logit("e", "earlyloc: Invalid module name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				if ( ( str = k_str() ) ) {
					if ( GetType(str, &Getlogo[nLogo].type) != 0 ) {
						logit("e", "earlyloc: Invalid message type name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				nLogo++;
				init[17] = 1;
			}
		 /* Unknown command */
			else {
				logit("e", "earlyloc: <%s> Unknown command in <%s>.\n", com, configfile);
				continue;
			}

		/* See if there were any errors processing the command */
			if ( k_err() ) {
				logit("e", "earlyloc: Bad <%s> command in <%s>; exiting!\n", com, configfile);
				exit(-1);
			}
		}
		nfiles = k_close();
	}

/* After all files are closed, check init flags for missed commands */
	nmiss = 0;
	for ( i = 0; i < ncommand; i++ )
		if ( !init[i] )
			nmiss++;
	if ( nmiss ) {
		logit("e", "earlyloc: ERROR, no ");
		if ( !init[0] )  logit("e", "<LogFile> "              );
		if ( !init[1] )  logit("e", "<MyModuleId> "           );
		if ( !init[2] )  logit("e", "<InRingName> "           );
		if ( !init[3] )  logit("e", "<OutRingName> "          );
		if ( !init[4] )  logit("e", "<HeartBeatInterval> "    );
		if ( !init[5] )  logit("e", "<ClusterTimeDiff> "      );
		if ( !init[6] )  logit("e", "<ClusterDist> "          );
		if ( !init[7] )  logit("e", "<PWaveBoundary> "        );
		if ( !init[8] )  logit("e", "<ShallowPWaveVel> "      );
		if ( !init[9] )  logit("e", "<ShallowPWaveGrad> "     );
		if ( !init[10] ) logit("e", "<DeepPWaveVel> "         );
		if ( !init[11] ) logit("e", "<DeepPWaveGrad> "        );
		if ( !init[12] ) logit("e", "<SWaveBoundary> "        );
		if ( !init[13] ) logit("e", "<ShallowSWaveVel> "      );
		if ( !init[14] ) logit("e", "<ShallowSWaveGrad> "     );
		if ( !init[15] ) logit("e", "<DeepSWaveVel> "         );
		if ( !init[16] ) logit("e", "<DeepSWaveGrad> "        );
		if ( !init[17] ) logit("e", "any <GetEventsFrom> "    );
		if ( !init[18] ) logit("e", "<SQLPort> "              );
		if ( !init[19] ) logit("e", "<SQLUser> "              );
		if ( !init[20] ) logit("e", "<SQLPassword> "          );
		if ( !init[21] ) logit("e", "<SQLDatabase> "          );
		if ( !init[22] ) logit("e", "any <SQLStationTable> "  );
		logit("e", "command(s) in <%s>; exiting!\n", configfile);
		exit(-1);
	}

	return;
}

/*
 * earlyloc_lookup() - Look up important info from earthworm.h tables
 */
static void earlyloc_lookup( void )
{
/* Look up keys to shared memory regions */
	if ( ( InRingKey = GetKey(InRingName) ) == -1 ) {
		fprintf(stderr, "earlyloc: Invalid ring name <%s>; exiting!\n", InRingName);
		exit(-1);
	}
	if ( ( OutRingKey = GetKey(OutRingName) ) == -1 ) {
		fprintf(stderr, "earlyloc: Invalid ring name <%s>; exiting!\n", OutRingName);
		exit(-1);
	}
/* Look up installations of interest */
	if ( GetLocalInst(&InstId) != 0 ) {
		fprintf(stderr, "earlyloc: Error getting local installation id; exiting!\n");
		exit(-1);
	}
/* Look up modules of interest */
	if ( GetModId(MyModName, &MyModId) != 0 ) {
		fprintf(stderr, "earlyloc: Invalid module name <%s>; exiting!\n", MyModName);
		exit(-1);
	}
/* Look up message types of interest */
	if ( GetType("TYPE_HEARTBEAT", &TypeHeartBeat) != 0 ) {
		fprintf(stderr, "earlyloc: Invalid message type <TYPE_HEARTBEAT>; exiting!\n");
		exit(-1);
	}
	if ( GetType("TYPE_ERROR", &TypeError) != 0) {
		fprintf(stderr, "earlyloc: Invalid message type <TYPE_ERROR>; exiting!\n");
		exit(-1);
	}
	if ( GetType("TYPE_EARLY_PICK", &TypePickInput) != 0 ) {
		fprintf(stderr, "earlyloc: Invalid message type <TYPE_EARLY_PICK>; exiting!\n");
		exit(-1);
	}
	if ( GetType("TYPE_EARLY_EVENT", &TypeHypoOutput) != 0 ) {
		fprintf(stderr, "earlyloc: Invalid message type <TYPE_EARLY_EVENT>; exiting!\n");
		exit(-1);
	}

	return;
}

/*
 * earlyloc_status() - Builds a heartbeat or error message & puts it into
 *                     shared memory.  Writes errors to log file & screen.
 */
static void earlyloc_status( uint8_t type, short ierr, char *note )
{
	MSG_LOGO    logo;
	char        msg[512];
	uint64_t    size;
	time_t      t;

/* Build the message */
	logo.instid = InstId;
	logo.mod    = MyModId;
	logo.type   = type;

	time(&t);

	if ( type == TypeHeartBeat ) {
		sprintf(msg, "%ld %ld\n", (long)t, (long)myPid);
	}
	else if ( type == TypeError ) {
		sprintf(msg, "%ld %hd %s\n", (long)t, ierr, note);
		logit("et", "earlyloc: %s\n", note);
	}

	size = strlen(msg);   /* don't include the null byte in the message */

/* Write the message to shared memory */
	if ( tport_putmsg(&InRegion, &logo, size, msg) != PUT_OK ) {
		if ( type == TypeHeartBeat ) {
			logit("et","earlyloc: Error sending heartbeat.\n");
		}
		else if ( type == TypeError ) {
			logit("et","earlyloc: Error sending error:%d.\n", ierr);
		}
	}

	return;
}

/*
 * earlyloc_end()  free all the local memory & close socket
 */
static void earlyloc_end( void )
{
	tport_detach(&InRegion);
	tport_detach(&OutRegion);
	el_loc_3dvelmod_free();

	return;
}

/*
 *
 */
static int thread_proc_trigger( void *arg )
{
	HYPO_STATE *const result = (HYPO_STATE *)arg;
	struct timespec waittime = { .tv_sec = 0, .tv_nsec = 100000000 }; /* 100 ms */
	PICK_STATE     *pick = NULL;
	PICK_STATE     *pick_in_pool = NULL;
	MSG_LOGO        logo = { 0 };
	char            _report_path[MAX_PATH_STR] = { 0 };
	int             max_valids = 0;
	int             pool_status = POOL_HAS_NEW_PICK;
	double          min_error = CONVERGE_CRITERIA;
	double          search_begin;
	double          time_last_hypo = el_misc_timenow();

/* Build the message */
	logo.instid = InstId;
	logo.mod    = MyModId;
	logo.type   = TypeHypoOutput;
/* Mark the real trigger time */
	result->trigger_time = time_last_hypo;
	result->rep_count = 0;
/* */
	if ( strlen(ReportPath) )
		mk_outdir_by_evt( _report_path, ReportPath, result->trigger_time, result->eid, OutputPostfix );
/* */
	do {
	/* */
		if ( pool_status == POOL_HAS_NEW_PICK ) {
		/* Initial guess */
			if ( result->ig_origin_time < 0.0 ) {
				el_loc_location_guess( result, &PWaveModel, &SWaveModel );
				el_loc_location_refine( result, &PWaveModel, &SWaveModel );
				save_to_init_guess( result );
			}
		/* Main hypo iteration process... */
			use_init_guess( result );
			search_begin = result->origin_time;
			do {
			/* The main locate process */
				if ( el_loc_primary_locate( result, &PWaveModel, &SWaveModel ) ) {
					logit("et", "earlyloc: Hypo(#%d) locating ERROR, skip this time!\n", result->eid);
					break;
				}
			/* Check if it converge or not */
				if ( HYPO_IS_CONVERGED( result ) ) {
				/* Adjust the origin time to reduce the overall residual */
					el_loc_origintime_adjust( result, &PWaveModel, &SWaveModel );
					el_loc_all_states_update( result, &PWaveModel, &SWaveModel );
				/* Finally, output the report... */
					if ( !ReportTermNum || result->rep_count < ReportTermNum ) {
						el_report_ring_output( result, &OutRegion, &logo, OutputPostfix, OutputRejectPick );
						if ( strlen(_report_path) )
							el_report_file_output( result, _report_path, OutputPostfix, OutputRejectPick );
					}
				/* Keep the best solution for next initial guess */
					if (
						result->pool.valids > max_valids ||
						(result->pool.valids == max_valids && result->avg_error < min_error)
					) {
						save_to_best_state( result );
						max_valids = result->pool.valids;
						min_error = result->avg_error;
					/* */
						if ( result->pool.valids > 12 )
							save_to_init_guess( result );
					}
				/* Increase the report # */
					result->rep_count++;
					break;
				}
			/* */
				use_init_guess( result );
				if ( !result->rep_count && (search_begin - result->origin_time) < 20.0 ) {
					result->origin_time -= 1.0;
					unmark_locmask_picks( &result->pool, 0.0 );
					el_loc_location_refine( result, &PWaveModel, &SWaveModel );
					save_to_init_guess( result );
					continue;
				}
				el_loc_all_states_update( result, &PWaveModel, &SWaveModel );
			/* Mask the pick with highest residual until it can't do it any more */
				if ( !mark_locmask_pick_max_residual( &result->pool ) )
					break;
			} while ( result->pool.valids >= MIN_LOCATE_PICKS );
		/* */
			if ( result->pool.valids < MIN_LOCATE_PICKS )
				printf("earlyloc: ##### NO REPORT ##### Hypo(#%d) can't reach the converging criteria!\n", result->eid);
		/* Debug null output */
		#ifdef _DEBUG
			if ( strlen(_report_path) )
				el_report_file_output( result, _report_path, OutputPostfix, OutputRejectPick );
		#endif
		/* Unmask all picks */
			unmark_locmask_picks( &result->pool, 0.0 );
		/* */
			if ( !result->rep_count || (result->pool.valids < 12 && result->avg_weight < 0.5) ) {
			/* */
				reset_init_guess( result );
			}
		/* */
			if ( result->rep_count > MIN_STABLE_REPCOUNT ) {
			/* Don't reject any picks before the stable result */
				mark_reject_picks( &result->pool, REJECT_CRITERIA );
			}
			else if ( result->rep_count > MAX_ALLOW_REPCOUNT ) {
			/* Finish the hypo process when over the report count limit */
				break;
			}
		/* */
			pool_status = POOL_ALREADY_USED;
			time_last_hypo = el_misc_timenow();
		}
		else if ( (el_misc_timenow() - time_last_hypo) >= HypoAliveTime ) {
			break;
		}
		else {
			thrd_yield();
			thrd_sleep(&waittime, NULL);
		}
	/* Looking for the new pick in the pool queue, include those picks added in the previous step */
		while ( result->flag != HYPO_IS_FINISHED ) {
		/* Using trylock instead lock to avoid the dead lock */
			if ( mtx_trylock(&result->queue_mutex) != thrd_success )
				break;
		/* Real process after we gor the mutex */
			pick = pop_pick_from_pool( &result->pick_queue );
			mtx_unlock(&result->queue_mutex);
			if ( pick ) {
			/* */
				pick_in_pool = insert_pick_to_pool( &result->pool, pick, 1 );
				free(pick);
			/* */
				if ( EL_PICK_VALID_LOCATE( pick_in_pool ) ) {
					pool_status = POOL_HAS_NEW_PICK;
				/* */
					if ( PickFetchStrategy == PICK_FETCH_STRATEGY_STEP ) {
						pick_in_pool = NULL;
						break;
					}
				}
			}
			else if ( pick_in_pool ) {
				pick_in_pool = NULL;
				thrd_yield();
			}
			else {
				break;
			}
		}
	} while ( result->flag != HYPO_IS_FINISHED );
/* End process */
	result->flag = HYPO_IS_FINISHED;
	logit("ot", "earlyloc: Finished hypo(#%d) at the end of hypo life.\n", result->eid);

	return 0;
}

/**
 * @brief
 *
 * @param target
 * @param usesnl
 * @return EARLY_PICK_MSG*
 */
static EARLY_PICK_MSG *fill_coor2epick( EARLY_PICK_MSG *target, const USE_SNL *usesnl )
{
/* */
	target->latitude = usesnl->latitude;
	target->longitude = usesnl->longitude;
	target->elevation = usesnl->elevation * -0.001;

	return target;
}


/**
 * @brief
 *
 * @param dest
 * @param src
 * @return PICK_STATE*
 */
static PICK_STATE *parse_epick2pickstate( PICK_STATE *dest, const EARLY_PICK_MSG *src )
{
/* */
	memcpy(&dest->observe, src, sizeof(EARLY_PICK_MSG));
/* */
	dest->flag = 0;
	dest->recvtime = el_misc_timenow();
	dest->trv_time = 0.0;
	dest->distance = 0.0;
	dest->residual = 0.0;
	dest->r_weight = 0.0;

	return dest;
}

/**
 * @brief Create a hypo state object
 *
 * @return HYPO_STATE*
 */
static HYPO_STATE *create_hypo_state( void )
{
	HYPO_STATE *result = (HYPO_STATE *)calloc(1, sizeof(HYPO_STATE));

	result->eid          = 0;
	result->tid          = 0;
	result->flag         = HYPO_IS_UNUSED;
	result->rep_count    = 0;
	result->latitude     = 0.0;
	result->longitude    = 0.0;
	result->depth        = 0.0;
	result->origin_time  = -1.0;
	result->trigger_time = -1.0;
	result->gap          = 360.0;
	result->avg_error    = REJECT_CRITERIA;
	result->avg_weight   = REJECT_CRITERIA;
	result->q            = HYPO_RESULT_QUALITY_D;
/* */
	result->ig_latitude    = 0.0;
	result->ig_longitude   = 0.0;
	result->ig_depth       = 0.0;
	result->ig_origin_time = -1.0;
/* */
	EL_PICK_POOL_INIT( result->pool );
	EL_PICK_POOL_INIT( result->pick_queue );
	if ( mtx_init(&result->queue_mutex, mtx_plain) != thrd_success )
		return NULL;
/* */
	if ( !(result->best = (HYPO_STATE *)calloc(1, sizeof(HYPO_STATE))) )
		return NULL;

	return result;
}

/**
 * @brief
 *
 * @param target
 * @return HYPO_STATE*
 */
static HYPO_STATE *save_to_best_state( HYPO_STATE *target )
{
/* */
	if ( target->best ) {
		target->best->eid          = target->eid;
		target->best->tid          = target->tid;
		target->best->flag         = target->flag;
		target->best->rep_count    = target->rep_count;
		target->best->latitude     = target->latitude;
		target->best->longitude    = target->longitude;
		target->best->depth        = target->depth;
		target->best->origin_time  = target->origin_time;
		target->best->trigger_time = target->trigger_time;
		target->best->gap          = target->gap;
		target->best->avg_error    = target->avg_error;
		target->best->avg_weight   = target->avg_weight;
		target->best->q            = target->q;
	/* */
		target->best->ig_latitude    = target->ig_latitude;
		target->best->ig_longitude   = target->ig_longitude;
		target->best->ig_depth       = target->ig_depth;
		target->best->ig_origin_time = target->ig_origin_time;
	}

	return target;
}

/**
 * @brief
 *
 * @param target
 * @return HYPO_STATE*
 */
static HYPO_STATE *use_init_guess( HYPO_STATE *target )
{
	target->latitude    = target->ig_latitude;
	target->longitude   = target->ig_longitude;
	target->depth       = target->ig_depth;
	target->origin_time = target->ig_origin_time;

	return target;
}

/**
 * @brief
 *
 * @param target
 * @return HYPO_STATE*
 */
static HYPO_STATE *save_to_init_guess( HYPO_STATE *target )
{
	target->ig_latitude    = target->latitude;
	target->ig_longitude   = target->longitude;
	target->ig_depth       = target->depth;
	target->ig_origin_time = target->origin_time;

	return target;
}

/**
 * @brief
 *
 * @param target
 * @return HYPO_STATE*
 */
static HYPO_STATE *reset_init_guess( HYPO_STATE *target )
{
	target->ig_latitude    = 0.0;
	target->ig_longitude   = 0.0;
	target->ig_depth       = 0.0;
	target->ig_origin_time = -1.0;

	return target;
}

/**
 * @brief
 *
 * @param pool
 * @param input_pick
 */
static void associate_pick_hypos( HYPOS_POOL *pool, PICK_STATE *input_pick )
{
	DL_NODE    *node;
	HYPO_STATE *hyp;

	DL_LIST_FOR_EACH_DATA( pool->entry, node, hyp ) {
		if ( hyp->flag == HYPO_IS_PROCESSING || hyp->flag == HYPO_IS_WAITING ) {
			if (
			/* The pick with the same SCNL & phase should not add in the pool again */
				!check_scnlp_exist_pool( &hyp->pool, input_pick ) &&
			/* Not only insert the clustered picks here, also insert co-sited picks */
				(PICK_IS_ASSOCIATED_WITH_POOL( &hyp->pool, input_pick ) ||
			/* Try to associate those pick with low time residual in the main pool */
				(hyp->rep_count && fabs(el_loc_residual_estimate( hyp->best, input_pick, &PWaveModel, &SWaveModel )) <= ROUGH_ASSOC_CRITERIA))
			) {
				EL_MARK_PICK_INUSE( input_pick );
				if ( hyp->flag == HYPO_IS_PROCESSING ) {
				/* When hypo processing, insert the pick to the queue of pool */
					mtx_lock(&hyp->queue_mutex);
					insert_pick_to_pool( &hyp->pick_queue, input_pick, 0 );
					mtx_unlock(&hyp->queue_mutex);
				}
				else {
				/* When hypo waiting, insert the pick to the pool */
					insert_pick_to_pool( &hyp->pool, input_pick, 1 );
				/* Again, mark the last pick insert time, not the real trigger time */
					hyp->trigger_time = el_misc_timenow();
				}
			}
		}
	}

	return;
}

/**
 * @brief
 *
 * @param pool
 * @param input
 * @return HYPO_STATE*
 */
static HYPO_STATE *insert_hypo_to_pool( HYPOS_POOL *pool, HYPO_STATE *input )
{
/* */
	pool->last = dl_node_append( &pool->entry, input );
	pool->totals++;

	return input;
}

/**
 * @brief
 *
 * @param pool
 */
static void bootstrap_hypo_in_pool( HYPOS_POOL *pool )
{
	DL_NODE    *node;
	HYPO_STATE *hyp;

	DL_LIST_FOR_EACH_DATA( pool->entry, node, hyp ) {
		if ( hyp->flag == HYPO_IS_WAITING && hyp->pool.valids >= TriggerPicks ) {
		/* Creating the thread for hypo... */
			logit("ot", "earlyloc: There is a trigger, creating a hypo(#%d) thread to process...\n", hyp->eid);
			if ( thrd_create(&hyp->tid, thread_proc_trigger, hyp) != thrd_success )
				logit("e", "earlyloc: Error starting trigger processing thread; skip it!\n");
			else if ( thrd_detach(hyp->tid) != thrd_success )
				logit("e", "earlyloc: Error detaching trigger processing thread for hypo(#%d); notice it!\n", hyp->eid);
			else
				hyp->flag = HYPO_IS_PROCESSING;
		}
	}
}

/**
 * @brief
 *
 * @param pool
 */
static void close_obsolete_hypos( HYPOS_POOL *pool )
{
	DL_NODE    *node, *safe;
	HYPO_STATE *hyp;
	double      time_now = el_misc_timenow();

	DL_LIST_FOR_EACH_DATA_SAFE( pool->entry, node, hyp, safe ) {
		if (
			hyp->flag == HYPO_IS_FINISHED ||
			(hyp->flag == HYPO_IS_WAITING && (time_now - hyp->trigger_time) >= ClusterTimeDiff)
		) {
			logit("ot", "earlyloc: Closed the finished(obsoleted) hypo(#%d)!\n", hyp->eid);
			destroy_pick_pool( &hyp->pick_queue );
			destroy_pick_pool( &hyp->pool );
			mtx_destroy(&hyp->queue_mutex);
			dl_node_delete( node, free );
		/* */
			if ( node == pool->entry )
				pool->entry = safe;
			if ( node == pool->last )
				pool->last = NULL;
		/* */
			pool->totals--;
		}
	}

	return;
}

/**
 * @brief
 *
 * @param pool
 * @return HYPOS_POOL*
 */
static HYPOS_POOL *destroy_hypo_pool( HYPOS_POOL *pool )
{
	DL_NODE    *node, *safe;
	HYPO_STATE *hyp;

	DL_LIST_FOR_EACH_DATA_SAFE( pool->entry, node, hyp, safe ) {
		destroy_pick_pool( &hyp->pick_queue );
		destroy_pick_pool( &hyp->pool );
		mtx_destroy(&hyp->queue_mutex);
		dl_node_delete( node, free );
	}
	EL_HYPO_POOL_INIT( *pool );

	return pool;
}

/**
 * @brief
 *
 * @param pool
 * @return PICKS_POOL*
 */
static PICKS_POOL *destroy_pick_pool( PICKS_POOL *pool )
{
	dl_list_destroy( &pool->entry, free );
	EL_PICK_POOL_INIT( *pool );

	return pool;
}

/**
 * @brief
 *
 * @param pool
 * @param input
 * @return int
 */
static int check_pick_exist_pool( const PICKS_POOL *pool, const PICK_STATE *input )
{
	DL_NODE    *node;
	PICK_STATE *pick;

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( !compare_pick( pick, input ) ) {
			return 1;
		}
	}

	return 0;
}

/**
 * @brief
 *
 * @param pool
 * @param input
 * @return int
 */
static int check_scnlp_exist_pool( const PICKS_POOL *pool, const PICK_STATE *input )
{
	DL_NODE    *node;
	PICK_STATE *pick;

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( !compare_scnl( pick, input ) && !strcmp( pick->observe.phase_name, input->observe.phase_name ) ) {
			return 1;
		}
	}

	return 0;
}

/**
 * @brief
 *
 * @param pool
 * @param input
 * @return int
 */
static int check_pick_cluster_pool( const PICKS_POOL *pool, const PICK_STATE *input )
{
	DL_NODE    *node;
	PICK_STATE *pick;
	int         result = 0;

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( PICKS_ARE_CLUSTER( pick, input ) ) {
			result++;
		/* */
			if ( pool->totals < ClusterPicks && result == ClusterPicks - 1 )
				result++;
		/* ... */
			if ( result >= ClusterPicks )
				break;
		}
	}

	return result >= ClusterPicks ? 1 : 0;
}

/**
 * @brief
 *
 * @param pool
 * @param input
 * @return int
 */
static int check_pick_cosite_pool( const PICKS_POOL *pool, const PICK_STATE *input )
{
	DL_NODE    *node;
	PICK_STATE *pick;

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		if ( PICKS_ARE_COSITE( pick, input ) ) {
			return 1;
		}
	}

	return 0;
}

/**
 * @brief
 *
 * @param pool
 * @param input
 * @param sort_pick_time
 * @return PICK_STATE*
 */
static PICK_STATE *insert_pick_to_pool( PICKS_POOL *pool, const PICK_STATE *input, const int sort_pick_time )
{
	DL_NODE    *node;
	PICK_STATE *pick  = NULL;
	PICK_STATE *_pick = NULL;
	uint8_t     cosite = 0;
	uint8_t     flag   = 0;

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( !compare_pick( pick, input ) ) {
			break;
		}
		else if ( !cosite ) {
			if ( PICKS_ARE_COSITE( pick, input ) && !(pick->flag & PICK_FLAG_COSITE) ) {
				flag = pick->flag;
			/* Switch the cosite flag if the new incoming pick is better than the pick existed in pool */
				if (
					(!PICK_IS_ON_SURFACE( input ) && PICK_IS_ON_SURFACE( pick )) ||
					((PICK_IS_ON_SURFACE( input ) == PICK_IS_ON_SURFACE( pick )) &&
					(input->observe.weight < pick->observe.weight || (input->observe.weight == pick->observe.weight)))
				) {
				/* We should unmark the locmask then mark the cosite */
					EL_MARK_PICK_COSITE( pick );
				}
				else {
				/* The same as above process on the temp flag */
					EL_MARK_UNMARK_PICK_FLAG( flag, PICK_FLAG_COSITE, PICK_FLAG_LOCMASK );
				}
				cosite = 1;
			}
		}
		pick = NULL;
	}
/* This input should be the new picking */
	if ( !pick ) {
	/* */
		pick = malloc(sizeof(PICK_STATE));
		memcpy(pick, input, sizeof(PICK_STATE));
	/* */
		if ( cosite )
			pick->flag = flag;
		else if ( pick->flag & PICK_FLAG_COSITE )
			pick->flag ^= PICK_FLAG_COSITE;
	/* */
		_pick = (PICK_STATE *)DL_NODE_GET_DATA( pool->last );
		if ( !sort_pick_time || !_pick || pick->observe.picktime >= _pick->observe.picktime ) {
			pool->last = dl_node_append( &pool->entry, pick );
		}
		else {
			DL_LIST_FOR_EACH_DATA( pool->entry, node, _pick ) {
				if ( pick->observe.picktime < _pick->observe.picktime ) {
					if ( (node = DL_NODE_GET_PREV( node )) )
						dl_node_insert( node, pick );
					else
						dl_node_push( &pool->entry, pick );
				/* */
					break;
				}
			}
		}
	/* */
		pool->totals++;
	/* */
		if ( cosite )
			pool->cosites++;
		else if ( EL_PICK_VALID_LOCATE( pick ) )
			pool->valids++;
	/* */
		if ( pick->flag & PICK_FLAG_REJECT )
			pool->rejects++;

	}

	return pick;
}

/**
 * @brief
 *
 * @param pool
 * @return PICK_STATE*
 */
static PICK_STATE *pop_pick_from_pool( PICKS_POOL *pool )
{
	PICK_STATE *result = NULL;

/* */
	if ( (result = (PICK_STATE *)dl_node_data_extract( dl_node_pop( &pool->entry ) )) ) {
		pool->totals--;
		if ( EL_PICK_VALID_LOCATE( result ) ) {
			pool->valids--;
		}
		else {
			if ( result->flag & PICK_FLAG_COSITE )
				pool->cosites--;
			if ( result->flag & PICK_FLAG_REJECT )
				pool->rejects--;
		}
	}
/* Properly close this pool when there is not any pick in pool */
	if ( !pool->entry )
		pool->last = pool->entry;

	return result;
}

/**
 * @brief
 *
 * @param pool
 * @return PICK_STATE*
 */
static PICK_STATE *mark_locmask_pick_max_residual( PICKS_POOL *pool )
{
	DL_NODE    *node;
	PICK_STATE *pick;
	PICK_STATE *target = NULL;
	double      max_residual = -1.0;

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( EL_PICK_VALID_LOCATE( pick ) && !(pick->flag & PICK_FLAG_PRIMARY) ) {
		/* */
			if ( fabs(pick->residual) > max_residual ) {
				max_residual = fabs(pick->residual);
				target = pick;
			}
		}
	}
/* */
	if ( target ) {
	/* Mask out the pick with maximum residual */
		EL_MARK_PICK_LOCMASK( target );
	/* Then, decrease the valid count by 1 */
		pool->valids--;
	}

	return target;
}

/**
 * @brief
 *
 * @param pool
 */
static void remove_obsolete_picks( PICKS_POOL *pool )
{
	DL_NODE    *node, *safe;
	PICK_STATE *pick = NULL;
	double      time_ref;

/* */
	time_ref = el_misc_timenow() - PickAliveTime;
/* */
	DL_LIST_FOR_EACH_DATA_SAFE( pool->entry, node, pick, safe ) {
		if ( pick->recvtime <= time_ref ) {
			if ( EL_PICK_VALID_LOCATE( pick ) ) {
				pool->valids--;
			}
			else {
				if ( pick->flag & PICK_FLAG_COSITE )
					pool->cosites--;
				if ( pick->flag & PICK_FLAG_REJECT )
					pool->rejects--;
			}
			pool->entry = dl_node_delete( node, free );
			pool->totals--;
		}
		else if ( pick->recvtime > time_ref ) {
			break;
		}
	}
/* Properly close this pool when there is not any pick in pool */
	if ( !pool->entry )
		pool->last = pool->entry;

	return;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param skip_inuse
 */
static void copy_associated_picks( PICKS_POOL *dest, const PICKS_POOL *src, const int skip_inuse )
{
	DL_NODE    *node;
	PICK_STATE *pick;

/* */
	DL_LIST_FOR_EACH_DATA( src->entry, node, pick ) {
	/* */
		if ( !skip_inuse || !(pick->flag & PICK_FLAG_INUSE) ) {
			if ( PICK_IS_ASSOCIATED_WITH_POOL( dest, pick ) ) {
				EL_MARK_PICK_INUSE( pick );
				insert_pick_to_pool( dest, pick, 1 );
			}
		}
	}

	return;
}

/**
 * @brief
 *
 * @param pool
 * @param residual
 */
static void mark_reject_picks( PICKS_POOL *pool, double residual )
{
	DL_NODE    *node;
	PICK_STATE *pick;

/* */
	residual = fabs(residual);
/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( !(pick->flag & PICK_FLAG_REJECT) && !(pick->flag & PICK_FLAG_PRIMARY) ) {
			if ( (fabs(pick->residual) > residual) ) {
				if ( EL_PICK_VALID_LOCATE( pick ) )
					pool->valids--;
			/* */
				EL_MARK_PICK_REJECT( pick );
			/* */
				pool->rejects++;
			}
		}
	}

	return;
}

/**
 * @brief
 *
 * @param pool
 * @param residual
 */
static void unmark_locmask_picks( PICKS_POOL *pool, double residual )
{
	DL_NODE    *node;
	PICK_STATE *pick;
	int         num_unmask = 0;

/* */
	residual = fabs(residual);
/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
	/* */
		if ( EL_PICK_LOCMASK_UNMARKABLE( pick ) ) {
			if ( residual < RESIDUAL_EPSILON || fabs(pick->residual) <= residual ) {
				EL_UNMARK_PICK_FLAG( pick->flag, PICK_FLAG_LOCMASK );
				num_unmask++;
			}
		}
	}
/* Then,  */
	pool->valids += num_unmask;

	return;
}

/**
 * @brief
 *
 * @param pool
 */
static void mark_primary_picks( PICKS_POOL *pool )
{
	DL_NODE    *node;
	PICK_STATE *pick;

/* */
	DL_LIST_FOR_EACH_DATA( pool->entry, node, pick ) {
		EL_MARK_PICK_PRIMARY( pick );
	}

	return;
}

/**
 * @brief
 *
 * @param hyp
 * @return int
 */
static int gen_hypo_eid( const HYPO_STATE *hyp )
{
	int         result = 0;
	double      time_now = el_misc_timenow();
	DL_NODE    *node;
	PICK_STATE *pick;
	char       *c = NULL;

	DL_LIST_FOR_EACH_DATA( hyp->pool.entry, node, pick ) {
		if ( pick->flag & PICK_FLAG_PRIMARY ) {
			for ( c = pick->observe.station; *c; c++ )
				result += (int)*c;
		}
	}
/* TODO: It still need some random factor? */
	result = (int)time_now % result;

	return result;
}

/**
 * @brief
 *
 * @param pick_a
 * @param pick_b
 * @return int
 */
static int cluster_pick( const PICK_STATE *pick_a, const PICK_STATE *pick_b )
{
	const double dtime = fabs(pick_a->observe.picktime - pick_b->observe.picktime);
	const double dist  = el_misc_geog2distf( pick_a->observe.longitude, pick_a->observe.latitude, pick_b->observe.longitude, pick_b->observe.latitude );

	if ( dtime < ClusterTimeDiff && dist > PICK_COSITE_DIST && dist < ClusterDist ) {
		if ( !strcmp(pick_a->observe.phase_name, "P") && !strcmp(pick_b->observe.phase_name, "P") ) {
			if ( dtime < (dist / PWaveModel.shallow_init) )
				return 1;
		}
		else if ( !strcmp(pick_a->observe.phase_name, "S") && !strcmp(pick_b->observe.phase_name, "S") ) {
			if ( dtime < (dist / SWaveModel.shallow_init) )
				return 1;
		}
	}

	return 0;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return int
 */
static int compare_pick( const void *a, const void *b )
{
	int result;
	PICK_STATE *tmpa = (PICK_STATE *)a;
	PICK_STATE *tmpb = (PICK_STATE *)b;

/* */
	if ( (result = compare_scnl( tmpa, tmpb )) )
		return result;
	if ( (result = strcmp( tmpa->observe.phase_name, tmpb->observe.phase_name )) )
		return result;
/* Maybe this check is not necessary */
	if ( fabs(tmpa->observe.picktime - tmpb->observe.picktime) > MIN_DIFF_PICK_TIME )
		result = tmpa->observe.picktime > tmpb->observe.picktime ? 1 : -1;

	return result;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return int
 */
static int compare_scnl( const void *a, const void *b )
{
	int result;
	PICK_STATE *tmpa = (PICK_STATE *)a;
	PICK_STATE *tmpb = (PICK_STATE *)b;

/* */
	if ( (result = strcmp( tmpa->observe.station, tmpb->observe.station )) )
		return result;
	if ( (result = strcmp( tmpa->observe.channel, tmpb->observe.channel )) )
		return result;
	if ( (result = strcmp( tmpa->observe.network, tmpb->observe.network )) )
		return result;

	return strcmp( tmpa->observe.location, tmpb->observe.location );
}

/**
 * @brief
 *
 * @param opath
 * @param parent_path
 * @param trigger_time
 * @param eid
 * @param postfix
 * @return int
 */
static int mk_outdir_by_evt( char *opath, const char *parent_path, const double trigger_time, const int eid, const char *postfix )
{
	char timestamp[EL_SIMPLE_TIMESTAMP_BUF_SIZE];

/* If the previous thread is still alived, kill it! */
	if ( postfix && strlen(postfix) ) {
		sprintf(
			opath, "%s%s_%d_%s/", parent_path, el_misc_simple_timestamp_gen( timestamp, EL_SIMPLE_TIMESTAMP_BUF_SIZE, trigger_time ), eid, postfix
		);
	}
	else {
		sprintf(
			opath, "%s%s_%d/", parent_path, el_misc_simple_timestamp_gen( timestamp, EL_SIMPLE_TIMESTAMP_BUF_SIZE, trigger_time ), eid
		);
	}
/* */
	if ( mkdir(opath, S_IRWXU | S_IRGRP | S_IROTH) ) {
		logit("e", "earlyloc: Cannot make the new directory for event (%s), skip it!\n", timestamp);
		return -1;
	}

	return 0;
}
