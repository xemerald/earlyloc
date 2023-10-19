/**
 * @file raytracing.c
 * @author 2009/07/26 by Yih-Min Wu created the 3D velocity earthquake location program
 * @anchor 2017/08/11 by Benjamin Ming Yang modified
 * @brief
 * @version 0.1
 * @date 2023-09-25
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <raytracing.h>

/* Used constant */
#define RT_HALF_PI   1.570796326794896619231321691640f
#define RT_PI        3.141592653589793238462643383279f
#define RT_PI2       6.283185307179586476925286766559f
#define RT_DEG2RAD   0.017453292519943295769236907685f
#define RT_RAD2DEG   57.29577951308232087679815481411f
#define RT_B2A_SQ    0.993305521
#define RT_EPS       1.0e-8  /* Maybe 1.0e-10 */
#define RT_RNULL     0.0e+10

/*
 * Parameters for calculation:
 * XFAC    = enhancement factor (see um and thurber, 1987)
 * NLOOP   = number of bending iterations
 * N1, N2  = min & max of ray segments
 * MINS    = min. length of segment (km)
*/
#define XFAC         1.9
#define N1           2
#define N2           RT_MAX_NODE
#define NLOOP        12800
#define FLIMIT       1.e-6
#define MINS         2.0

/**
 * @brief
 *
 */
typedef struct {
	double x;
	double y;
	double z;
} NODE_INFO;

/**
 * @brief
 *
 */
typedef struct {
	double vel[8];
} VEL_GRID;


#define POLAR2CARTESIAN( _R_MACRO, _THETA_MACRO, _PHI_MACRO, _X_MACRO, _Y_MACRO, _Z_MACRO ) \
		__extension__({ \
			(_X_MACRO) = (_R_MACRO) * sin((_THETA_MACRO)); \
			(_Y_MACRO) = (_X_MACRO) * sin((_PHI_MACRO)); \
			(_X_MACRO) *= cos((_PHI_MACRO)); \
			(_Z_MACRO) = (_R_MACRO) * cos((_THETA_MACRO)); \
		})

/* For 3D-velocity model maping */
static double *Lat_c, *Lon_c, *Dep_c;
static double *Vel_p, *Vel_s;
static VEL_GRID *Velgrid_p, *Velgrid_s;
static VEL_GRID * volatile Vel_grid;

static int *Ilon_c, *Ilat_c, *Idep_c;

/* Constant */
static double Ro, Rs;
static int Ilat1_c, Ilon1_c, Idep1_c;
static int Ilonmax, Ilatmax, Idepmax;
static int Ibld3, Ibld4;
static int Nlat_c, Nlon_c, Ndep_c;
static int Nxyz_c, Nxy_c, Nx_c;

/* */
static void raytracing_pb( double, double, double, double, double, double, RAY_INFO *, int *, double * );
static void step_ray_node( RAY_INFO *, const RAY_INFO *, const double, const double );
static double get_ray_traveltime( const RAY_INFO *, const int );
static double get_vel_ray( const RAY_INFO *, const double );
static double get_vel_geog( const double, const double, const double );
static int intmap_3d( int *, int *, int * );
static int input_vel_model( const char * );
static int bldmap( const double, const double );
static int vel_point2grid( const double *, VEL_GRID *, const int, const int, const int );

static double geog2geoc( const double );
static double geoc2geog( const double );
static double get_earth_radius( double );

/**
 * @brief
 *
 * @param ray_out
 * @param np
 * @param travel_time
 * @param evla
 * @param evlo
 * @param evdp
 * @param stla
 * @param stlo
 * @param stdp
 * @param vel_type
 * @return int
 */
int rt_main( RAY_INFO *ray_out, int *np, double *travel_time, double evla, double evlo, double evdp, double stla, double stlo, double stdp, const int vel_type )
{
/* Check coordinates */
	if ( evla < Lat_c[0] || evla > Lat_c[Nlat_c - 1] ) {
		fprintf(stderr, "rt_main: Latitude(%lf) of source is out of range!\n", evla);
		return -1;
	}
	if ( stla < Lat_c[0] || stla > Lat_c[Nlat_c - 1] ) {
		fprintf(stderr, "rt_main: Latitude(%lf) of station is out of range!\n", stla);
		return -1;
	}
	if ( evlo < Lon_c[0] || evlo > Lon_c[Nlon_c - 1] ) {
		fprintf(stderr, "rt_main: Longitude(%lf) of source is out of range!\n", evlo);
		return -1;
	}
	if ( stlo < Lon_c[0] || stlo > Lon_c[Nlon_c - 1] ) {
		fprintf(stderr, "rt_main: Longitude(%lf) of station is out of range!\n", stlo);
		return -1;
	}
	if ( evdp < Dep_c[0] || evdp > Dep_c[Ndep_c - 1] ) {
		fprintf(stderr, "rt_main: Depth(%lf) of source is out of range!\n", evdp);
		return -1;
	}
	if ( stdp < Dep_c[0] || stdp > Dep_c[Ndep_c - 1] ) {
		fprintf(stderr, "rt_main: Depth(%lf) of station is out of range!\n", stdp);
		return -1;
	}
/* */
	stla = geog2geoc( stla );
	evla = geog2geoc( evla );
/* Define the velocity type, P or S */
	Vel_grid = vel_type == RT_P_WAVE_VELOCITY ? Velgrid_p : Velgrid_s;
/* */
	raytracing_pb( evla, evlo, evdp, stla, stlo, stdp, ray_out, np, travel_time );
/* Show full information */
#ifdef _DEBUG
	printf("# nodes = %d, travel time = %lf\n", np, tt);
	printf("takeoff = %lf, azimuth = %lf\n", *tko, *azi);
#endif

	return 0;
}

/**
 * @brief
 *
 * @param model_path
 * @return int
 */
int rt_velmod_load( const char *model_path )
{
/* */
	if ( input_vel_model( model_path ) )
		return -1;
/* */
	Velgrid_p = calloc(Nxyz_c, sizeof(VEL_GRID));
	Velgrid_s = calloc(Nxyz_c, sizeof(VEL_GRID));
	vel_point2grid( Vel_p, Velgrid_p, Nxyz_c, Nxy_c, Nx_c );
	vel_point2grid( Vel_s, Velgrid_s, Nxyz_c, Nxy_c, Nx_c );
	free(Vel_p);
	free(Vel_s);
/* Default to P-wave velocity */
	Vel_grid = Velgrid_p;

	return 0;
}

/**
 * @brief
 *
 */
void rt_velmod_free( void )
{
	free(Lat_c);
	free(Lon_c);
	free(Dep_c);
	free(Velgrid_p);
	free(Velgrid_s);
	free(Ilon_c);
	free(Ilat_c);
	free(Idep_c);

	return;
}

/**
 * @brief
 *
 * @param ray
 * @param np
 * @param tko
 * @param azi
 */
void rt_tko_azi_cal( const RAY_INFO *ray, const int np, double *tko, double *azi )
{
	double _tko, _azi;
	double tmp1, tmp2;
	double x1, y1, z1;
	double x2, y2, z2;

/* */
	if ( np < 2 )
		return;
/* Get coordinate of epc. */
	POLAR2CARTESIAN( ray[0].r, ray[0].a, ray[0].b, x1, y1, z1 );
/* Get coordinate of next point of ray */
	POLAR2CARTESIAN( ray[1].r, ray[1].a, ray[1].b, x2, y2, z2 );
	x2 -= x1;
	y2 -= y1;
	z2 -= z1;
	tmp1 = sqrt(x2 * x2 + y2 * y2 + z2 * z2 + RT_EPS);
	x1 = -x1;
	y1 = -y1;
	z1 = -z1;
/* */
	_tko = acos((x1 * x2 + y1 * y2 + z1 * z2) / (tmp1 * ray[0].r));
/* */
	tmp1 = ray[1].a - ray[0].a;
	tmp2 = ray[1].b - ray[0].b;
	tmp1 = tmp1 < 0.0 ? tan(fabs(tmp1)) : -tan(fabs(tmp1));
	tmp2 = tmp2 > 0.0 ? tan(fabs(tmp2)) : -tan(fabs(tmp2));
	_azi = atan2(tmp2 * sin(ray[0].a), tmp1);
	if ( _azi < 0.0 )
		_azi += RT_PI2;
/* */
	*tko = _tko * RT_RAD2DEG;
	*azi = _azi * RT_RAD2DEG;

	return;
}

/**
 * @brief
 *
 * @param ray
 * @param np
 * @param dx
 * @param dy
 * @param dz
 */
void rt_drvt_cal( const RAY_INFO *ray, const int np, double *dx, double *dy, double *dz )
{
	double x1, y1, z1;
	double x2, y2, z2;
	double us;

/* */
	if ( np < 2 )
		return;
/* Get coordinate of epc. */
	POLAR2CARTESIAN( ray[0].r, ray[0].a, ray[0].b, x1, y1, z1 );
/* Get coordinate of next point of ray */
	POLAR2CARTESIAN( ray[1].r, ray[1].a, ray[1].b, x2, y2, z2 );
/* */
	*dx = x2 - x1;
	*dy = y2 - y1;
	*dz = z2 - z1;
/* */
	us = -1.0 / ray[0].v;
	us /= sqrt((*dx) * (*dx) + (*dy) * (*dy) + (*dz) * (*dz));
/* */
	*dx *= us;
	*dy *= us;
	*dz *= us;

	return;
}

/**
 * @brief Return coordinates to the origin
 *
 * @param ray
 * @param np
 * @return RAY_INFO*
 */
RAY_INFO *rt_ray_rad2deg( RAY_INFO *ray, const int np, double stlo, double evlo )
{
	RAY_INFO *ray_now  = ray;
	double    shiftlon = 0.0;
	double    diff;

/* */
	if ( stlo < 0.0 )
		stlo += 360.0;
	if ( evlo < 0.0 )
		evlo += 360.0;
/* */
	if ( (diff = fabs(evlo - stlo)) < 180.0 ) {
		if ( evlo < stlo )
			shiftlon = evlo - (180.0 - diff) * 0.5;
		else
			shiftlon = stlo - (180.0 - diff) * 0.5;
	}
	else {
		diff = 360.0 - diff;
		if ( evlo < stlo )
			shiftlon = evlo - (90.0 + diff * 0.5);
		else
			shiftlon = stlo - (90.0 + diff * 0.5);
	}
/* */
	for ( int i = 0; i < np; i++, ray_now++ ) {
	/* */
		ray_now->r  = Ro - ray_now->r;
		ray_now->a *= RT_RAD2DEG;
		ray_now->b *= RT_RAD2DEG;
	/* */
		ray_now->a = geoc2geog( 90.0 - ray_now->a );
		if ( (ray_now->b += shiftlon) < 0.0 )
			ray_now->b += 360.0;
	}

	return ray;
}

/**
 * @brief
 *
 * @param evla
 * @param evlo
 * @param evdp
 * @param stla
 * @param stlo
 * @param stel
 * @param record
 * @param np
 * @param tk
 */
static void raytracing_pb( double evla, double evlo, double evdp, double stla, double stlo, double stel, RAY_INFO *ray, int *np, double *tk )
{
	int ni, i, j, k, l;

	RAY_INFO _ray[RT_MAX_NODE + 1];

	RAY_INFO *ray_now = NULL;
	RAY_INFO *ray_prev = NULL;
	RAY_INFO *ray_next = NULL;
	RAY_INFO *ray_end = NULL;
	RAY_INFO ray_mid, calc_tmp;

	NODE_INFO terminal[2];

	double xfac;
	double shiftlo;
	double to, tn;

	double acosa, sina, rsina;

	double x1, y1, z1;
	double x2, y2, z2;
	double x3, y3, z3;

	double dn, dr, da, db;
	double dseg, ddseg;
	double upz, dwz;
	double vr, vb, va;
	double rvs, cc, rcur;

/* Parameters for calculation initialization */
	xfac = XFAC;
/* ni : number of ray segments */
	ni = N1;
/* Random algorism used */
	//srand(time(NULL));
/*
 * Longitude and latitude range from 0 to 180. This program does not work with angles greater than 180.
 * Pass from latitude to colatitude
 */
	evla = (90.0 - evla) * RT_DEG2RAD;
	stla = (90.0 - stla) * RT_DEG2RAD;
/* da = bre, db = bso, dr = dlo */
	if ( stlo < 0.0 )
		da = 360.0 + stlo;
	else
		da = stlo;

	if ( evlo < 0.0 )
		db = 360.0 + evlo;
	else
		db = evlo;
/* */
	dr = fabs(db - da);
	shiftlo = 0.0;
	if ( dr < 180.0 ) {
		if ( db < da ) {
			evlo = (180.0 - dr) * 0.5;
			stlo = evlo + dr;
			shiftlo = db - evlo;
		}
		else {
			stlo = (180.0 - dr) * 0.5;
			evlo = stlo + dr;
			shiftlo = da - stlo;
		}
	}
	else {
		dr = 360.0 - dr;
		if ( db < da ) {
			evlo = 90.0 + 0.5 * dr;
			stlo = evlo - dr;
			shiftlo = db - evlo;
		}
		else {
			stlo = 90.0 + 0.5 * dr;
			evlo = stlo - dr;
			shiftlo = da - stlo;
		}
	}
/* */
	evlo *= RT_DEG2RAD;
	stlo *= RT_DEG2RAD;
/* */
	evdp = Ro - evdp;
	stel = Ro - stel;
/* Initial straight ray */
/* Epc. coordinates */
	POLAR2CARTESIAN( evdp, evla, evlo, terminal[0].x, terminal[0].y, terminal[0].z );
/* Rec. coordinates */
	POLAR2CARTESIAN( stel, stla, stlo, terminal[1].x, terminal[1].y, terminal[1].z );

	dr = (terminal[1].x - terminal[0].x) / ni;
	da = (terminal[1].y - terminal[0].y) / ni;
	db = (terminal[1].z - terminal[0].z) / ni;
	dseg = sqrt(dr * dr + da * da + db * db + RT_EPS);

	ray_now = _ray;
	ray_end = _ray + ni;

	do {
		i = ray_now - _ray;

		x1 = terminal[0].x + dr * i;
		y1 = terminal[0].y + da * i;
		z1 = terminal[0].z + db * i;

		dn = x1 * x1 + y1 * y1 + z1 * z1;
		ray_now->r = sqrt(dn + RT_EPS);

		acosa = z1 / ray_now->r;
		if ( acosa <= -1.0 )
			ray_now->a = RT_PI;
		else if ( acosa >= 1.0 )
			ray_now->a = 0.0;
		else
			ray_now->a = acos(acosa);
/* */
		sina = sin(ray_now->a);

		acosa = x1 / (sina * ray_now->r);
		if ( acosa <= -1.0 )
			ray_now->b = RT_PI;
		else if ( acosa >= 1.0 )
			ray_now->b = 0.0;
		else
			ray_now->b = acos(acosa);
		if ( y1 < 0.0 )
			ray_now->b = RT_PI2 - ray_now->b;
/* */
		ray_now->v = get_vel_ray( ray_now, shiftlo );
	} while ( ray_now++ < ray_end );
/* */
	tn = get_ray_traveltime( _ray, ni );

/* interation loop */
	while ( ni <= N2 ) {
		xfac = XFAC;
		l = ni - 1;
		for ( k = 0; k < NLOOP; k++ ) {
			if ( ni > 2 || k == 0 ) {
				for( j = 0; j < l; j++ ) {
				/* See Um & Thurber (1987) p.974. */
					if ( !(j & 0x01) )
						ray_now = _ray + (j >> 1) + 1;
					else
						ray_now = ray_end - ((j + 1) >> 1);
				/* */
					ray_prev = ray_now - 1;
					ray_next = ray_now + 1;
				/* Check if the first point of ray */
					if ( ray_prev > _ray ) {
						POLAR2CARTESIAN( ray_prev->r, ray_prev->a, ray_prev->b, x1, y1, z1 );
					}
					else {
						x1 = terminal[0].x;
						y1 = terminal[0].y;
						z1 = terminal[0].z;
					}
				/* Check if the last point of ray */
					if ( ray_next < ray_end ) {
						POLAR2CARTESIAN( ray_next->r, ray_next->a, ray_next->b, x3, y3, z3 );
					}
					else {
						x3 = terminal[1].x;
						y3 = terminal[1].y;
						z3 = terminal[1].z;
					}

				/*
				 * Check if the first time of loop & the first two perturbation,
				 * if ture, we can just skip the calculation of cartesian coordinates.
				 */
					if ( k > 0 || j >= 2 ) {
						x2 = x1 + x3;
						y2 = y1 + y3;
						z2 = z1 + z3;
					/* */
						ray_mid.r = sqrt(x2 * x2 + y2 * y2 + z2 * z2 + RT_EPS);
						acosa = z2 / ray_mid.r;
						if ( acosa <= -1.0 )
							ray_mid.a = RT_PI;
						else if ( acosa >= 1.0 )
							ray_mid.a = 0.0;
						else
							ray_mid.a = acos(acosa);
					/* */
						sina = sin(ray_mid.a);
						rsina = sina * ray_mid.r;
					/* cosa = cos(ray_mid.a); */
						acosa = x2 / rsina;
						if ( acosa <= -1.0 )
							ray_mid.b = RT_PI;
						else if ( acosa >= 1.0 )
							ray_mid.b = 0.0;
						else
							ray_mid.b = acos(acosa);
						if ( y2 < 0.0 )
							ray_mid.b = RT_PI2 - ray_mid.b;
					/* */
						ray_mid.r *= 0.5;
						rsina *= 0.5;

					/* Determine velocity at 3 points */
						/* v1 = ray_prev->v; */
						ray_mid.v = get_vel_ray( &ray_mid, shiftlo );
						/* v3 = ray_next->v; */
					}
					else {
						ray_mid = *ray_now;
						sina = sin(ray_mid.a);
						rsina = sina * ray_mid.r;
					}
				/* */
					dr = x3 - x1;
					da = y3 - y1;
					db = z3 - z1;
					dn = dr * dr + da * da + db * db;
					dseg = sqrt(dn + RT_EPS);
					dr = (ray_next->r - ray_prev->r); /* Orig. dr = (ray_next->r - ray_prev->r)/dseg */
					da = (ray_next->a - ray_prev->a); /* Orig. da = (ray_next->a - ray_prev->a)/dseg */
					db = (ray_next->b - ray_prev->b); /* Orig. db = (ray_next->b - ray_prev->b)/dseg */
				/*
				 * Now ddseg will be a distance to find dV along the coordinates
				 * Begin find the gradients and velocities first find the length of segment
				 */
					dseg *= 0.5;
					ddseg = dseg * 0.5;
				/* Begin to determine coordinates of pints surroundibg point a2, b2, r2 at the distance ddseg */
					upz = ray_mid.r + ddseg;
					dwz = ray_mid.r - ddseg;
					if ( upz > Rs ) {
					/* I guess it should be Ro + 10.0(Rs) */
						upz = Rs;
						dwz = upz - dseg;
					}
					if ( dwz <= RT_EPS ) {
						dwz = RT_EPS;
						upz = dwz + dseg;
					/* Set to Ro, the old mistake? */
						//upz = Ro;
					}

				/*
				 * The following if-endif is just for P & S, thus comment out for SKS & PKP!
				 * This gives the lowermost mantle Vp in the outer core
				 */
					calc_tmp.a = ray_mid.a;
					calc_tmp.b = ray_mid.b;
					calc_tmp.r = upz;
					vr = get_vel_ray( &calc_tmp, shiftlo );
					calc_tmp.r = dwz;
					vr -= get_vel_ray( &calc_tmp, shiftlo );
					//vr = calc_tmp.v; /* Orig. vr = calc_tmp.v/dseg */
					step_ray_node( &calc_tmp, &ray_mid, ddseg, RT_RNULL );
					vb = get_vel_ray( &calc_tmp, shiftlo );
					step_ray_node( &calc_tmp, &ray_mid, -ddseg, RT_RNULL );
					vb -= get_vel_ray( &calc_tmp, shiftlo );
					//vb = calc_tmp.v; /* Orig. vb = calc_tmp.v/dseg */
					step_ray_node( &calc_tmp, &ray_mid, RT_RNULL, ddseg );
					va = get_vel_ray( &calc_tmp, shiftlo );
					step_ray_node( &calc_tmp, &ray_mid, RT_RNULL, -ddseg );
					va -= get_vel_ray( &calc_tmp, shiftlo );
					//va = calc_tmp.v; /* Orig. va = calc_tmp.v/dseg */
				/*
				 * spherical velocity gradient:
				 * va = va / r2
				 * vb = vb / r2 / sina
				 * (tangential vector) = (slowness vector) / s
				 */
				/* dr *= 1 */
					da *= ray_mid.r;
					db *= rsina;
					rcur = dr * vr + da * va + db * vb;
					rcur /= dn;
				/* */
					dr = vr - rcur * dr;
					da = va - rcur * da;
					db = vb - rcur * db;
					rvs = dr * dr + da * da + db * db;
				/* */
					if ( rvs <= RT_EPS ) {
					/* Special condition: velocity gradient equal to zero */
						*ray_now = ray_mid;
					}
					else {
						rvs = 1.0 / sqrt(rvs); /* sqrt(rvs + RT_EPS) */
						cc = 0.5 / ray_prev->v + 0.5 / ray_next->v; /* cc = (1.0/v1 + 1.0/v3)/2.0 */
						rcur = vr * dr + va * da + vb * db;
						rcur /= dseg;
						rcur *= rvs;
					/*
					 * Tut esli rcur < 0.0 proishodit hernia poetomu postavlen abs. Ne yasno mozhno li eto delat
					 * ili net no rabotaet. Obichno oshibka poyavliaetsia ochen redko v nekotorih tochkah v etom sluchae
					 * abs prosto ne daet oshibki y posledniaya iteraciya uzhe ne imeet rcur negativnim
					 * y podgoniaet normalno reshenie (mozhet bit)
					 */
						if ( rcur < 0.0 ) {
							fprintf(stderr, "raytracing_pb: Got negative Rc!\n");
							rcur = fabs(rcur);
						}
					/* */
						ray_mid.v *= cc;
						rcur = -(ray_mid.v + 1.0) / (4.0 * cc * rcur);
						rcur += sqrt(rcur * rcur + dn / (8.0 * ray_mid.v) + RT_EPS);
						rcur *= rvs;
					/* */
						dr *= rcur;
						rcur /= rsina;
						da *= rcur * sina;
						db *= rcur;
						dr += ray_mid.r;
						da += ray_mid.a;
						db += ray_mid.b;
					/* if ray_now->r > 6371 then force it to the surface. */
						if ( (ray_now->r = (dr - ray_now->r) * xfac + ray_now->r) > Rs )
							ray_now->r = Rs;
						ray_now->a = (da - ray_now->a) * xfac + ray_now->a;
						ray_now->b = (db - ray_now->b) * xfac + ray_now->b;
						ray_now->v = get_vel_ray( ray_now, shiftlo );
					}
				}
			}
			else {
			/* ray_now = ray + (ni >> 1); */
			/* if ray_now->r > 6371 then force it to the surface. */
				if ( (ray_now->r = (dr - ray_now->r) * xfac + ray_now->r) > Rs )
					ray_now->r = Rs;
				ray_now->a = (da - ray_now->a) * xfac + ray_now->a;
				ray_now->b = (db - ray_now->b) * xfac + ray_now->b;
				ray_now->v = get_vel_ray( ray_now, shiftlo );
			}
/* */
			to = tn;
			tn = get_ray_traveltime( _ray, ni );
			if ( fabs(to - tn) <= to * FLIMIT )
				break;
		/* Random algorism, it will be faster but the result is not stable! */
			//xfac = (double)(rand()%10 + 10)/10.0;
			if ( (xfac *= 0.98) < 1.0 )
				xfac = 1.0;
			/* printf("%d %lf\n", k, xfac); */
		}

	/* Skip increasing of segment number if minimum length of segment is exceed or maximum number of segments was reached */
		if ( dseg < MINS || ni >= N2 )
			break;
		/* igood = 1; */

	/* Double the number of points. */
		ray_next = ray_end + 1;
		memcpy(ray_next, _ray + 1, sizeof(RAY_INFO) * ni);
		for ( i = 0; i < ni; i++ )
			_ray[(i + 1) << 1] = ray_next[i];
	/* */
		ni <<= 1;
		ray_end = _ray + ni;
	/* */
		x1 = terminal[0].x;
		y1 = terminal[0].y;
		z1 = terminal[0].z;
		ray_now = _ray + 1;
		ray_next = ray_now + 1;
		do {
			if ( ray_next < ray_end ) {
				POLAR2CARTESIAN( ray_next->r, ray_next->a, ray_next->b, x3, y3, z3 );
			}
			else {
				x3 = terminal[1].x;
				y3 = terminal[1].y;
				z3 = terminal[1].z;
			}
		/* */
			x2 = x3 + x1;
			y2 = y3 + y1;
			z2 = z3 + z1;
			ray_now->r = sqrt(x2 * x2 + y2 * y2 + z2 * z2 + RT_EPS);
		/* */
			acosa = z2 / ray_now->r;
			if ( acosa <= -1.0 )
				ray_now->a = RT_PI;
			else if ( acosa >= 1.0 )
				ray_now->a = 0.0;
			else
				ray_now->a = acos(acosa);
		/* */
			sina = sin(ray_now->a);
			acosa = x2 / (sina * ray_now->r);
			if ( acosa <= -1.0 )
				ray_now->b = RT_PI;
			else if ( acosa >= 1.0 )
				ray_now->b = 0.0;
			else
				ray_now->b = acos(acosa);
			if ( y2 < 0.0 )
				ray_now->b = RT_PI2 - ray_now->b;
		/* */
			ray_now->r *= 0.5;
			ray_now->v = get_vel_ray( ray_now, shiftlo );
		/* */
			x1 = x3;
			y1 = y3;
			z1 = z3;
			ray_now += 2;
			ray_next += 2;
		} while ( ray_next <= ray_end );
	/* */
		to = tn;
		tn = get_ray_traveltime( _ray, ni );
	/* */
		if ( fabs(to - tn) <= to * FLIMIT )
			break;
		/* igood = 1; */
	}

/* Output the result */
	*tk = tn;
	*np = ni + 1;
	memcpy(ray, _ray, sizeof(RAY_INFO) * (*np));

	return;
}

/**
 * @brief This subroutine calculate position of new point in polar coordinates basing on the coordinates
 *        of main point in radians (la is colatitude) and dx and dy in kilometers
 *
 * @param new
 * @param origin
 * @param dx
 * @param dy
 */
static void step_ray_node( RAY_INFO *new, const RAY_INFO *origin, const double dx, const double dy )
{
	const int flag_x = fabs(dx - RT_RNULL) > RT_RNULL;
	const int flag_y = fabs(dy - RT_RNULL) > RT_RNULL;

/* */
	new->b = origin->b;
	new->a = origin->a;
/* */
	if ( flag_y ) {
		new->a += atan2(dy, origin->r);
		if ( new->a > RT_PI ) {
			new->a = RT_PI2 - new->a;
			new->b += RT_PI;
		}
		if ( new->a < 0.0 ) {
			new->a = fabs(new->a);
			new->b += RT_PI;
		}
	}
/* */
	if ( flag_x ) {
		new->b += atan2(dx, origin->r * sin(origin->a));
		if ( new->b < 0.0 )
			new->b += RT_PI2;
	}
	if ( new->b > RT_PI2 )
		new->b -= RT_PI2;

	new->r = sqrt(origin->r * origin->r + (flag_x ? dx * dx : 0.0) + (flag_y ? dy * dy : 0.0) + RT_EPS);

	return;
}

/**
 * @brief
 *
 * @param ray
 * @param m
 * @return double
 */
static double get_ray_traveltime( const RAY_INFO *ray, const int np )
{
	const RAY_INFO * const ray_end = ray + np;

	double r1, r2;
	double rv1, rv2;
	double r1sq, r2sq;
	double sin1, sin2;
	double cos1, cos2, cosa;
	double result = 0.0;

	r1   = ray->r;
	rv1  = 0.5 / ray->v; /* 1.0/v[0]/2.0 */
	r1sq = r1 * r1;
	sin1 = sin(ray->a);
	cos1 = cos(ray->a);
	cosa = ray->b;

	while ( ray++ < ray_end ) {
		r2   = ray->r;
		rv2  = 0.5 / ray->v; /* 1.0/v[i]/2.0 */
		r2sq = r2 * r2;
		sin2 = sin(ray->a);
		cos2 = cos(ray->a);
		cosa = sin1 * sin2 * cos(cosa - ray->b) + cos1 * cos2;

		result += sqrt(r1sq + r2sq - 2.0 * r1 * r2 * cosa + RT_EPS) * (rv1 + rv2);

		r1   = r2;
		rv1  = rv2;
		r1sq = r2sq;
		sin1 = sin2;
		cos1 = cos2;
		cosa = ray->b;
	}

	return result;
}

/**
 * @brief
 *
 * @param ray
 * @param shiftlon
 * @return double
 */
static double get_vel_ray( const RAY_INFO *ray, const double shiftlon )
{
	double lat, lon, dep;

	lat = geoc2geog(90.0 - ray->a * RT_RAD2DEG);
	lon = ray->b * RT_RAD2DEG + shiftlon;
	dep = Ro - ray->r;

	return get_vel_geog(lon, lat, dep);
}

/**
 * @brief
 *
 * @param lon
 * @param lat
 * @param dep
 * @return double
 */
static double get_vel_geog( const double lon, const double lat, const double dep )
{
	double lonf, latf, depf;
	double wv[8];
	int ip, jp, kp;

	VEL_GRID velg;

	ip = (int)(lon * Ibld3);
	jp = (int)(lat * Ibld3);
	kp = (int)(dep * Ibld4);

/* If the ray point out of range, give it the center velocity */
	if ( intmap_3d(&ip, &jp, &kp) )
		return RT_EPS;
/* */
	lonf = Lon_c[ip];
	lonf = (lon - lonf) / (Lon_c[ip + 1] - lonf);
	latf = Lat_c[jp];
	latf = (lat - latf) / (Lat_c[jp + 1] - latf);
	depf = Dep_c[kp];
	depf = (dep - depf) / (Dep_c[kp + 1] - depf);

	//memcpy(&velg, Vel_grid + ip + jp*Nx_c + kp*Nxy_c, sizeof(VEL_GRID));
	velg = *(Vel_grid + ip + jp * Nx_c + kp * Nxy_c);

/*
 * lonf1 = 1.0 - lonf
 * latf1 = 1.0 - latf
 * depf1 = 1.0 - depf
 */
	wv[6] = wv[7] = lonf * latf;         /* lonf*latf; 1,1,0 & 1,1,1 */
	wv[4] = wv[5] = lonf - wv[7];        /* lonf*latf1; 1,0,0 & 1,0,1 */
	wv[2] = wv[3] = latf - wv[7];        /* lonf1*latf; 0,1,0 & 0,1,1 */
	wv[0] = wv[1] = 1.0 - lonf - wv[3];  /* lonf1*latf1; 0,0,0 & 0,0,1 */

	wv[0] *= velg.vel[0];
	wv[2] *= velg.vel[2];
	wv[4] *= velg.vel[4];
	wv[6] *= velg.vel[6];
	wv[0] += wv[2] + wv[4] + wv[6];

	if ( depf <= RT_EPS ) {
		wv[1] = 0.0;
	}
	else {
		wv[1] *= velg.vel[1];
		wv[3] *= velg.vel[3];
		wv[5] *= velg.vel[5];
		wv[7] *= velg.vel[7];

		wv[1] += wv[3] + wv[5] + wv[7];
		wv[1] = wv[0] - wv[1];
		wv[1] *= depf;

		wv[0] -= wv[1];
	}

	return wv[0];
}

/**
 * @brief
 *
 * @param ip
 * @param jp
 * @param kp
 */
static int intmap_3d( int *ip, int *jp, int *kp )
{
	const int lon = *ip;
	const int lat = *jp;
	const int dep = *kp;

/*
 * *ip = (int)((lon + Lon1_c)/Bld3-1.0);
 * *jp = (int)((lat + Lat1_c)/Bld3-1.0);
 * *kp = (int)((dep + Dep1_c)/Bld4-1.0);
 */
	*ip += Ilon1_c;
	*jp += Ilat1_c;
	*kp += Idep1_c;
/* Checking process */
	if ( *ip < 0 || *jp < 0 || *kp < 0 ||
		*ip >= (Ilonmax - 2) || *jp >= (Ilatmax - 2) || *kp >= (Idepmax - 2)
	) {
		fprintf(stderr, "intmap_3d: NOTICE!!! lon, lat or dep out of range! Force back to the nearest bolder.\n" );
		fprintf(stderr, "lon = %lf, lat = %lf, dep = %lf\n", (double)lon / Ibld3, (double)lat / Ibld3, (double)dep / Ibld4);
		fprintf(stderr, "ip = %d, jp = %d, kp = %d\n", *ip, *jp, *kp);
		return -1;
	}
/* */
	*ip = Ilon_c[*ip];
	*jp = Ilat_c[*jp];
	*kp = Idep_c[*kp];

	return 0;
}

/**
 * @brief
 *
 * @param modelfile
 * @return int
 */
static int input_vel_model( const char *modelfile )
{
	double *ptrtmp = NULL;
	double bld3, bld4;
	double average;
	char fracline[1024] = { 0 };
	FILE *fp = NULL;

/* */
	if ( (fp = fopen(modelfile, "r")) == NULL ) {
		fprintf(stderr, "input_vel_model: Opening %s file ERROR; exiting!\n", modelfile);
		return -1;
	}
	if ( fscanf(fp, "%lf %lf %d %d %d\n", &bld3, &bld4, &Nlon_c, &Nlat_c, &Ndep_c) != 5 ) {
		fprintf(stderr, "input_vel_model: Reading VpVs Model header ERROR; exiting!\n" );
		return -1;
	}
/* */
	Lon_c = calloc(Nlon_c, sizeof(double));
	Lat_c = calloc(Nlat_c, sizeof(double));
	Dep_c = calloc(Ndep_c, sizeof(double));

	Nx_c = Nlon_c;
	Nxy_c = Nx_c * Nlat_c;
	Nxyz_c = Nxy_c * Ndep_c;

	Vel_p = calloc(Nxyz_c, sizeof(double));
	Vel_s = calloc(Nxyz_c, sizeof(double));
/* */
	if ( fgets(fracline, sizeof(fracline) - 1, fp) != NULL ) {
		ptrtmp = Lon_c;
		for ( int i = 0; i < Nlon_c; i++ ) {
			if ( i < Nlon_c - 1 ) {
				if ( sscanf(fracline, " %lf %[^\n]", ptrtmp, fracline) != 2 ) {
					fprintf(stderr, "input_vel_model: Reading VpVs Model lon_c ERROR; exiting!\n" );
					return -1;
				}
			}
			else {
				if ( sscanf(fracline, "%lf", ptrtmp) != 1 ) {
					fprintf(stderr, "input_vel_model: Reading VpVs Model lon_c ERROR; exiting!\n" );
					return -1;
				}
			}
			ptrtmp++;
		}
	}
/* */
	if ( fgets( fracline, sizeof(fracline) - 1, fp ) != NULL ) {
		ptrtmp = Lat_c;
		for ( int i = 0; i < Nlat_c; i++ ) {
			if ( i < Nlat_c - 1 ) {
				if ( sscanf( fracline, "%lf %[^\n]", ptrtmp, fracline ) != 2 ) {
					printf("Reading VpVs Model lat_c error; exiting!\n");
					return -1;
				}
			}
			else {
				if ( sscanf( fracline, "%lf", ptrtmp ) != 1 ) {
					fprintf(stderr, "input_vel_model: Reading VpVs Model lat_c ERROR; exiting!\n" );
					return -1;
				}
			}
			ptrtmp++;
		}
	}
/* */
	if ( fgets( fracline, sizeof(fracline) - 1, fp ) != NULL ) {
		ptrtmp = Dep_c;
		for ( int i = 0; i < Ndep_c; i++ ) {
			if ( i < Ndep_c - 1 ) {
				if ( sscanf( fracline, "%lf %[^\n]", ptrtmp, fracline ) != 2 ) {
					fprintf(stderr, "input_vel_model: Reading VpVs Model dep_c ERROR; exiting!\n" );
					return -1;
				}
			}
			else {
				if ( sscanf( fracline, "%lf", ptrtmp ) != 1 ) {
					fprintf(stderr, "input_vel_model: Reading VpVs Model dep_c ERROR; exiting!\n" );
					return -1;
				}
			}
			ptrtmp++;
		}
	}
/* Read P velocity model */
	ptrtmp = Vel_p;
	for ( int k = 0; k < Ndep_c; k++ ) {
		for ( int j = 0; j < Nlat_c; j++ ) {
			if ( fgets(fracline, sizeof(fracline) - 1, fp) != NULL ) {
				for ( int i = 0; i < Nlon_c; i++ ) {
					if ( i < Nlon_c - 1 ) {
						if ( sscanf(fracline, "%lf %[^\n]", ptrtmp, fracline) != 2 ) {
							fprintf(stderr, "input_vel_model: Reading VpVs Model Vel_p ERROR; exiting!\n" );
							return -1;
						}
					}
					else {
						if ( sscanf(fracline, "%lf", ptrtmp) != 1 ) {
							fprintf(stderr, "input_vel_model: Reading VpVs Model Vel_p ERROR; exiting!\n" );
							return -1;
						}
					}
					ptrtmp++;
				}
			}
		}
	}
/* Read S velocity model */
	ptrtmp = Vel_s;
	for ( int k = 0; k < Ndep_c; k++ ) {
		for ( int j = 0; j < Nlat_c; j++ ) {
			if ( fgets( fracline, sizeof(fracline) - 1, fp ) != NULL ) {
				for ( int i = 0; i < Nlon_c; i++ ) {
					if ( i < Nlon_c - 1 ) {
						if ( sscanf( fracline, "%lf %[^\n]", ptrtmp, fracline ) != 2 ) {
							fprintf(stderr, "input_vel_model: Reading VpVs Model Vel_s ERROR; exiting!\n" );
							return -1;
						}
					}
					else {
						if ( sscanf( fracline, "%lf", ptrtmp ) != 1 ) {
							fprintf(stderr, "input_vel_model: Reading VpVs Model Vel_s ERROR; exiting!\n" );
							return -1;
						}
					}
					ptrtmp++;
				}
			}
		}
	}
/* */
	fclose(fp);
/* */
	if ( bldmap( bld3, bld4 ) )
		return -1;

/*
 * Nx2_c = Nlon_c - 2;
 * Nxy2_c = Nx2_c * (Nlat_c-2);
 * Nxyz2_c = Nxy2_c * (Ndep_c-2);
 */
	average = 0.0;
	for ( int i = 0; i < Nlat_c; i++ )
		average += Lat_c[i];
	average /= (double)Nlat_c;
	Ro = get_earth_radius(average);
	Rs = Ro + 10.0;

	return 0;
}

/**
 * @brief
 *
 * @param bld3
 * @param bld4
 * @return int
 */
static int bldmap( const double bld3, const double bld4 )
{
	int matrix_index;
/* For crustal velocity */
	const double lon1 = -Lon_c[0];
	const double lat1 = -Lat_c[0];
	const double dep1 = -Dep_c[0];

	Ilonmax = (int)(RT_EPS + (Lon_c[Nlon_c-1] + lon1) / bld3);
	Ilatmax = (int)(RT_EPS + (Lat_c[Nlat_c-1] + lat1) / bld3);
	Idepmax = (int)(RT_EPS + (Dep_c[Ndep_c-1] + dep1) / bld4);

	if ( Ilonmax > RT_MAX_NODE || Ilatmax > RT_MAX_NODE || Idepmax > RT_MAX_NODE ) {
		fprintf(stderr, "bldmap: ERROR!! Model dimension too big (Max. %d)!\n", RT_MAX_NODE);
		return -1;
	}

	Ilon_c = calloc(Ilonmax, sizeof(int));
	Ilat_c = calloc(Ilatmax, sizeof(int));
	Idep_c = calloc(Idepmax, sizeof(int));

	matrix_index = 0;
	for ( int i = 0; i < Ilonmax; i++ ) {
		if ( (i * bld3 - lon1) >= Lon_c[matrix_index + 1] )
			matrix_index++;
		Ilon_c[i] = matrix_index;
	}

	matrix_index = 0;
	for ( int i = 0; i < Ilatmax; i++ ) {
		if ( (i * bld3 - lat1) >= Lat_c[matrix_index + 1] )
			matrix_index++;
		Ilat_c[i] = matrix_index;
	}

	matrix_index = 0;
	for ( int i = 0; i < Idepmax; i++ ) {
		if ( (i * bld4 - dep1) >= Dep_c[matrix_index + 1] )
			matrix_index++;
		Idep_c[i] = matrix_index;
	}

/* Change coordinates to integer */
	Ilon1_c = (int)(lon1 / bld3);
	Ilat1_c = (int)(lat1 / bld3);
	Idep1_c = (int)(dep1 / bld4 + RT_EPS);
	Ibld3 = (int)(1.0 / bld3 + RT_EPS);
	Ibld4 = (int)(1.0 / bld4 + RT_EPS);

	return 0;
}

/**
 * @brief
 *
 * @param vel_p
 * @param vel_g
 * @param nxyz
 * @param nxy
 * @param nx
 * @return int
 */
static int vel_point2grid( const double *vel_p, VEL_GRID *vel_g, const int nxyz, const int nxy, const int nx )
{
	const double *velptr;
	const double *velprt_end = vel_p + nxyz;
	VEL_GRID *velgptr;

/* */
	for ( int i = 0; i < nxyz; i++ ) {
		velptr = vel_p + i;
		velgptr = vel_g + i;

		velgptr->vel[0] = *(velptr);
		if ( velptr + nx > velprt_end )
			velgptr->vel[2] = *(velptr);
		else
			velgptr->vel[2] = *(velptr + nx);
		if ( velptr + 1 > velprt_end )
			velgptr->vel[4] = *(velptr);
		else
			velgptr->vel[4] = *(velptr + 1);
		if ( velptr + 1 + nx > velprt_end )
			velgptr->vel[6] = *(velptr);
		else
			velgptr->vel[6] = *(velptr + 1 + nx);

		if ( velptr + nxy > velprt_end ) {
			velgptr->vel[1] = velgptr->vel[0];
			velgptr->vel[3] = velgptr->vel[2];
			velgptr->vel[5] = velgptr->vel[4];
			velgptr->vel[7] = velgptr->vel[6];
		}
		else {
			velptr += nxy;
			velgptr->vel[1] = *(velptr);
			if ( velptr + nx > velprt_end )
				velgptr->vel[3] = *(velptr);
			else
				velgptr->vel[3] = *(velptr + nx);
			if ( velptr + 1 > velprt_end )
				velgptr->vel[5] = *(velptr);
			else
				velgptr->vel[5] = *(velptr + 1);
			if ( velptr + 1 + nx > velprt_end )
				velgptr->vel[7] = *(velptr);
			else
				velgptr->vel[7] = *(velptr + 1 + nx);
		}
	}

	return 0;
}

/**
 * @brief
 *
 * @param xlat
 * @return double
 */
static double geog2geoc( const double xlat )
{
	return atan(RT_B2A_SQ * tan(xlat * RT_DEG2RAD)) * RT_RAD2DEG;
}

/**
 * @brief
 *
 * @param xlat
 * @return double
 */
static double geoc2geog( const double xlat )
{
	return atan(tan(xlat * RT_DEG2RAD) / RT_B2A_SQ) * RT_RAD2DEG;
}

/**
 * @brief
 *
 * @param xlat
 * @return double
 */
static double get_earth_radius( double xlat )
{
	double dlt1;
	const double re = 6378.163;
	const double ell = 298.26;

/*
 * xlat = xlat * RT_DEG2RAD;
 * dlt1 = atan(0.99330647 * tan(xlat));
 */
	xlat *= RT_DEG2RAD;
/* conversion factor for latitude, SDC. conversion factor is 0.99330647 */
	dlt1 = atan(0.99330647 * tan(xlat));
/* dlt1 = atan(RT_B2A_SQ * tan(xlat * RT_DEG2RAD / 60.0)); */
	dlt1 = sin(dlt1);
	dlt1 *= dlt1;

	return re * (1.0 - dlt1 / ell);
}
