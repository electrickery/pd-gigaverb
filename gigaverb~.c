/*
  Copyright (C) 1999 Juhana Sadeharju
  kouhia at nic.funet.fi

  Max/MSP port 2004 by Olaf Matthes, <olaf.matthes@gmx.de>
  Pd port 2015 by Marco Matteo Markidis, <mm.markidis@gmail.com>
  Copyright (C) 2015 Marco Matteo Markidis

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include <math.h>
#include <string.h>
#include "m_pd.h"

#define FDNORDER 4

/* Convert a value in dB's to a coefficent */
#define DB_CO(g) ((g) > -90.0f ? pow(10.0f, (g) * 0.05f) : 0.0f)
/* and back to dB */
#define CO_DB(g) ((g) != 0.0f ? 20.0f/log(10) * log((g)) : -90.0f)
#define IS_DENORM_FLOAT(v) ((((*(unsigned long *)&(v))&0x7f800000)==0)&&((v)!=0.f))
#define IS_DENORM_DOUBLE(v) ((((((unsigned long *)&(v))[1])&0x7fe00000)==0)&&((v)!=0.))			

#define IS_NAN_FLOAT(v)	(((*(unsigned long *)&(v))&0x7f800000)==0x7f800000)
#define IS_NAN_DOUBLE(v) (((((unsigned long *)&(v))[1])&0x7fe00000)==0x7fe00000) 
#define FIX_DENORM_NAN_FLOAT(v)	((v)=IS_DENORM_NAN_FLOAT(v)?0.f:(v))
#define IS_DENORM_NAN_FLOAT(v)	(IS_DENORM_FLOAT(v)||IS_NAN_FLOAT(v))
#define IS_DENORM_NAN_DOUBLE(v)	(IS_DENORM_DOUBLE(v)||IS_NAN_DOUBLE(v))
#define CLIP(a, lo, hi) ( (a)>(lo)?( (a)<(hi)?(a):(hi) ):(lo) )

typedef struct {
  int size;
  int idx;
  float *buf;
} ty_fixeddelay;

typedef struct {
  int size;
  float coeff;
  int idx;
  float *buf;
} ty_diffuser;

typedef struct {
  float damping;
  float delay;
} ty_damper;

typedef struct
{
  t_object obj;
  int bypass;
  int rate;
  float inputbandwidth;
  float drylevel;
  float wetlevel;
  float taillevel;
  float earlylevel;
  ty_damper *inputdamper;
  float maxroomsize;
  float roomsize;
  float revtime;
  float maxdelay;
  float largestdelay;
  ty_fixeddelay **fdndels;
  float *fdngains;
  int *fdnlens;
  ty_damper **fdndamps; 
  float fdndamping;
  ty_diffuser **ldifs;
  ty_diffuser **rdifs;
  ty_fixeddelay *tapdelay;
  int *taps;
  float *tapgains;
  float *d;
  float *u;
  float *f;
  double alpha;
  float x_f;
} ty_gverb;

void *gigaverb_class;

ty_gverb *gigaverb_new(t_symbol *s, short argc, t_atom *argv);
void gigaverb_free(ty_gverb *);
void gigaverb_flush(ty_gverb *);
static void gigaverb_do(ty_gverb *, float, float *, float *);
static void gigaverb_set_roomsize(ty_gverb *, t_floatarg);
static void gigaverb_set_revtime(ty_gverb *, t_floatarg);
static void gigaverb_set_damping(ty_gverb *, t_floatarg);
static void gigaverb_set_inputbandwidth(ty_gverb *, t_floatarg);
static void gigaverb_set_drylevel(ty_gverb *, t_floatarg);
static void gigaverb_set_wetlevel(ty_gverb *, t_floatarg);
static void gigaverb_set_earlylevel(ty_gverb *, t_floatarg);
static void gigaverb_set_taillevel(ty_gverb *, t_floatarg);
ty_diffuser *diffuser_make(int, float);
void diffuser_free(ty_diffuser *);
void diffuser_flush(ty_diffuser *);
ty_damper *damper_make(float);
void damper_free(ty_damper *);
void damper_flush(ty_damper *);
ty_fixeddelay *fixeddelay_make(int);
void fixeddelay_free(ty_fixeddelay *);
void fixeddelay_flush(ty_fixeddelay *);
int isprime(int);
int nearest_prime(int, float);
int ff_round(float f);
int ff_trunc(float f);

static inline float diffuser_do(ty_diffuser *p, float x)
{
  float y,w;

  w = x - p->buf[p->idx]*p->coeff;
  FIX_DENORM_NAN_FLOAT(w);
  y = p->buf[p->idx] + w*p->coeff;
  p->buf[p->idx] = w;
  p->idx = (p->idx + 1) % p->size;
  return(y);
}

static inline float fixeddelay_read(ty_fixeddelay *p, int n)
{
  int i;

  i = (p->idx - n + p->size) % p->size;
  return(p->buf[i]);
}

static inline void fixeddelay_write(ty_fixeddelay *p, float x)
{
  FIX_DENORM_NAN_FLOAT(x);
  p->buf[p->idx] = x;
  p->idx = (p->idx + 1) % p->size;
}

static inline void damper_set(ty_damper *p, float damping)
{ 
  p->damping = damping;
} 
  
static inline float damper_do(ty_damper *p, float x)
{ 
  float y;

  y = x*(1.0-p->damping) + p->delay*p->damping;
  p->delay = y;
  return(y);
}

/*
 * This FDN reverb can be made smoother by setting matrix elements at the
 * diagonal and near of it to zero or nearly zero. By setting diagonals to zero
 * means we remove the effect of the parallel comb structure from the
 * reverberation.  A comb generates uniform impulse stream to the reverberation
 * impulse response, and thus it is not good. By setting near diagonal elements
 * to zero means we remove delay sequences having consequtive delays of the
 * similar lenths, when the delays are in sorted in length with respect to
 * matrix element index. The matrix described here could be generated by
 * differencing Rocchesso's circulant matrix at max diffuse value and at low
 * diffuse value (approaching parallel combs).
 *
 * Example 1:
 * Set a(k,k), for all k, equal to 0.
 *
 * Example 2:
 * Set a(k,k), a(k,k-1) and a(k,k+1) equal to 0.
 *
 * Example 3: The transition to zero gains could be smooth as well.
 * a(k,k-1) and a(k,k+1) could be 0.3, and a(k,k-2) and a(k,k+2) could
 * be 0.5, say.
 */

static inline void gigaverb_fdnmatrix(float *a, float *b)
{
  const float dl0 = a[0], dl1 = a[1], dl2 = a[2], dl3 = a[3];

  b[0] = 0.5f*(+dl0 + dl1 - dl2 - dl3);
  b[1] = 0.5f*(+dl0 - dl1 - dl2 + dl3);
  b[2] = 0.5f*(-dl0 + dl1 - dl2 + dl3);
  b[3] = 0.5f*(+dl0 + dl1 + dl2 + dl3);
}

static inline void gigaverb_do(ty_gverb *p, float x, float *yl, float *yr)
{
  float z;
  unsigned int i;
  float lsum,rsum,sum,sign;

  if(IS_NAN_FLOAT(x) || IS_DENORM_FLOAT(x) || fabsf(x) > 100000.0f) {
    x = 0.0f;
  }

  z = damper_do(p->inputdamper, x);

  z = diffuser_do(p->ldifs[0],z);

  for(i = 0; i < FDNORDER; i++) {
    p->u[i] = p->tapgains[i]*fixeddelay_read(p->tapdelay,p->taps[i]);
  }
  fixeddelay_write(p->tapdelay,z);

  for(i = 0; i < FDNORDER; i++) {
    p->d[i] = damper_do(p->fdndamps[i],
			p->fdngains[i]*fixeddelay_read(p->fdndels[i],
						       p->fdnlens[i]));
  }

  sum = 0.0f;
  sign = 1.0f;
  for(i = 0; i < FDNORDER; i++) {
    sum += sign*(p->taillevel*p->d[i] + p->earlylevel*p->u[i]);
    sign = -sign;
  }
  sum += x*p->earlylevel;
  lsum = sum;
  rsum = sum;

  gigaverb_fdnmatrix(p->d,p->f);

  for(i = 0; i < FDNORDER; i++) {
    fixeddelay_write(p->fdndels[i],p->u[i]+p->f[i]);
  }

  lsum = diffuser_do(p->ldifs[1],lsum);
  lsum = diffuser_do(p->ldifs[2],lsum);
  lsum = diffuser_do(p->ldifs[3],lsum);
  rsum = diffuser_do(p->rdifs[1],rsum);
  rsum = diffuser_do(p->rdifs[2],rsum);
  rsum = diffuser_do(p->rdifs[3],rsum);

  *yl = lsum;
  *yr = rsum;
}

static inline void gigaverb_set_roomsize(ty_gverb *p, t_floatarg a)
{
  int i;
  
  if(a <= 1.0 || IS_NAN_FLOAT(a)) {
    p->roomsize = 1.0;
  }
  else {
    p->roomsize = CLIP(a, 1.0f, p->maxroomsize);
  }
  p->largestdelay = p->rate * p->roomsize * 0.00294f;

  p->fdnlens[0] = ff_round(1.000000f*p->largestdelay);
  p->fdnlens[1] = ff_round(0.816490f*p->largestdelay);
  p->fdnlens[2] = ff_round(0.707100f*p->largestdelay);
  p->fdnlens[3] = ff_round(0.632450f*p->largestdelay);
  for(i = 0; i < FDNORDER; i++)
    {
      p->fdngains[i] = -powf((float)p->alpha, p->fdnlens[i]);
    }

  p->taps[0] = 5+ff_round(0.410f*p->largestdelay);
  p->taps[1] = 5+ff_round(0.300f*p->largestdelay);
  p->taps[2] = 5+ff_round(0.155f*p->largestdelay);
  p->taps[3] = 5+ff_round(0.000f*p->largestdelay);

  for(i = 0; i < FDNORDER; i++) {
    p->tapgains[i] = powf((float)p->alpha, p->taps[i]);
  }
}

static inline void gigaverb_set_revtime(ty_gverb *p, t_floatarg a)
{
  float ga,gt;
  double n;
  unsigned int i;

  p->revtime = CLIP(a, 0.1f, 360.0f);

  ga = 60.0f;
  gt = p->revtime;
  ga = powf(10.0f,-ga/20.0f);
  n = p->rate*gt;
  p->alpha = (double)powf(ga,1.0f/n);

  for(i = 0; i < FDNORDER; i++) {
    p->fdngains[i] = -powf((float)p->alpha, p->fdnlens[i]);
  }
}

static inline void gigaverb_set_damping(ty_gverb *p, t_floatarg a)
{
  unsigned int i;

  p->fdndamping = CLIP(a, 0.0f, 1.0f);
  for(i = 0; i < FDNORDER; i++) {
    damper_set(p->fdndamps[i],p->fdndamping);
  }
}

static inline void gigaverb_set_inputbandwidth(ty_gverb *p, t_floatarg a)
{
  p->inputbandwidth = CLIP(a, 0.0f, 1.0f);
  damper_set(p->inputdamper,1.0f - p->inputbandwidth);
}

static inline void gigaverb_set_drylevel(ty_gverb *p, t_floatarg a)
{
  a = CLIP(a, -90.0f, 0.0f);
  p->drylevel = DB_CO(a);
}

static inline void gigaverb_set_wetlevel(ty_gverb *p, t_floatarg a)
{
  a = CLIP(a, -90.f, 0.f);
  p->wetlevel = DB_CO(a);
}

static inline void gigaverb_set_earlylevel(ty_gverb *p, t_floatarg a)
{
  a = CLIP(a, -90.0f, 0.0f);
  p->earlylevel = DB_CO(a);
}

static inline void gigaverb_set_taillevel(ty_gverb *p, t_floatarg a)
{
  a = CLIP(a, -90.0f, 0.0f);
  p->taillevel = DB_CO(a);
}

ty_diffuser *diffuser_make(int size, float coeff)
{
  ty_diffuser *p;
  int i;

  p = (ty_diffuser *)t_getbytes(sizeof(ty_diffuser));
  if(!p) return (NULL);
  p->size = size;
  p->coeff = coeff;
  p->idx = 0;
  p->buf = (float *)t_getbytes(size*sizeof(float));
  if(!p->buf) return (NULL);
  for (i = 0; i < size; i++) p->buf[i] = 0.0;
  return(p);
}

void diffuser_free(ty_diffuser *p)
{
  t_freebytes(p->buf, p->size*sizeof(float));
  t_freebytes(p, sizeof(ty_diffuser));
}

void diffuser_flush(ty_diffuser *p)
{
  memset(p->buf, 0, p->size * sizeof(float));
}

ty_damper *damper_make(float damping)
{
  ty_damper *p;

  p = (ty_damper *)t_getbytes(sizeof(ty_damper));
  if(!p) return (NULL);
  p->damping = damping;
  p->delay = 0.0;
  return(p);
}

void damper_free(ty_damper *p)
{
  t_freebytes(p, sizeof(ty_damper));
}

void damper_flush(ty_damper *p)
{
  p->delay = 0.0f;
}

void fixeddelay_flush(ty_fixeddelay *p)
{
  memset(p->buf, 0, p->size * sizeof(float));
}

ty_fixeddelay *fixeddelay_make(int size)
{
  ty_fixeddelay *p;
  int i;

  p = (ty_fixeddelay *)t_getbytes(sizeof(ty_fixeddelay));
  if(!p) return (NULL);
  p->size = size;
  p->idx = 0;
  p->buf = (float *)t_getbytes(size*sizeof(float));
  if(!p->buf) return (NULL);
  for (i = 0; i < size; i++)
    p->buf[i] = 0.0;
  return(p);
}

void fixeddelay_free(ty_fixeddelay *p)
{
  t_freebytes(p->buf, p->size*sizeof(float));
  t_freebytes(p, sizeof(ty_diffuser));
}

int isprime(int n)
{
  unsigned int i;
  const unsigned int lim = (int)sqrtf((float)n);

  if (n == 2) return(1);
  if ((n & 1) == 0) return(0);
  for(i = 3; i <= lim; i += 2)
    if ((n % i) == 0) return(0);
  return(1);
}

int nearest_prime(int n, float rerror)
/* relative error; new prime will be in range
 * [n-n*rerror, n+n*rerror];
 */
{
  int bound,k;

  if (isprime(n)) return(n);
  /* assume n is large enough and n*rerror enough smaller than n */
  bound = n*rerror;
  for(k = 1; k <= bound; k++) {
    if (isprime(n+k)) return(n+k);
    if (isprime(n-k)) return(n-k);
  }
  return(-1);
}

/* Truncate float to int */
int ff_trunc(float f) {
  f -= 0.5f;
  f += (3<<22);
  return *((int*)&f) - 0x4b400000;
}

/* Round float to int (faster than f_trunc) */
int ff_round(float f) {
  f += (3<<22);
  return *((int*)&f) - 0x4b400000;
}

char *version = "gigaverb~ 1.0test3: \n1999 Juhana Sadeharju \n2004 Olaf Matthes \n2015 Marco Matteo Markidis";

ty_gverb *gigaverb_new(t_symbol *s, short argc, t_atom *argv)
{
  float maxroomsize = 300.0f;
  float roomsize = 50.0f;
  float revtime = 7.0f;
  float damping = 0.5f;
  float spread = 15.0f;
  float inputbandwidth = 0.5f;
  float drylevel = 1.0f; //-1.9832f;
  float wetlevel = 0.f;
  float earlylevel = 0.0f; //-1.9832f;
  float taillevel = 0.0f;

  float ga,gb,gt;
  int i,n;
  float r;
  float diffscale;
  int a,b,c,cc,d,dd,e;
  float spread1,spread2;

  if(argc >= 1)
    maxroomsize = CLIP(argv[0].a_w.w_float, 0.1f, 10000.0f);
  else if(argc >= 2)
    spread = CLIP(argv[1].a_w.w_float, 0.0f, 100.0f);

  ty_gverb *p = (ty_gverb *)pd_new(gigaverb_class);

  outlet_new(&p->obj, gensym("signal"));
  outlet_new(&p->obj, gensym("signal"));

  p->rate = sys_getsr();
  p->fdndamping = damping;
  p->maxroomsize = maxroomsize;
  p->roomsize = CLIP(roomsize, 0.1f, maxroomsize);
  p->revtime = revtime;
  p->drylevel = drylevel;
  p->wetlevel = wetlevel;
  p->earlylevel = earlylevel;
  p->taillevel = taillevel;

  p->maxdelay = p->rate*p->maxroomsize/340.0;
  p->largestdelay = p->rate*p->roomsize/340.0;

  p->bypass = 0;

  if(p->maxroomsize != 300.0f)
    post("gigaverb~: maximum roomsize: %f", p->maxroomsize);

  /* Input damper */

  p->inputbandwidth = inputbandwidth;
  p->inputdamper = damper_make(1.0 - p->inputbandwidth);

  /* FDN section */

  p->fdndels = (ty_fixeddelay **)t_getbytes(FDNORDER*sizeof(ty_fixeddelay *));
  if(!p->fdndels) {
    error("gigaverb~: out of memory");
    return (NULL);
  }
  for(i = 0; i < FDNORDER; i++) {
    p->fdndels[i] = fixeddelay_make((int)p->maxdelay+1000);
    if(!p->fdndels[i]) {
      error("gigaverb~: out of memory");
      return (NULL);
    }
  }
  p->fdngains = (float *)t_getbytes(FDNORDER*sizeof(float));
  p->fdnlens = (int *)t_getbytes(FDNORDER*sizeof(int));
  if(!p->fdngains || !p->fdnlens) {
    error("gigaverb~: out of memory");
    return (NULL);
  }

  p->fdndamps = (ty_damper **)t_getbytes(FDNORDER*sizeof(ty_damper *));
  if(!p->fdndamps) {
    error("gigaverb~: out of memory");
    return (NULL);
  }
  for(i = 0; i < FDNORDER; i++) {
    p->fdndamps[i] = damper_make(p->fdndamping);
    if(!p->fdndamps[i]) {
      error("gigaverb~: out of memory");
      return (NULL);
    }
  }

  ga = 60.0;
  gt = p->revtime;
  ga = pow(10.0,-ga/20.0);
  n = p->rate*gt;
  p->alpha = pow((double)ga,(double)1.0/(double)n);

  gb = 0.0;
  for(i = 0; i < FDNORDER; i++) {
    if (i == 0) gb = 1.000000*p->largestdelay;
    if (i == 1) gb = 0.816490*p->largestdelay;
    if (i == 2) gb = 0.707100*p->largestdelay;
    if (i == 3) gb = 0.632450*p->largestdelay;
    
#if 0
    p->fdnlens[i] = nearest_prime((int)gb, 0.5);
#else
    p->fdnlens[i] = (int)gb;
#endif
      // p->fdngains[i] = -pow(p->alpha,(double)p->fdnlens[i]);
    p->fdngains[i] = -powf((float)p->alpha,p->fdnlens[i]);
  }

  p->d = (float *)t_getbytes(FDNORDER*sizeof(float));
  p->u = (float *)t_getbytes(FDNORDER*sizeof(float));
  p->f = (float *)t_getbytes(FDNORDER*sizeof(float));
  if(!p->d || !p->u || !p->f) {
    error("gigaverb~: out of memory");
    return (NULL);
  }

  /* Diffuser section */

  diffscale = (float)p->fdnlens[3]/(210+159+562+410);
  spread1 = spread;
  spread2 = 3.0*spread;

  b = 210;
  r = 0.125541f;
  a = spread1*r;
  c = 210+159+a;
  cc = c-b;
  r = 0.854046f;
  a = spread2*r;
  d = 210+159+562+a;
  dd = d-c;
  e = 1341-d;

  p->ldifs = (ty_diffuser **)t_getbytes(4*sizeof(ty_diffuser *));
  if(!p->ldifs) {
    error("gigaverb~: out of memory");
    return (NULL);
  }
  p->ldifs[0] = diffuser_make((int)(diffscale*b),0.75);
  p->ldifs[1] = diffuser_make((int)(diffscale*cc),0.75);
  p->ldifs[2] = diffuser_make((int)(diffscale*dd),0.625);
  p->ldifs[3] = diffuser_make((int)(diffscale*e),0.625);
  if(!p->ldifs[0] || !p->ldifs[1] || !p->ldifs[2] || !p->ldifs[3]) {
    error("gigaverb~: out of memory");
    return (NULL);
  }

  b = 210;
  r = -0.568366f;
  a = spread1*r;
  c = 210+159+a;
  cc = c-b;
  r = -0.126815f;
  a = spread2*r;
  d = 210+159+562+a;
  dd = d-c;
  e = 1341-d;

  p->rdifs = (ty_diffuser **)t_getbytes(4*sizeof(ty_diffuser *));
  if(!p->rdifs) {
    error("gigaverb~: out of memory");
    return (NULL);
  }
  p->rdifs[0] = diffuser_make((int)(diffscale*b),0.75);
  p->rdifs[1] = diffuser_make((int)(diffscale*cc),0.75);
  p->rdifs[2] = diffuser_make((int)(diffscale*dd),0.625);
  p->rdifs[3] = diffuser_make((int)(diffscale*e),0.625);
  if(!p->rdifs[0] || !p->rdifs[1] || !p->rdifs[2] || !p->rdifs[3]) {
    error("gigaverb~: out of memory");
    return (NULL);
  }

  /* Tapped delay section */

  p->tapdelay = fixeddelay_make(44000);
  p->taps = (int *)t_getbytes(FDNORDER*sizeof(int));
  p->tapgains = (float *)t_getbytes(FDNORDER*sizeof(float));
  if(!p->tapdelay || !p->taps || !p->tapgains) {
    error("gigaverb~: out of memory");
    return (NULL);
  }

  p->taps[0] = 5+0.410*p->largestdelay;
  p->taps[1] = 5+0.300*p->largestdelay;
  p->taps[2] = 5+0.155*p->largestdelay;
  p->taps[3] = 5+0.000*p->largestdelay;

  for(i = 0; i < FDNORDER; i++) {
    p->tapgains[i] = pow(p->alpha,(double)p->taps[i]);
  }

  s = NULL; /* compiler annoying */
  
  return(p);
}

void gigaverb_free(ty_gverb *p)
{
  int i;

  damper_free(p->inputdamper);
  for (i = 0; i < FDNORDER; i++) {
    fixeddelay_free(p->fdndels[i]);
    damper_free(p->fdndamps[i]);
    diffuser_free(p->ldifs[i]);
    diffuser_free(p->rdifs[i]);
  }
  t_freebytes(p->fdndels, FDNORDER*sizeof(ty_fixeddelay *));
  t_freebytes(p->fdngains, FDNORDER*sizeof(float));
  t_freebytes(p->fdnlens, FDNORDER*sizeof(int));
  t_freebytes(p->fdndamps, FDNORDER*sizeof(ty_damper *));
  t_freebytes(p->d, FDNORDER*sizeof(float));
  t_freebytes(p->u, FDNORDER*sizeof(float));
  t_freebytes(p->f, FDNORDER*sizeof(float));
  t_freebytes(p->ldifs, 4*sizeof(ty_diffuser *));
  t_freebytes(p->rdifs, 4*sizeof(ty_diffuser *));
  t_freebytes(p->taps, FDNORDER*sizeof(int));
  t_freebytes(p->tapgains, FDNORDER*sizeof(float));
  fixeddelay_free(p->tapdelay);
}

t_int *gigaverb_perform(t_int *w)
{
  ty_gverb *p = (ty_gverb *)(w[1]);
  t_float *in = (t_float *)(w[2]);
  t_float *out1 = (t_float *)(w[3]);
  t_float *out2 = (t_float *)(w[4]);
  int n = (int)(w[5]);
  t_float outL, outR, input;
  float dry = p->drylevel;
  float wet = p->wetlevel;
		
  if (p->bypass) {
    /* Bypass, so just copy input to output */
      while(n--) {
	input = *in++;
	*out1++ = input;
	*out2++ = input;
      }
  }
  else {
    /* DSP loop */
    while (n--) {
      input = *in++;
      gigaverb_do(p, input, &outL, &outR);
      *out1++ = outL * wet + input * dry;
      *out2++ = outR * wet + input * dry;
    }
  }

  return (w+6);
}

void gigaverb_dsp(ty_gverb *p, t_signal **sp, short *count) /* count not used */
{
  if (p->rate != sp[0]->s_sr) {
    p->rate = sp[0]->s_sr;
  }
	
  dsp_add(gigaverb_perform, 5, p, sp[0]->s_vec, sp[1]->s_vec,
	  sp[2]->s_vec, sp[0]->s_n);
}

void gigaverb_set_bypass(ty_gverb *p, t_floatarg a)
{
  p->bypass = CLIP(a, 0, 1);
}

void gigaverb_flush(ty_gverb *p)
{
  int i;

  damper_flush(p->inputdamper);
  for(i = 0; i < FDNORDER; i++) {
    fixeddelay_flush(p->fdndels[i]);
    damper_flush(p->fdndamps[i]);
    diffuser_flush(p->ldifs[i]);
    diffuser_flush(p->rdifs[i]);
  }
  memset(p->d, 0, FDNORDER * sizeof(float));
  memset(p->u, 0, FDNORDER * sizeof(float));
  memset(p->f, 0, FDNORDER * sizeof(float));
  fixeddelay_flush(p->tapdelay);
}

/* clear the delay lines and other stuff */
void gigaverb_clear(ty_gverb *p)
{
  int i, k;

  for(i = 0; i < FDNORDER; i++) {
    for(k = 0; k < p->fdnlens[i]; k++) {
      fixeddelay_write(p->fdndels[i], 0);
    }
  }
}

/* print internal values */
void gigaverb_print(ty_gverb *p)
{
  post("gigaverb~ 1.0:");
  post("    roomsize: %0.0f meters (%0.0f maximum)", p->roomsize, p->maxroomsize);
  post("    reverbtime: %0.02f seconds", p->revtime);
  post("    damping: %0.02f", p->fdndamping);
  post("    input bandwidth: %02.02f", p->inputbandwidth);
  post("    dry signal level: %02.02f dB", CO_DB(p->drylevel));
  post("    wet signal level: %02.02f dB", CO_DB(p->wetlevel));
  post("    early reflection level: %02.02f dB", CO_DB(p->earlylevel));
  post("    reverb tail level: %02.02f dB", CO_DB(p->taillevel));
  post("    bypass: %d (Off=0,On=1)", p->bypass);
}

void gigaverb_tilde_setup(void)
{
  t_class *c;
  gigaverb_class = class_new(gensym("gigaverb~"), (t_newmethod)gigaverb_new,
			     (t_method)gigaverb_free, sizeof(ty_gverb),
			     0, A_GIMME, 0);
  CLASS_MAINSIGNALIN(gigaverb_class, ty_gverb, x_f);
  c = gigaverb_class;
  class_addmethod(c,(t_method)gigaverb_dsp, gensym("dsp"), 0);
  class_addmethod(c,(t_method)gigaverb_set_roomsize,
		  gensym("roomsize"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_set_revtime,
		  gensym("revtime"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_set_damping,
		  gensym("damping"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_set_inputbandwidth,
		  gensym("bandwidth"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_set_drylevel, gensym("dry"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_set_wetlevel, gensym("wet"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_set_earlylevel,
		  gensym("early"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_set_taillevel,
		  gensym("tail"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_set_bypass,
		  gensym("bypass"), A_FLOAT, 0);
  class_addmethod(c,(t_method)gigaverb_flush, gensym("clear"), 0);
  class_addmethod(c,(t_method)gigaverb_print, gensym("print"), 0);

  post(version);
}
