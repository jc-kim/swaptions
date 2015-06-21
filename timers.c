#include <stdio.h>
#include <sys/time.h>

#define NUM_TIMERS    64

static double t_start[NUM_TIMERS];
static double t_elapsed[NUM_TIMERS];
static unsigned t_count[NUM_TIMERS];


static inline double get_time()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)1.0e-6*tv.tv_usec;
}


void timer_init() {
  int i;
  for (i = 0; i < NUM_TIMERS; i++) {
    t_start[i] = 0.0;
    t_elapsed[i] = 0.0;
    t_count[i] = 0;
  }
}


void timer_clear(int i) {
  t_elapsed[i] = 0.0;
  t_count[i] = 0;
}


void timer_start(int i) {
  t_start[i] = get_time();
}


void timer_stop(int i) {
  t_elapsed[i] += (get_time() - t_start[i]);
  t_count[i]++;
}


double timer_read(int i) {
  return t_elapsed[i];
}


unsigned timer_count(int i) {
  return t_count[i];
}
