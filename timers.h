#ifndef __TIMERS_H__
#define __TIMERS_H__

/* Timer functions */
void timer_init();
void timer_clear(int i);
void timer_start(int i);
void timer_stop(int i);
double   timer_read(int i);
unsigned timer_count(int i);

#endif //__TIMERS_H__
