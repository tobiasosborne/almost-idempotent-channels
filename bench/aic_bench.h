/* aic_bench.h — header-only benchmark timing helper (beads aic-73v,
 * "T-harness"). Performance is a first-class gate (CLAUDE.md Law 4: every
 * algorithm candidate is auditioned against correctness AND a performance
 * benchmark; the su2-fft PROFILING.md discipline). This is how candidate routes
 * (matrix-sign vs Newton, projection Route A/B/C, ...) get compared on speed.
 *
 * What it measures. AIC_BENCH(name, n, body) runs `body` n times inside a loop,
 * brackets the loop with clock_gettime(CLOCK_MONOTONIC) — the steady,
 * not-wall-adjusted monotonic clock, the right choice for elapsed-time
 * measurement (CLOCK_REALTIME can jump on NTP steps) — and prints one stable,
 * greppable, machine-readable line:
 *
 *     BENCH <name> n=<N> ns/op=<x> ops/s=<y>
 *
 * Concrete numbers only, no adjectives (CLAUDE.md Rule 12). ns/op is total
 * elapsed nanoseconds / n; ops/s is n / total-elapsed-seconds. The caller is
 * responsible for choosing n large enough that the per-op cost dominates timer
 * granularity (CLOCK_MONOTONIC resolution is ~1-100ns on Linux); a one-line
 * warm-up before the measured loop is the caller's job, kept out of the macro
 * to avoid hidden behaviour.
 *
 * This is deliberately minimal: no statistics framework, no config (CLAUDE.md
 * "no over-build"). A bench file includes this header, sets up its inputs, and
 * calls AIC_BENCH on the hot loop. Compile with _POSIX_C_SOURCE for
 * clock_gettime under -std=c11 (the Makefile passes it).
 */
#ifndef AIC_BENCH_H
#define AIC_BENCH_H

#include <stdio.h>
#include <time.h>

/* Elapsed seconds between two CLOCK_MONOTONIC timespecs (double, sub-ns is
 * below the clock's resolution so double is exact enough for the report). */
static inline double aic_bench_elapsed(const struct timespec *t0,
                                       const struct timespec *t1)
{
    return (double) (t1->tv_sec - t0->tv_sec)
         + (double) (t1->tv_nsec - t0->tv_nsec) * 1e-9;
}

/* Run `body` n times, measure wall-clock with CLOCK_MONOTONIC, print one
 * BENCH line. `name` is a string literal/expr; `n` is the loop count (long);
 * `body` is a statement run each iteration. The loop variable _bi is scoped to
 * the macro. If n <= 0 the report prints ns/op=nan to stay greppable rather
 * than divide by zero. */
#define AIC_BENCH(name, n, body)                                              \
    do {                                                                      \
        const long _bn = (long) (n);                                          \
        struct timespec _bt0, _bt1;                                           \
        clock_gettime(CLOCK_MONOTONIC, &_bt0);                                \
        for (long _bi = 0; _bi < _bn; _bi++) { body; }                        \
        clock_gettime(CLOCK_MONOTONIC, &_bt1);                                \
        double _bsec = aic_bench_elapsed(&_bt0, &_bt1);                       \
        double _bnsop = (_bn > 0) ? (_bsec * 1e9 / (double) _bn) : (0.0/0.0); \
        double _bops = (_bsec > 0) ? ((double) _bn / _bsec) : (0.0/0.0);      \
        printf("BENCH %s n=%ld ns/op=%.3f ops/s=%.3f\n",                      \
               (name), _bn, _bnsop, _bops);                                   \
    } while (0)

#endif /* AIC_BENCH_H */
