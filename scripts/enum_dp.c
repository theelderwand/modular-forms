/* DP-guided eta-quotient enumerator.
   Suffix counts cnt[i][sum+OFF][a1][a2] = # assignments of divisors[i:]
   with that sum and Newman residues. DFS follows only count>0 branches.
   Usage: enum96b N weight max_abs
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXD 16
#define OFF 64
#define SZ 129

static int N, W, MAXA, ND;
static int divs[MAXD];
static int r[MAXD];
static long long count_out = 0;
/* cnt[i][s][a][b]: i in 0..ND, s in 0..2*OFF (sum+OFF), a,b in 0..23.
   Stored as a flat array indexed [i][s][a][b], capped at 2 (we only need >0,
   but counts help nothing else; cap at 2 to avoid overflow). */
static unsigned char cnt[MAXD + 1][SZ][24][24];

static int igcd(int a, int b) { return b ? igcd(b, a % b) : a; }

static int newman_ok(void) {
    long s1 = 0, s2 = 0;
    for (int i = 0; i < ND; i++) {
        s1 += (long)r[i] * divs[i];
        s2 += (long)r[i] * (N / divs[i]);
    }
    return (s1 % 24 == 0) && (s2 % 24 == 0);
}

static int ligozat_ok(void) {
    for (int ci = 0; ci < ND; ci++) {
        int c = divs[ci];
        int g = igcd(c * c, N);
        long num = 0;
        for (int i = 0; i < ND; i++) {
            int gd = igcd(divs[i], c);
            num += (long)r[i] * gd * gd * (N / divs[i]);
        }
        if (num < 0) return 0;
    }
    return 1;
}

static void emit(void) {
    int first = 1;
    for (int i = 0; i < ND; i++) {
        if (r[i] != 0) {
            printf("%s%d:%d", first ? "" : ",", divs[i], r[i]);
            first = 0;
        }
    }
    printf("\n");
    count_out++;
}

static void rec(int i, int rem, int a1, int a2) {
    if (i == ND) {
        if (rem == 0 && a1 == 0 && a2 == 0 && newman_ok() && ligozat_ok()) emit();
        return;
    }
    int d = divs[i];
    int nd = N / d;
    for (int v = -MAXA; v <= MAXA; v++) {
        int rem2 = rem - v;
        int b1 = (a1 + v * d) % 24; if (b1 < 0) b1 += 24;
        int b2 = (a2 + v * nd) % 24; if (b2 < 0) b2 += 24;
        int need1 = (24 - b1) % 24, need2 = (24 - b2) % 24;
        int sidx = rem2 + OFF;
        if (sidx < 0 || sidx >= SZ) continue;
        if (!cnt[i + 1][sidx][need1][need2]) continue;
        r[i] = v;
        rec(i + 1, rem2, b1, b2);
    }
    r[i] = 0;
}

int main(int argc, char **argv) {
    N = atoi(argv[1]); W = atoi(argv[2]); MAXA = atoi(argv[3]);
    ND = 0;
    for (int d = 1; d <= N; d++) if (N % d == 0) divs[ND++] = d;
    memset(r, 0, sizeof r);
    memset(cnt, 0, sizeof cnt);
    cnt[ND][OFF][0][0] = 1;
    for (int i = ND - 1; i >= 0; i--) {
        int d = divs[i], nd = N / divs[i];
        for (int s = 0; s < SZ; s++)
            for (int a = 0; a < 24; a++)
                for (int b = 0; b < 24; b++) {
                    if (!cnt[i + 1][s][a][b]) continue;
                    for (int v = -MAXA; v <= MAXA; v++) {
                        int s2 = s + v;
                        if (s2 < 0 || s2 >= SZ) continue;
                        int a2 = (a + v * d) % 24; if (a2 < 0) a2 += 24;
                        int b2 = (b + v * nd) % 24; if (b2 < 0) b2 += 24;
                        if (cnt[i][s2][a2][b2] < 2) cnt[i][s2][a2][b2]++;
                    }
                }
    }
    rec(0, 2 * W, 0, 0);
    fprintf(stderr, "# %lld edicts\n", count_out);
    return 0;
}
