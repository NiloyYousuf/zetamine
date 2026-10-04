# =====================================================================================
# ZetaMine vs. FP-Growth / Eclat / LCM  --  Windows-native benchmark  (v2)
# Needs: g++ on PATH, Python 3, Java (for SPMF). Optional: pyfim package.
# Outputs land in .\zm_bench\ (relative to the current directory).
#
# What changed vs v1
#   1. Timing is split into parse / prep / mine for ZetaMine; headline time `t` now
#      INCLUDES file loading and preprocessing, and a second metric `wall` (whole
#      subprocess) is reported for every method.
#   2. Two like-for-like comparisons: OUTPUT mode (every method emits all itemsets) and
#      COUNT-only mode (ZetaMine-count vs PyFIM report='#'). No more "ZM counts, others list".
#   3. PyFIM no longer runs the minsup-1 workaround or the overhead subtraction; optional
#      native baselines (Borgelt CLI, Uno's LCM) can be plugged in via EXTERNAL_BASELINES.
#   4. The "PyFIM boundary quirk" was a misdiagnosis. Borgelt's fim omits every itemset
#      made up ONLY of items that occur in all transactions (2^c - 1 itemsets for c such
#      items). Counts are corrected by that amount and verified itemset by itemset.
# =====================================================================================
import os, sys, re, json, math, time, gzip, shutil, platform, subprocess, urllib.request
import importlib, collections, itertools
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IS_WINDOWS = platform.system() == "Windows"
NULL_DEV   = "NUL" if IS_WINDOWS else "/dev/null"
EXE_SUFFIX = ".exe" if IS_WINDOWS else ""

# ------------------------------- USER SETTINGS ---------------------------------------
PRESET        = os.environ.get("ZM_PRESET", "full")
REPS          = int(os.environ.get("ZM_REPS", "3"))
TIMEOUT       = int(os.environ.get("ZM_TIMEOUT", "600"))
SLOW_S        = 30
USE_SPMF      = True
RUN_ABLATIONS = os.environ.get("ZM_ABLATIONS", "1") == "1"
UNIT_TESTS    = int(os.environ.get("ZM_UNIT_TESTS", "40"))
REAL_CHECK_MAX = int(os.environ.get("ZM_REAL_CHECK_MAX", "1500000"))  # skip itemset-level check above this many itemsets (RAM)
LIST_MAX       = int(os.environ.get("ZM_LIST_MAX", "20000000"))      # skip PyFIM list mode (RAM) above this many itemsets
LIST_MAX      = int(os.environ.get("ZM_LIST_MAX", "20000000"))  # skip PyFIM list mode above this many itemsets (RAM: ~150-200 B each)
BRUTE_TESTS   = 15          # ZetaMine vs exhaustive enumeration (independent of every library)
ZM_MAIN = dict(kmax=20, ksmall=8, theta=0.2)
ABLATIONS = {"ZM-noKernel": dict(kmax=0,  ksmall=0, theta=0.0),
             "ZM-k16":      dict(kmax=16, ksmall=8, theta=0.2),
             "ZM-k22":      dict(kmax=22, ksmall=8, theta=0.2)}
PRESETS = {
 "smoke": {"mushroom": [0.30], "chess": [0.70], "retail": [0.01]},
 "mid":   {"mushroom": [0.10], "chess": [0.50], "retail": [0.001], "T10I4D100K": [0.001]},
 "quick": {"mushroom": [0.20, 0.10], "chess": [0.50, 0.35], "connect": [0.90, 0.80],
           "retail": [0.005, 0.001], "T10I4D100K": [0.005, 0.001]},
 "full":  {"mushroom": [0.30, 0.20, 0.10, 0.05], "chess": [0.60, 0.50, 0.40, 0.30],
           "connect": [0.95, 0.90, 0.85, 0.80], "pumsb": [0.90, 0.85, 0.80, 0.75],
           "accidents": [0.60, 0.50, 0.40, 0.30], "retail": [0.01, 0.005, 0.002, 0.001],
           "T10I4D100K": [0.01, 0.005, 0.002, 0.001], "T40I10D100K": [0.05, 0.03, 0.02, 0.01],
           "kosarak": [0.01, 0.005, 0.003, 0.002]},
}

# Optional NATIVE baselines (compile them yourself; not run unless filled in).
#   group="count"  -> compared against ZetaMine-count     group="output" -> against ZetaMine
#   {data} / {minsup} are substituted. count_re extracts the itemset count from stdout+stderr.
#   borgelt=True applies the same "full-support items" correction as PyFIM (same C code base).
# Command lines are UNTESTED templates -- check them against each tool's usage text.
EXTERNAL_BASELINES = {
    # "Borgelt-FPGrowth": dict(cmd=[r"C:\tools\fpgrowth.exe", "-ts", "-s-{minsup}", "{data}"],
    #                          group="count", count_re=r"\[(\d+) set\(s\)\]", borgelt=True),
    # "Borgelt-Eclat":    dict(cmd=[r"C:\tools\eclat.exe", "-ts", "-s-{minsup}", "{data}"],
    #                          group="count", count_re=r"\[(\d+) set\(s\)\]", borgelt=True),
    # "Uno-LCM":          dict(cmd=[r"C:\tools\lcm.exe", "F", "{data}", "{minsup}", NULL_DEV],
    #                          group="output", count_re=None, borgelt=False),
}

VERIFIED_MIRRORS = {
    "chess":       ["https://raw.githubusercontent.com/mrinmoy29a/Apriori-Algorithm/master/chess.dat"],
    "mushroom":    ["https://raw.githubusercontent.com/PoChuanHuang/FPgrowth_Algorithm/master/DataSet/mushrooms.txt"],
    "retail":      ["https://raw.githubusercontent.com/xikuun/Data_Mining/master/retail.dat"],
    "T10I4D100K":  ["https://raw.githubusercontent.com/Sujit26/DataMining-Algorithms/main/T10I4D100K.dat.txt"],
    "T40I10D100K": ["https://raw.githubusercontent.com/Sujit26/DataMining-Algorithms/main/T40I10D100K.dat.txt"],
}
USER_MIRRORS = {}

# --------------------------------------------------------------------------------------
WORK = os.path.abspath("zm_bench")
DATA = os.path.join(WORK, "data")
os.makedirs(DATA, exist_ok=True)
os.chdir(WORK)

# ------------------------------- progress helpers ------------------------------------
_T0 = time.time()
def log(msg):
    print(f"[{time.time() - _T0:7.1f}s] {msg}", flush=True)

USER_AGENT = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
              "(KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36")

def sh(cmd, timeout=None):
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                          encoding="utf-8", errors="replace")

def timed_sh(cmd, timeout=None):
    t0 = time.perf_counter()
    p = sh(cmd, timeout)
    return p, time.perf_counter() - t0

def fetch(url, timeout=120):
    req = urllib.request.Request(url, headers={
        "User-Agent": USER_AGENT, "Accept": "*/*", "Accept-Encoding": "gzip, identity"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()
        enc = (resp.headers.get("Content-Encoding") or "").lower()
        if "gzip" in enc or url.endswith(".gz"):
            try: raw = gzip.decompress(raw)
            except OSError: pass
        return raw

# ---------------------------- 1. ZetaMine C++ source ----------------------------------
# Reports parse_ms (file -> vectors), prep_ms (rank/sort/dedup/allocate), mine_ms (DFS+kernel+output).
CPP = r'''
#include <bits/stdc++.h>
using namespace std;
typedef uint32_t u32; typedef uint64_t u64;
static const u32 NONE = 0xFFFFFFFFu;
struct Rows { vector<u32> it, off{0}, w; size_t n() const { return w.size(); } };
static u32 MINSUP; static int KMAX, KSMALL; static double THETA;
static u64 TOTAL = 0, KCALLS = 0, KSETS = 0, DNODES = 0;
static vector<u32> cntA, locIdx, hist, rankToItem, prefix_;
static FILE* OUT = nullptr; static vector<char> obuf;
static inline void flush_obuf() { if (OUT && !obuf.empty()) { fwrite(obuf.data(), 1, obuf.size(), OUT); obuf.clear(); } }
static inline void put_u32(u32 v) { char t[12]; int n = 0; do { t[n++] = '0' + v % 10; v /= 10; } while (v); while (n) obuf.push_back(t[--n]); }
static void emit_set(const vector<u32>& extra, u32 sup) {
  for (u32 r : prefix_) { put_u32(rankToItem[r]); obuf.push_back(' '); }
  for (u32 r : extra)   { put_u32(rankToItem[r]); obuf.push_back(' '); }
  obuf.push_back('('); put_u32(sup); obuf.push_back(')'); obuf.push_back('\n');
  if (obuf.size() > (1u << 20)) flush_obuf();
}
static void zeta(u32* h, int k) {
  const int B = min(k, 13); const u32 N = 1u << k, blk = 1u << B;
  for (u32 base = 0; base < N; base += blk)
    for (int i = 0; i < B; i++) { u32 s = 1u << i; u32* p = h + base;
      for (u32 j0 = 0; j0 < blk; j0 += 2 * s) for (u32 j = j0; j < j0 + s; j++) p[j] += p[j + s]; }
  for (int i = B; i < k; i++) { u32 s = 1u << i;
    for (u32 j0 = 0; j0 < N; j0 += 2 * s) for (u32 j = j0; j < j0 + s; j++) h[j] += h[j + s]; }
}
static void rec(const Rows& db) {
  DNODES++;
  vector<u32> touched; u64 W = 0; const size_t R = db.n();
  for (size_t r = 0; r < R; r++) { u32 w = db.w[r]; W += w;
    for (u32 q = db.off[r]; q < db.off[r + 1]; q++) { u32 e = db.it[q]; if (cntA[e] == 0) touched.push_back(e); cntA[e] += w; } }
  vector<pair<u32,u32>> fl; u64 sumsup = 0;
  for (u32 e : touched) { if (cntA[e] >= MINSUP) { fl.push_back({e, cntA[e]}); sumsup += cntA[e]; } cntA[e] = 0; }
  const int k = (int)fl.size(); if (k == 0) return;
  sort(fl.begin(), fl.end());
  const bool use_kernel = (k <= KSMALL) || (k <= KMAX && (double)sumsup >= THETA * (double)k * (double)W);
  if (use_kernel) {
    KCALLS++;
    for (int j = 0; j < k; j++) locIdx[fl[j].first] = j;
    const u32 N = 1u << k; memset(hist.data(), 0, sizeof(u32) * (size_t)N);
    for (size_t r = 0; r < R; r++) { u32 mask = 0;
      for (u32 q = db.off[r]; q < db.off[r + 1]; q++) { u32 l = locIdx[db.it[q]]; if (l != NONE) mask |= 1u << l; }
      hist[mask] += db.w[r]; }
    for (int j = 0; j < k; j++) locIdx[fl[j].first] = NONE;
    zeta(hist.data(), k); hist[0] = 0;
    if (!OUT) { u64 c = 0; const u32 ms = MINSUP; const u32* h = hist.data();
      for (u32 m = 1; m < N; m++) c += (h[m] >= ms); TOTAL += c; KSETS += c; }
    else { vector<u32> ex;
      for (u32 m = 1; m < N; m++) if (hist[m] >= MINSUP) { ex.clear();
        for (int b = 0; b < k; b++) if ((m >> b) & 1) ex.push_back(fl[b].first);
        emit_set(ex, hist[m]); TOTAL++; KSETS++; } }
    return;
  }
  for (int j = 0; j < k; j++) locIdx[fl[j].first] = j;
  Rows D; D.it.reserve(db.it.size()); D.off.reserve(R + 1); D.w.reserve(R);
  for (size_t r = 0; r < R; r++) { size_t before = D.it.size();
    for (u32 q = db.off[r]; q < db.off[r + 1]; q++) if (locIdx[db.it[q]] != NONE) D.it.push_back(db.it[q]);
    if (D.it.size() > before) { D.off.push_back((u32)D.it.size()); D.w.push_back(db.w[r]); } }
  const size_t R2 = D.n();
  vector<u32> bstart(k + 1, 0);
  for (size_t r = 0; r < R2; r++) for (u32 q = D.off[r]; q < D.off[r + 1]; q++) bstart[locIdx[D.it[q]] + 1]++;
  for (int j = 0; j < k; j++) bstart[j + 1] += bstart[j];
  vector<u32> eRow(bstart[k]), ePos(bstart[k]), fillp(bstart.begin(), bstart.end() - 1);
  for (size_t r = 0; r < R2; r++) for (u32 q = D.off[r]; q < D.off[r + 1]; q++) {
    u32 j = locIdx[D.it[q]]; u32 t = fillp[j]++; eRow[t] = (u32)r; ePos[t] = q; }
  for (int j = 0; j < k; j++) locIdx[fl[j].first] = NONE;
  for (int j = 0; j < k; j++) {
    Rows ch; ch.w.reserve(bstart[j + 1] - bstart[j]);
    for (u32 t = bstart[j]; t < bstart[j + 1]; t++) { u32 r = eRow[t], q = ePos[t];
      for (u32 z = q + 1; z < D.off[r + 1]; z++) ch.it.push_back(D.it[z]);
      if (ch.it.size() > ch.off.back()) { ch.off.push_back((u32)ch.it.size()); ch.w.push_back(D.w[r]); } }
    prefix_.push_back(fl[j].first);
    if (OUT) { vector<u32> none; emit_set(none, fl[j].second); }
    TOTAL++;
    if (ch.n()) rec(ch);
    prefix_.pop_back();
  }
}
int main(int argc, char** argv) {
  if (argc < 7) { fprintf(stderr, "usage: zetamine file minsup kmax ksmall theta out\n"); return 1; }
  MINSUP = (u32)atoll(argv[2]); KMAX = atoi(argv[3]); KSMALL = atoi(argv[4]); THETA = atof(argv[5]);
  string outp = argv[6]; if (KMAX > 26) KMAX = 26; if (KMAX < 0) KMAX = 0; if (KSMALL > KMAX) KSMALL = KMAX;
  auto T0 = chrono::steady_clock::now();
  vector<vector<u32>> raw;
  { ifstream f(argv[1]); string line;
    while (getline(f, line)) { vector<u32> t; const char* s = line.c_str(); char* e;
      while (*s) { while (*s == ' ' || *s == '\t' || *s == '\r') s++; if (!*s) break;
        unsigned long v = strtoul(s, &e, 10); if (e == s) break; t.push_back((u32)v); s = e; }
      if (!t.empty()) { sort(t.begin(), t.end()); t.erase(unique(t.begin(), t.end()), t.end()); raw.push_back(move(t)); } } }
  auto T1 = chrono::steady_clock::now();
  u32 maxid = 0; for (auto& t : raw) maxid = max(maxid, t.back());
  vector<u32> g(maxid + 1, 0); for (auto& t : raw) for (u32 e : t) g[e]++;
  vector<u32> items; for (u32 e = 0; e <= maxid; e++) if (g[e] >= MINSUP) items.push_back(e);
  sort(items.begin(), items.end(), [&](u32 a, u32 b) { return g[a] != g[b] ? g[a] < g[b] : a < b; });
  const u32 m = (u32)items.size(); rankToItem = items;
  vector<u32> rk(maxid + 1, NONE); for (u32 r = 0; r < m; r++) rk[items[r]] = r;
  vector<vector<u32>> rows;
  for (auto& t : raw) { vector<u32> v; for (u32 e : t) if (rk[e] != NONE) v.push_back(rk[e]);
    if (!v.empty()) { sort(v.begin(), v.end()); rows.push_back(move(v)); } }
  vector<vector<u32>>().swap(raw);
  sort(rows.begin(), rows.end());
  Rows root;
  for (size_t i = 0; i < rows.size();) { size_t j = i; while (j < rows.size() && rows[j] == rows[i]) j++;
    for (u32 e : rows[i]) root.it.push_back(e);
    root.off.push_back((u32)root.it.size()); root.w.push_back((u32)(j - i)); i = j; }
  vector<vector<u32>>().swap(rows);
  cntA.assign(m + 1, 0); locIdx.assign(m + 1, NONE); hist.assign((size_t)1 << KMAX, 0);
  if (outp != "none") OUT = fopen(outp.c_str(), "w");
  auto T1b = chrono::steady_clock::now();
  rec(root);
  flush_obuf(); if (OUT) fclose(OUT);
  auto T2 = chrono::steady_clock::now();
  printf("RESULT count=%llu parse_ms=%.3f prep_ms=%.3f mine_ms=%.3f kernel_calls=%llu kernel_sets=%llu dfs_nodes=%llu\n",
         (unsigned long long)TOTAL,
         chrono::duration<double, milli>(T1 - T0).count(),
         chrono::duration<double, milli>(T1b - T1).count(),
         chrono::duration<double, milli>(T2 - T1b).count(),
         (unsigned long long)KCALLS, (unsigned long long)KSETS, (unsigned long long)DNODES);
  return 0;
}
'''

# ---------------------------- 2. PyFIM worker -----------------------------------------
# mode "list"  : report='a' -> materialises every itemset as a Python tuple (its native output)
# mode "count" : report='#' -> C-level histogram of (size, support) -> count only, no list
# Timed region = the fim call only (the DB is already a Python list). No overhead subtraction,
# no minsup-1 workaround: supp=-minsup is exact. nfull = #items present in ALL transactions;
# Borgelt's code omits the 2^nfull - 1 itemsets made only of such items (corrected by the caller).
WORKER = r'''
import sys, time, json, gc, fim
path, minsup, algo, mode = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
tr = []
with open(path) as f:
    for line in f:
        t = line.split()
        if t: tr.append(sorted(set(int(x) for x in t)))
n = len(tr); cnt = {}
for t in tr:
    for i in t: cnt[i] = cnt.get(i, 0) + 1
nfull = sum(1 for c in cnt.values() if c == n)
fn = {"fpgrowth": fim.fpgrowth, "eclat": fim.eclat}[algo]
gc.collect(); gc.disable()
t0 = time.perf_counter()
if mode == "count":
    r = fn(tr, target="s", supp=-minsup, zmin=1, report="#")
    dt = time.perf_counter() - t0
    c = int(round(sum(r.values())))
else:
    r = fn(tr, target="s", supp=-minsup, zmin=1, report="a")
    dt = time.perf_counter() - t0
    c = len(r)
print(json.dumps({"count": c, "t": dt, "n": n, "nfull": nfull}))
'''

# ---------------------------- 3. Setup -------------------------------------------------
log("== Setup ==")
log(f"platform: {platform.system()} {platform.release()}  python: {sys.version.split()[0]}  preset: {PRESET}")

gpp_bin = shutil.which("g++")
if gpp_bin is None:
    raise SystemExit("g++ not found on PATH. Install MSYS2 mingw-w64-x86_64-gcc and add C:\\msys64\\mingw64\\bin to PATH.")
log(f"g++: {gpp_bin}")

log("writing zetamine.cpp and worker_pyfim.py")
open("zetamine.cpp", "w").write(CPP)
open("worker_pyfim.py", "w").write(WORKER)

log("compiling zetamine.cpp ...")
ZM_BIN = os.path.join(WORK, "zm" + EXE_SUFFIX)
gpp = sh([gpp_bin, "-O3", "-std=c++17", "-DNDEBUG", "-o", ZM_BIN, "zetamine.cpp"])
if gpp.returncode != 0:
    print(gpp.stderr[:3000]); raise SystemExit("C++ compile failed")
log(f"compiled zetamine -> {ZM_BIN}")

HAVE_FIM = False
try:
    log("checking PyFIM ...")
    sh([sys.executable, "-m", "pip", "install", "-q", "pyfim"])
    importlib.invalidate_caches()
    import fim
    HAVE_FIM = True
    log("PyFIM OK")
except Exception as e:
    log(f"!! PyFIM unavailable ({e}). Will use SPMF for FP-Growth/Eclat baselines.")

HAVE_SPMF = False
SPMF_JAR = os.path.join(WORK, "spmf.jar")
SPMF_URLS = [
    "https://raw.githubusercontent.com/loopGod/maximal-sequential-patterns-mining/master/spmf.jar",
    "https://raw.githubusercontent.com/Sandy4321/SequencedbSPMF/main/spmf.jar",
    "https://www.philippe-fournier-viger.com/spmf/spmf.jar",
]
if USE_SPMF and shutil.which("java"):
    if not os.path.exists(SPMF_JAR):
        log("downloading SPMF jar ...")
        for url in SPMF_URLS:
            try:
                blob = fetch(url, timeout=180)
                if len(blob) > 100_000 and blob[:2] == b"PK":
                    open(SPMF_JAR, "wb").write(blob)
                    log(f"  SPMF jar from {url} ({len(blob)/1024/1024:.1f} MB)")
                    break
            except Exception as e:
                log(f"  mirror failed: {url} ({e})")
    HAVE_SPMF = os.path.exists(SPMF_JAR) and os.path.getsize(SPMF_JAR) > 100_000
log(f"SPMF baselines: {'on' if HAVE_SPMF else 'off'}")
for k, v in EXTERNAL_BASELINES.items():
    log(f"external baseline enabled: {k} ({v['group']})")

# ---------------------------- 4. Runners ----------------------------------------------
# Every runner returns: status, count, t (headline), wall (whole subprocess, incl. startup).
#   ZetaMine : t = parse + prep + mine  (file loading INCLUDED -> conservative for ZetaMine)
#   SPMF     : t = SPMF's own "Total time" (also includes reading the file + JIT warm-up)
#   PyFIM    : t = fim call only on an already-loaded Python list (loading EXCLUDED -> favours PyFIM)
def run_zm(path, minsup, cfg, outfile="none", timeout=TIMEOUT):
    cmd = [ZM_BIN, path, str(minsup), str(cfg["kmax"]), str(cfg["ksmall"]), str(cfg["theta"]), outfile]
    try: p, wall = timed_sh(cmd, timeout)
    except subprocess.TimeoutExpired: return dict(status="timeout")
    m = re.search(r"RESULT count=(\d+) parse_ms=([\d.]+) prep_ms=([\d.]+) mine_ms=([\d.]+) "
                  r"kernel_calls=(\d+) kernel_sets=(\d+) dfs_nodes=(\d+)", p.stdout)
    if p.returncode != 0 or not m: return dict(status="error", msg=(p.stderr or p.stdout)[-300:])
    tp, tq, tm = float(m[2]) / 1000, float(m[3]) / 1000, float(m[4]) / 1000
    return dict(status="ok", count=int(m[1]), t=tp + tq + tm, wall=wall,
                t_parse=tp, t_prep=tq, t_mine=tm,
                kernel_calls=int(m[5]), kernel_sets=int(m[6]), dfs_nodes=int(m[7]))

def run_pyfim(path, minsup, algo, mode, timeout=TIMEOUT):
    try: p, wall = timed_sh([sys.executable, "worker_pyfim.py", path, str(minsup), algo, mode], timeout)
    except subprocess.TimeoutExpired: return dict(status="timeout")
    try: d = json.loads(p.stdout.strip().splitlines()[-1])
    except Exception: return dict(status="error", msg=(p.stderr or p.stdout)[-300:])
    missing = (2 ** d["nfull"] - 1) if d["nfull"] <= 30 else 0   # itemsets Borgelt's fim omits
    return dict(status="ok", count=d["count"] + missing, count_raw=d["count"],
                t=max(d["t"], 1e-6), wall=wall)

def run_spmf(path, minsup, n, algo, timeout=TIMEOUT):
    pct = (minsup - 0.5) / n * 100.0
    try:
        p, wall = timed_sh(["java", "-Xmx4g", "-jar", SPMF_JAR, "run", algo, path, NULL_DEV, f"{pct:.8f}%"], timeout)
    except subprocess.TimeoutExpired: return dict(status="timeout")
    t = re.search(r"(?i)total time\s*~?\s*:?\s*(\d+)\s*ms", p.stdout)
    c = re.search(r"(?i)itemsets?\s+count\s*:?\s*(\d+)", p.stdout)
    if not (t and c): return dict(status="error", msg=p.stdout[-300:])
    return dict(status="ok", count=int(c[1]), t=int(t[1]) / 1000.0, wall=wall)

def run_external(spec, path, minsup, nfull, timeout=TIMEOUT):
    cmd = [a.format(data=path, minsup=minsup) for a in spec["cmd"]]
    try: p, wall = timed_sh(cmd, timeout)
    except subprocess.TimeoutExpired: return dict(status="timeout")
    if p.returncode != 0: return dict(status="error", msg=(p.stderr or p.stdout)[-300:])
    cnt = None
    if spec.get("count_re"):
        m = re.search(spec["count_re"], (p.stdout or "") + (p.stderr or ""))
        if m:
            cnt = int(m[1])
            if spec.get("borgelt") and nfull <= 30: cnt += 2 ** nfull - 1
    return dict(status="ok", count=cnt, t=wall, wall=wall)   # native tools: only wall-clock available

# group names: which baselines are like-for-like with which ZetaMine mode
ZM_OUT, ZM_CNT = "ZetaMine", "ZetaMine-count"
def make_methods(n, nfull):
    M = [(ZM_OUT, lambda p, ms: run_zm(p, ms, ZM_MAIN, outfile=NULL_DEV)),   # emits every itemset (text, to NUL)
         (ZM_CNT, lambda p, ms: run_zm(p, ms, ZM_MAIN))]                    # count only
    if RUN_ABLATIONS:                                                         # count-only, compare vs ZetaMine-count
        for nm_, cfg in ABLATIONS.items(): M.append((nm_, lambda p, ms, cfg=cfg: run_zm(p, ms, cfg)))
    if HAVE_FIM:
        M += [("PyFIM-FPGrowth",       lambda p, ms: run_pyfim(p, ms, "fpgrowth", "list")),
              ("PyFIM-Eclat",          lambda p, ms: run_pyfim(p, ms, "eclat", "list")),
              ("PyFIM-FPGrowth-count", lambda p, ms: run_pyfim(p, ms, "fpgrowth", "count")),
              ("PyFIM-Eclat-count",    lambda p, ms: run_pyfim(p, ms, "eclat", "count"))]
    if HAVE_SPMF:
        M += [("SPMF-LCM",      lambda p, ms: run_spmf(p, ms, n, "LCMFreq")),
              ("SPMF-FPGrowth", lambda p, ms: run_spmf(p, ms, n, "FPGrowth_itemsets"))]
        if not HAVE_FIM:
            M += [("SPMF-Eclat", lambda p, ms: run_spmf(p, ms, n, "Eclat"))]
    for nm_, spec in EXTERNAL_BASELINES.items():
        M.append((nm_, lambda p, ms, spec=spec: run_external(spec, p, ms, nfull)))
    return M

OUT_BASE = ["PyFIM-FPGrowth", "PyFIM-Eclat", "SPMF-LCM", "SPMF-FPGrowth", "SPMF-Eclat"] + \
           [k for k, v in EXTERNAL_BASELINES.items() if v["group"] == "output"]
CNT_BASE = ["PyFIM-FPGrowth-count", "PyFIM-Eclat-count"] + \
           [k for k, v in EXTERNAL_BASELINES.items() if v["group"] == "count"]
ALL_BASE = OUT_BASE + CNT_BASE

# ---------------------------- 5. Correctness tests ------------------------------------
def read_zm_out(fn):
    d, nl = {}, 0
    with open(fn) as f:
        for line in f:
            a, b = line.rsplit(" (", 1); nl += 1
            d[frozenset(int(x) for x in a.split())] = int(b.strip().rstrip(")"))
    return d, nl

def pyfim_raw(tr, minsup):
    """PyFIM exactly as shipped: exact absolute threshold, no tricks."""
    raw = fim.fpgrowth(tr, target="s", supp=-minsup, zmin=1, report="a")
    return {frozenset(i): int(round(s)) for i, s in raw}

def full_support_closure(tr, minsup):
    """All non-empty subsets of items that occur in EVERY transaction (support = n).
    Borgelt's fim does not report these; they are frequent by definition."""
    n = len(tr); cnt = collections.Counter(i for t in tr for i in t)
    full = [i for i, c in cnt.items() if c == n]
    if n < minsup or not full or len(full) > 20: return {}
    return {frozenset(c): n for r in range(1, len(full) + 1) for c in itertools.combinations(full, r)}

def brute_force(tr, minsup):
    """Exhaustive enumeration, independent of every library. Tiny inputs only."""
    items = sorted({i for t in tr for i in t}); sets = [frozenset(t) for t in tr]; out = {}
    level = [(frozenset([i]),) for i in items]
    for k in range(1, len(items) + 1):
        found = False
        for c in itertools.combinations(items, k):
            s = frozenset(c); sup = sum(1 for t in sets if s <= t)
            if sup >= minsup: out[s] = sup; found = True
        if not found: break
    return out

def load_tr(path):
    tr = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if t: tr.append(sorted(set(int(x) for x in t)))
    return tr

TEST_CFGS = [dict(kmax=0, ksmall=0, theta=0.0), dict(kmax=3, ksmall=2, theta=0.0),
             dict(kmax=12, ksmall=6, theta=0.0), ZM_MAIN]

def zm_matches(path, minsup, ref, tag):
    """Strict equality (itemsets AND supports AND line count AND count-only count) for all configs."""
    ok = True
    for cfg in TEST_CFGS:
        outf = os.path.join(WORK, "tmp_out.txt")
        r = run_zm(path, minsup, cfg, outfile=outf)
        rc = run_zm(path, minsup, cfg)             # count-only path must agree too
        if r["status"] != "ok" or rc["status"] != "ok":
            print(f"     {tag} cfg={cfg} -> {r['status']}/{rc['status']}"); ok = False; continue
        got, nl = read_zm_out(outf)
        if got == ref and nl == len(ref) and r["count"] == len(ref) and rc["count"] == len(ref): continue
        ok = False
        print(f"     MISMATCH {tag} cfg={cfg}: zm={len(got)} lines={nl} count={r['count']}/{rc['count']} ref={len(ref)}")
        ex_z, ex_r = set(got) - set(ref), set(ref) - set(got)
        if ex_z: print(f"       ZM only ({len(ex_z)}): {sorted(map(sorted, ex_z))[:3]}")
        if ex_r: print(f"       ref only ({len(ex_r)}): {sorted(map(sorted, ex_r))[:3]}")
    return ok

def check_vs_pyfim(path, tr, minsup, tag):
    """Reference = PyFIM output + the documented omission (full-support closure). Strict equality.
    Anything else that differs is a real failure."""
    raw = pyfim_raw(tr, minsup); exp = full_support_closure(tr, minsup)
    added = [k for k in exp if k not in raw]
    ref = dict(raw); ref.update({k: exp[k] for k in added})
    if added:
        print(f"     NOTE {tag}: PyFIM omitted {len(added)} itemset(s) made only of items present in all "
              f"{len(tr)} transactions, e.g. {sorted(map(sorted, added))[:3]} (support = n). Added to the reference.")
    return zm_matches(path, minsup, ref, tag), len(ref), len(added)

CORRECT = {"brute_tests_passed": None, "brute_tests_total": 0,
           "unit_tests_passed": None, "unit_tests_total": 0,
           "pyfim_omitted_itemsets_total": 0, "real_data_verified": {}}

log(f"== Correctness (1/2): {BRUTE_TESTS} tiny DBs vs exhaustive enumeration (no library involved) ==")
rng = np.random.default_rng(777); bp = 0
for t in range(BRUTE_TESTS):
    n = int(rng.integers(8, 60)); m = int(rng.integers(3, 9))
    X = rng.random((n, m)) < rng.uniform(0.2, 0.95, size=m)
    if t % 3 == 0: X[:, 0] = True                      # exercise the "item in every transaction" case
    path = os.path.join(WORK, "tmp_brute.dat")
    with open(path, "w") as f:
        for r in X:
            idx = np.flatnonzero(r)
            if len(idx): f.write(" ".join(str(int(i) + 1) for i in idx) + "\n")
    tr = load_tr(path)
    if not tr: continue
    ms = int(rng.integers(1, max(2, len(tr) // 2)))
    ok = zm_matches(path, ms, brute_force(tr, ms), f"brute{t}")
    CORRECT["brute_tests_total"] += 1; bp += int(ok)
CORRECT["brute_tests_passed"] = bp
log(f"brute-force tests: {bp}/{CORRECT['brute_tests_total']}")
if bp != CORRECT["brute_tests_total"]: raise SystemExit("ZetaMine disagrees with exhaustive enumeration.")

if HAVE_FIM:
    log(f"== Correctness (2/2): {UNIT_TESTS} random DBs vs PyFIM FP-Growth (+ documented omission) ==")
    rng = np.random.default_rng(12345); passed = 0; tot = 0
    for t in range(UNIT_TESTS):
        n = int(rng.integers(30, 400)); m = int(rng.integers(4, 19))
        X = rng.random((n, m)) < rng.uniform(0.1, 0.95, size=m)
        if t % 4 == 0: X[:, 0] = True                  # one always-present item
        if t % 8 == 0: X[:, 1] = True                  # two always-present items
        path = os.path.join(WORK, "tmp_db.dat")
        with open(path, "w") as f:
            for r in X:
                idx = np.flatnonzero(r)
                if len(idx): f.write(" ".join(str(int(i) + 1) for i in idx) + "\n")
        tr = load_tr(path)
        if not tr:
            print(f"  [{t+1:>2}/{UNIT_TESTS}] skipped (empty db)"); continue
        ms = int(rng.integers(2, max(3, len(tr) // 2)))
        ok, refn, nadd = check_vs_pyfim(path, tr, ms, f"unit{t}")
        tot += 1; passed += int(ok); CORRECT["pyfim_omitted_itemsets_total"] += nadd
        print(f"  [{tot:>2}/{UNIT_TESTS}] n={len(tr):4d} m={m:2d} minsup={ms:4d}  "
              f"itemsets={refn:6d}  pyfim_omitted={nadd}  {'PASS' if ok else 'FAIL'}", flush=True)
    CORRECT.update(unit_tests_passed=passed, unit_tests_total=tot)
    log(f"unit tests: {passed}/{tot}")
    if passed != tot: raise SystemExit("Correctness FAILED -- do not trust timing below.")

# ---------------------------- 6. Datasets ---------------------------------------------
def candidates(name):
    out = list(VERIFIED_MIRRORS.get(name, []))
    for base in ("http://fimi.uantwerpen.be/data/", "https://fimi.uantwerpen.be/data/"):
        out += [base + name + ".dat", base + name + ".dat.gz"]
    out.append("https://www.philippe-fournier-viger.com/spmf/datasets/" + name + ".txt")
    out += USER_MIRRORS.get(name, [])
    return out

def valid_file(p):
    try:
        with open(p, encoding="utf-8", errors="ignore") as f:
            for _ in range(3):
                line = f.readline()
                if line and not all(x.isdigit() for x in line.split()): return False
        return os.path.getsize(p) > 1000
    except Exception: return False

def get_dataset(name):
    path = os.path.join(DATA, name + ".dat")
    if valid_file(path):
        log(f"  {name}: already present"); return path
    for url in candidates(name):
        try:
            raw = fetch(url, timeout=120)
            open(path, "wb").write(raw)
            if valid_file(path):
                log(f"  {name}: downloaded from {url}"); return path
            os.remove(path)
        except Exception:
            continue
    log(f"  {name}: !! could not download")
    return None

def ds_stats(path):
    n = tot = 0; cnt = collections.Counter()
    with open(path) as f:
        for line in f:
            t = set(line.split())
            if t: n += 1; tot += len(t); cnt.update(t)
    nfull = sum(1 for c in cnt.values() if c == n)
    return n, len(cnt), tot / n, nfull

log("== Datasets ==")
DS = {}
for name in PRESETS[PRESET]:
    p = get_dataset(name)
    if p:
        n, ni, al, nfull = ds_stats(p)
        DS[name] = dict(path=p, n=n, items=ni, avg_len=al, density=al / ni, n_full=nfull)
        log(f"  {name:12s} n={n:7d} items={ni:6d} avg_len={al:6.1f} density={al/ni:6.1%} "
            f"items_in_every_transaction={nfull}")

if not DS:
    log("!! No datasets downloaded. Generating synthetic (NOT valid for a paper).")
    def synth_zipf(nm, n, m, avg_len, seed):
        r = np.random.default_rng(seed)
        w = 1.0 / (np.arange(1, m + 1) ** 1.1); w /= w.sum()
        p = os.path.join(DATA, nm + ".dat")
        with open(p, "w") as f:
            for _ in range(n):
                L = max(1, min(m, int(r.poisson(avg_len))))
                idx = r.choice(m, size=L, replace=False, p=w)
                f.write(" ".join(str(int(i) + 1) for i in sorted(idx)) + "\n")
        n_, ni, al, nf = ds_stats(p); DS[nm] = dict(path=p, n=n_, items=ni, avg_len=al, density=al / ni, n_full=nf)
    synth_zipf("synthetic-zipf", 50000, 500, 15, 11)
    PRESETS[PRESET]["synthetic-zipf"] = [0.10, 0.05, 0.02]

if HAVE_FIM:
    log("== Correctness on real data (full itemset sets at highest support) ==")
    for name, d in DS.items():
        rel = max(PRESETS[PRESET][name]); ms = math.ceil(rel * d["n"] - 1e-9)
        c = run_zm(d["path"], ms, ZM_MAIN)
        if c["status"] == "ok" and c["count"] <= REAL_CHECK_MAX:
            log(f"  {name} @ {rel:.3%} (minsup={ms}): ZM found {c['count']} itemsets, checking ...")
            ok, cnt, nadd = check_vs_pyfim(d["path"], load_tr(d["path"]), ms, name)
            CORRECT["real_data_verified"][name] = dict(identical=bool(ok), itemsets=cnt, pyfim_omitted=nadd,
                                                       items_in_every_transaction=d["n_full"])
            CORRECT["pyfim_omitted_itemsets_total"] += nadd
            log(f"  {name}: {cnt} itemsets, identical to PyFIM(+omission): {ok}  "
                f"[PyFIM omitted {nadd}; dataset has {d['n_full']} item(s) in every transaction -> expected {2**d['n_full']-1 if d['n_full'] else 0}]")
            if not ok: raise SystemExit(f"Correctness FAILED on {name}.")
        else:
            log(f"  {name}: skipped (count={c.get('count')})")

# ---------------------------- 7. Benchmark --------------------------------------------
def med(reps, key):
    v = [r[key] for r in reps if r.get(key) is not None]
    return float(np.median(v)) if v else None

log(f"== Benchmark (median of up to {REPS} runs) ==")
RES = []; dead = {}
for ds_idx, (name, d) in enumerate(DS.items(), start=1):
    methods = make_methods(d["n"], d["n_full"]); dead.setdefault(name, set())
    log(f"[{ds_idx}/{len(DS)}] dataset: {name}  n={d['n']}  methods={len(methods)}")
    for rel in sorted(PRESETS[PRESET][name], reverse=True):
        ms = max(1, math.ceil(rel * d["n"] - 1e-9)); rows = {}
        log(f"   support {rel*100:g}%  (absolute minsup = {ms})")
        pre = run_zm(d["path"], ms, ZM_MAIN)            # cheap count-only pre-run, used only for the RAM guard
        too_big = pre["status"] == "ok" and pre["count"] > LIST_MAX
        if too_big: log(f"     {pre['count']:,} itemsets > ZM_LIST_MAX={LIST_MAX:,}: PyFIM list mode skipped (RAM)")
        for mname, fn in methods:
            if too_big and mname in ("PyFIM-FPGrowth", "PyFIM-Eclat"):
                rows[mname] = dict(status="skipped"); continue
            if mname in dead[name]:
                rows[mname] = dict(status="skipped"); log(f"     {mname:22s} skipped (previous timeout)"); continue
            if mname in ("PyFIM-FPGrowth", "PyFIM-Eclat") and rows.get(ZM_CNT, {}).get("status") == "ok" \
                    and rows[ZM_CNT]["count"] > LIST_MAX:
                rows[mname] = dict(status="skipped")
                log(f"     {mname:22s} skipped (list mode would hold {rows[ZM_CNT]['count']:,} itemsets in RAM; ZM_LIST_MAX={LIST_MAX:,})")
                continue
            reps, last = [], None
            for rep in range(REPS):
                r = fn(d["path"], ms); last = r
                if r["status"] != "ok":
                    log(f"     {mname:22s} rep {rep+1}/{REPS}  -> {r['status']} {r.get('msg','')[:120]}"); break
                reps.append(r)
                log(f"     {mname:22s} rep {rep+1}/{REPS}  t={r['t']:.4f}s wall={r['wall']:.3f}s (count={r['count']})")
                if r["t"] > SLOW_S:
                    log(f"       slow run > {SLOW_S}s; not repeating"); break
            if last["status"] == "ok":
                last = dict(last)
                for key in ("t", "wall", "t_parse", "t_prep", "t_mine"):
                    if med(reps, key) is not None: last[key] = med(reps, key)
                log(f"     {mname:22s} median t={last['t']:.4f}s wall={last['wall']:.3f}s (count={last['count']})")
            elif last["status"] == "timeout":
                dead[name].add(mname)
            rows[mname] = last
        ref = rows[ZM_OUT]["count"] if rows[ZM_OUT]["status"] == "ok" else next(
            (r["count"] for r in rows.values() if r["status"] == "ok" and r.get("count") is not None), None)
        for mname, r in rows.items():
            mismatch = (r["status"] == "ok" and ref is not None and r.get("count") is not None and r["count"] != ref)
            if mismatch: log(f"     !! COUNT MISMATCH {mname}: {r['count']} vs reference {ref}")
            RES.append(dict(dataset=name, rel_support=rel, minsup=ms, method=mname, status=r["status"],
                            count=r.get("count"), count_raw=r.get("count_raw"),
                            t_median=r.get("t"), wall_median=r.get("wall"),
                            t_parse=r.get("t_parse"), t_prep=r.get("t_prep"), t_mine=r.get("t_mine"),
                            count_mismatch=mismatch, kernel_calls=r.get("kernel_calls"),
                            kernel_sets=r.get("kernel_sets"), dfs_nodes=r.get("dfs_nodes")))
        parts = [f"{m}={rows[m]['t']:.3f}s" if rows[m]["status"] == "ok" else f"{m}={rows[m]['status']}"
                 for m in [ZM_OUT, ZM_CNT] + [b for b in ALL_BASE if b in rows]]
        log(f"   -> {name} sup={rel:<7.4f} n={ref if ref else '?':>10}  " + "  ".join(parts))

df = pd.DataFrame(RES); df.to_csv("results.csv", index=False)
log(f"wrote results.csv ({len(df)} rows)")
PRESENT = set(df.method.unique())

# ---------------------------- 8. Analysis ---------------------------------------------
def tget(ds, rel, m, col="t_median"):
    s = df[(df.dataset == ds) & (df.rel_support == rel) & (df.method == m)]
    if len(s) and s.iloc[0].status == "ok" and not s.iloc[0].count_mismatch and pd.notna(s.iloc[0][col]):
        return float(s.iloc[0][col])
    return None

cfgs = df[["dataset", "rel_support"]].drop_duplicates().values.tolist()
def speedups(zm, base_names, label, col="t_median"):
    base_names = [b for b in base_names if b in PRESENT]
    if not base_names: print(f"  {label}: no baselines available"); return None
    sp, wins, tos = [], 0, 0
    for ds, rel in cfgs:
        z = tget(ds, rel, zm, col)
        bs = [b for b in (tget(ds, rel, x, col) for x in base_names) if b is not None]
        if z is None: continue
        if not bs:
            if len(df[(df.dataset == ds) & (df.rel_support == rel) & (df.method.isin(base_names)) &
                      (df.status.isin(["timeout", "skipped"]))]): tos += 1
            continue
        s = min(bs) / z; sp.append((ds, rel, s)); wins += int(s > 1)
    if not sp: print(f"  {label}: no comparable configs"); return None
    v = np.array([x[2] for x in sp])
    print(f"  {label:44s} configs={len(v):3d}  ZM faster in {wins:3d}  slower in {len(v)-wins:3d}  "
          f"geo-mean={np.exp(np.log(v).mean()):6.2f}x  median={np.median(v):6.2f}x  min={v.min():.2f}x  max={v.max():.2f}x"
          + (f"  [+{tos} configs where baselines timed out]" if tos else ""))
    return sp

log("== Speedup summary (speedup = fastest like-for-like baseline / ZetaMine; >1 means ZetaMine faster) ==")
print("  -- OUTPUT mode: every method emits all itemsets (ZM/SPMF: text to null device; PyFIM: in-memory list)")
SP_OUT  = speedups(ZM_OUT, OUT_BASE, "t    (load+prep+mine+output)  vs output baselines", "t_median")
speedups(ZM_OUT, OUT_BASE, "wall (whole process)              vs output baselines", "wall_median")
print("  -- COUNT-ONLY mode: ZetaMine-count vs PyFIM report='#' (and native tools if configured)")
SP_CNT  = speedups(ZM_CNT, CNT_BASE, "t    (load+prep+mine)         vs count baselines", "t_median")
speedups(ZM_CNT, CNT_BASE, "wall (whole process)              vs count baselines", "wall_median")
for ds in DS:
    sub = [(d_, r_, s_) for d_, r_, s_ in (SP_OUT or []) if d_ == ds]
    if sub:
        v = np.array([s for _, _, s in sub])
        log(f"    [output] {ds:12s} geo-mean {np.exp(np.log(v).mean()):6.2f}x   per-support: " +
            ", ".join(f"{r:g}:{s:.2f}x" for _, r, s in sub))
    sub = [(d_, r_, s_) for d_, r_, s_ in (SP_CNT or []) if d_ == ds]
    if sub:
        v = np.array([s for _, _, s in sub])
        log(f"    [count ] {ds:12s} geo-mean {np.exp(np.log(v).mean()):6.2f}x   per-support: " +
            ", ".join(f"{r:g}:{s:.2f}x" for _, r, s in sub))

log("== ZetaMine time breakdown (ms): where does the time go? ==")
print(f"  {'dataset':12s} {'sup':>7s} {'parse':>9s} {'prep':>9s} {'mine':>9s} {'prep share of (prep+mine)':>26s}")
for ds, rel in cfgs:
    s = df[(df.dataset == ds) & (df.rel_support == rel) & (df.method == ZM_CNT) & (df.status == "ok")]
    if len(s):
        r = s.iloc[0]
        print(f"  {ds:12s} {rel*100:6g}% {r.t_parse*1e3:9.1f} {r.t_prep*1e3:9.1f} {r.t_mine*1e3:9.1f} "
              f"{100*r.t_prep/max(r.t_prep + r.t_mine, 1e-12):25.0f}%")

if RUN_ABLATIONS:
    log("== Ablation: slowdown relative to ZetaMine-count (same mode, same metric t) ==")
    for ds, rel in cfgs:
        z = tget(ds, rel, ZM_CNT)
        if z:
            print(f"  {ds:12s} {rel*100:6g}%  " + "  ".join(
                f"{a}={tget(ds, rel, a)/z:6.2f}x" if tget(ds, rel, a) else f"{a}=n/a" for a in ABLATIONS))

# ---------------------------- 9. LaTeX + plots + environment --------------------------
nm = {ZM_OUT: "ZetaMine", ZM_CNT: "ZetaMine (count)",
      "PyFIM-FPGrowth": "FP-Growth", "PyFIM-Eclat": "Eclat",
      "PyFIM-FPGrowth-count": "FP-Growth (count)", "PyFIM-Eclat-count": "Eclat (count)",
      "SPMF-LCM": "LCM (SPMF)", "SPMF-FPGrowth": "FP-Growth (SPMF)", "SPMF-Eclat": "Eclat (SPMF)"}
for k in EXTERNAL_BASELINES: nm.setdefault(k, k)

def fmt(ds, rel, m, best, col="t_median"):
    s = df[(df.dataset == ds) & (df.rel_support == rel) & (df.method == m)]
    if not len(s): return "--"
    r = s.iloc[0]
    if r.status == "timeout": return "TO"
    if r.status != "ok" or pd.isna(r[col]): return "--"
    txt = f"{r[col]:.3f}" + (r"$^\dagger$" if r.count_mismatch else "")
    return r"\textbf{" + txt + "}" if (best is not None and abs(r[col] - best) < 1e-12 and not r.count_mismatch) else txt

def write_table(zm, bases, fname, label, caption_extra):
    shown = [m for m in [zm] + bases if m in PRESENT]
    quirk = [f"{d}: {v['n_full']} item(s) -> {2**v['n_full']-1} itemset(s)" for d, v in DS.items() if v["n_full"]]
    qtxt = (r" Borgelt's \texttt{fim} (PyFIM) does not report itemsets consisting only of items present in every "
            r"transaction (" + "; ".join(quirk) + r"); its counts are corrected by that amount and verified "
            r"itemset-by-itemset against ZetaMine and exhaustive enumeration.") if quirk else ""
    tex = [r"\begin{table*}[t]", r"\centering", r"\small",
           r"\caption{Runtime in seconds (median of up to " + str(REPS) + r" runs, single thread, best per row in bold). "
           + caption_extra +
           r" ZetaMine and SPMF times include reading the input file; PyFIM times cover the library call on an "
           r"already-loaded in-memory database, which favours PyFIM." + qtxt + r"}",
           r"\label{" + label + r"}", r"\begin{tabular}{llr" + "r" * len(shown) + r"r}", r"\toprule",
           "Dataset & minsup & \\#itemsets & " + " & ".join(nm[m] for m in shown) + r" & Speedup \\", r"\midrule"]
    for ds in DS:
        first = True
        for rel in sorted(PRESETS[PRESET][ds], reverse=True):
            tv = {m: tget(ds, rel, m) for m in shown}
            okv = [v for v in tv.values() if v is not None]
            best = min(okv) if okv else None
            c = df[(df.dataset == ds) & (df.rel_support == rel) & (df.method == zm) & (df.status == "ok")]["count"]
            c = f"{int(c.iloc[0]):,}".replace(",", r"\,") if len(c) else "--"
            z = tv.get(zm); bs = [tv[b] for b in bases if b in tv and tv[b] is not None]
            sp = f"{min(bs)/z:.2f}$\\times$" if (z and bs) else "--"
            tex.append(f"{ds if first else ''} & {rel*100:g}\\% & {c} & " +
                       " & ".join(fmt(ds, rel, m, best) for m in shown) + f" & {sp} \\\\")
            first = False
        tex.append(r"\midrule")
    tex[-1] = r"\bottomrule"; tex += [r"\end{tabular}", r"\end{table*}"]
    open(fname, "w").write("\n".join(tex)); log(f"wrote {fname}")

write_table(ZM_OUT, OUT_BASE, "table_runtime_output.tex", "tab:runtime-output",
            r"\textbf{Output mode}: every method enumerates all frequent itemsets with supports "
            r"(ZetaMine, SPMF: text written to the null device; PyFIM: itemsets materialised as an in-memory list).")
write_table(ZM_CNT, CNT_BASE, "table_runtime_count.tex", "tab:runtime-count",
            r"\textbf{Count-only mode}: no itemsets are materialised or written "
            r"(PyFIM \texttt{report='\#'}).")

def make_plot(methods_, fname):
    ncol = 3; nr = math.ceil(len(DS) / ncol)
    fig, axes = plt.subplots(nr, ncol, figsize=(5 * ncol, 3.6 * nr), squeeze=False)
    style = {ZM_OUT: ("o-", "C3"), ZM_CNT: ("o--", "C3"), "PyFIM-FPGrowth": ("s-", "C0"), "PyFIM-Eclat": ("^-", "C2"),
             "PyFIM-FPGrowth-count": ("s--", "C0"), "PyFIM-Eclat-count": ("^--", "C2"),
             "SPMF-LCM": ("d-", "C1"), "SPMF-FPGrowth": ("v-", "C4"), "SPMF-Eclat": ("P-", "C5"),
             "ZM-noKernel": ("x--", "gray")}
    for ax, ds in zip(axes.flatten(), DS):
        for m in [k for k in methods_ if k in PRESENT]:
            pts = [(rel, tget(ds, rel, m)) for rel in sorted(PRESETS[PRESET][ds], reverse=True)]
            pts = [(a, b) for a, b in pts if b is not None]
            st = style.get(m, ("*-", "k"))
            if pts: ax.plot([a * 100 for a, _ in pts], [b for _, b in pts], st[0], color=st[1], label=nm.get(m, m))
        ax.set_yscale("log"); ax.set_title(ds); ax.set_xlabel("minsup (%)"); ax.set_ylabel("time (s)")
        ax.invert_xaxis(); ax.grid(alpha=.3)
    for ax in axes.flatten()[len(DS):]: ax.axis("off")
    axes.flatten()[0].legend(fontsize=8); plt.tight_layout()
    plt.savefig(fname + ".png", dpi=300); plt.savefig(fname + ".pdf"); plt.close(fig)
    log(f"wrote {fname}.png, {fname}.pdf")

make_plot([ZM_OUT] + OUT_BASE, "runtime_output")
make_plot([ZM_CNT, "ZM-noKernel"] + CNT_BASE, "runtime_count")

def cpu_name():
    try: return platform.processor() or os.environ.get("PROCESSOR_IDENTIFIER", "unknown")
    except Exception: return "unknown"
try: gpp_ver = sh([gpp_bin, "--version"]).stdout.splitlines()[0]
except Exception: gpp_ver = "unknown"
env = dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), cpu=cpu_name(), cores=os.cpu_count(),
           os=platform.platform(), python=sys.version.split()[0], gpp=gpp_ver,
           compile_flags="-O3 -std=c++17 -DNDEBUG", preset=PRESET, reps=REPS, timeout=TIMEOUT,
           zm_main=ZM_MAIN,
           timing_definitions=dict(
               t_ZetaMine="parse + prep + mine (+ formatted output to null device in output mode)",
               t_SPMF="SPMF-reported total time (includes file read and JIT warm-up)",
               t_PyFIM="fim call only, DB already in memory (list mode materialises tuples; count mode uses report='#')",
               wall="whole subprocess wall-clock incl. process/JVM/interpreter start-up and file parsing"),
           pyfim_note="fim omits all itemsets consisting only of items present in every transaction "
                      "(2^c-1 for c such items); corrected and verified, not a float-boundary effect",
           datasets={k: {kk: vv for kk, vv in v.items() if kk != "path"} for k, v in DS.items()},
           correctness=CORRECT, spmf=HAVE_SPMF, pyfim=HAVE_FIM, external=list(EXTERNAL_BASELINES))
json.dump(env, open("environment.json", "w"), indent=2, default=str)
log("wrote environment.json")

shutil.make_archive(os.path.join(os.path.dirname(WORK), "zetamine_results"), "zip", WORK, ".")
log(f"Done. Results in: {WORK}")