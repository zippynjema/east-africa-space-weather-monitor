#!/usr/bin/env python3
"""
GOES X-ray flare alert dashboard + flare -> geomagnetic impact latency
=====================================================================

What it does
------------
1. Pulls GOES primary XRS long-channel (0.1-0.8 nm) flux from NOAA SWPC.
2. Detects flares (A/B/C/M/X) and fires alerts (console, log file, optional
   Slack-style webhook) when a flare first crosses M1, M5, X1, X10 (tunable).
3. Pulls real-time L1 solar wind, Kp and Dst from SWPC.
4. Kenya / East Africa HF-radio view: is the region in daylight when a flare peaks, and what
   highest HF frequency is likely degraded (NOAA D-RAP method). Also shown in alert text.
5. Computes latency from flare peak to geomagnetic impact two ways:
     A) DONKI-linked chain  FLR -> CME -> (IPS) -> GST   (best attribution)
     B) Heuristic timing    flare peak -> L1 shock / Kp >= N / Dst <= X
6. Writes a self-refreshing HTML dashboard + CSV reports.

Install
-------
    pip install -r requirements.txt        # numpy, pandas, requests, plotly, python-dotenv (optional)
    Optional settings (email / webhook alerts) go in a .env file next to this script: copy .env.example to .env.

Run
---
    python goes_flare_alert_dashboard.py                 # live mode (default): refreshes every 60 s, opens the browser
    python goes_flare_alert_dashboard.py --once          # single pass, then exit
    python goes_flare_alert_dashboard.py --webhook https://hooks.slack.com/...   # also push alerts to Slack
    # The page reloads itself every minute, but the script must stay running to produce fresh data.

Data sources
------------
    SWPC GOES X-ray JSON : https://services.swpc.noaa.gov/json/goes/primary/xrays-{6-hour,1-day,3-day,7-day}.json
    SWPC solar wind      : https://services.swpc.noaa.gov/json/rtsw/rtsw_{wind,mag}_1m.json (~24 h; cached locally)
                           (old products/solar-wind/* files were withdrawn, SWPC notice 26-21)
    SWPC Kp / Dst        : .../json/planetary_k_index_1m.json, .../products/noaa-planetary-k-index.json,
                           .../products/kyoto-dst.json
    NASA CCMC DONKI      : https://ccmc.gsfc.nasa.gov/DONKI-API/get/{FLR,CME,IPS,GST}  (case-sensitive, no key)

Physics caveat (important for interpreting "latency")
-----------------------------------------------------
Flare X-rays reach Earth in ~8.3 min and cause immediate radio blackouts, but
*geomagnetic storms* come from CMEs (or high-speed streams), which take
roughly 15 h to 4+ days, and most flares have no Earth-directed CME at all.
So "flare -> geomagnetic impact" is only meaningful for flares that actually
launched an Earth-directed CME. DONKI links (method A) are the reliable
attribution; the heuristic (method B) can mis-associate and reports how many
competing flares share the same arrival window.
"""

import argparse
import json
import os
import re
import smtplib
import ssl
import sys
import time
import webbrowser
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from email.message import EmailMessage
from pathlib import Path

import numpy as np
import pandas as pd
import requests

try:                                # optional: read settings from a local .env file (see .env.example)
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).resolve().parent / ".env")   # .env next to this script
    load_dotenv()                                            # or in the current working folder
except ImportError:                 # python-dotenv not installed: plain environment variables still work
    pass

SWPC = "https://services.swpc.noaa.gov"
DONKI = "https://ccmc.gsfc.nasa.gov/DONKI-API/get"   # moved 2026-09-30; no API key needed
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "flare-alert-dashboard/1.0 (research script)"

CLASS_BASE = {"A": 1e-8, "B": 1e-7, "C": 1e-6, "M": 1e-5, "X": 1e-4}
R_SCALE = [(2e-3, "R5"), (1e-3, "R4"), (1e-4, "R3"), (5e-5, "R2"), (1e-5, "R1")]
FLARE_COLORS = {"A": "#9aa0a6", "B": "#9aa0a6", "C": "#b7950b", "M": "#e67e22", "X": "#c0392b"}


def log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S}Z] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def get_json(url, params=None, retries=3, timeout=(10, 45)):
    for i in range(retries):
        try:
            r = SESSION.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, ValueError) as exc:
            if i == retries - 1:
                log(f"FETCH FAILED {url}: {exc}")
                return None
            time.sleep(2 ** i)


def to_frame(raw) -> pd.DataFrame:
    """SWPC serves either [header_row, *rows] or a list of dicts."""
    if not raw:
        return pd.DataFrame()
    if isinstance(raw[0], (list, tuple)):
        return pd.DataFrame(raw[1:], columns=raw[0])
    return pd.DataFrame(raw)


def flare_class(flux: float) -> str:
    for letter, lo, hi in (("A", 0, 1e-7), ("B", 1e-7, 1e-6), ("C", 1e-6, 1e-5),
                           ("M", 1e-5, 1e-4), ("X", 1e-4, np.inf)):
        if flux < hi:
            return f"{letter}{flux / CLASS_BASE[letter]:.1f}"
    return "?"


def class_to_flux(cls: str) -> float:
    m = re.match(r"^([ABCMX])(\d+(?:\.\d+)?)$", cls.strip().upper())
    if not m:
        raise ValueError(f"Bad flare class '{cls}' (expected e.g. M1.0)")
    return CLASS_BASE[m.group(1)] * float(m.group(2))


def r_scale(flux: float) -> str:
    return next((r for lim, r in R_SCALE if flux >= lim), "-")


def hours(a, b) -> float:
    return round((b - a).total_seconds() / 3600, 2) if pd.notna(a) and pd.notna(b) else np.nan


def first_time(s: pd.Series, cond, t0, t1):
    w = s[(s.index >= t0) & (s.index <= t1)]
    hit = w[cond(w)]
    return hit.index[0] if len(hit) else pd.NaT


# --------------------------------------------------------------------------- #
# Fetchers
# --------------------------------------------------------------------------- #
def fetch_xray(window: str) -> pd.DataFrame:
    df = to_frame(get_json(f"{SWPC}/json/goes/primary/xrays-{window}.json"))
    if df.empty:
        return df
    if "energy" in df:
        df = df[df["energy"] == "0.1-0.8nm"]
    col = "flux" if "flux" in df else "observed_flux"
    df = df.assign(
        time=pd.to_datetime(df["time_tag"], utc=True, errors="coerce"),
        flux=pd.to_numeric(df[col], errors="coerce"),
    )
    df = df.dropna(subset=["time", "flux"])
    df = df[df["flux"] > 0].sort_values("time").drop_duplicates("time").set_index("time")
    keep = ["flux"] + (["satellite"] if "satellite" in df else [])
    return df[keep]


def _fetch_rtsw(name: str) -> pd.DataFrame:
    """New RTSW feeds (list of dict records, newest first, several spacecraft mixed)."""
    for base in (SWPC, "https://noaa-swpc-pds.s3.amazonaws.com"):      # S3 mirror as fallback
        d = to_frame(get_json(f"{base}/json/rtsw/{name}", retries=2))
        if not d.empty:
            return d
    return pd.DataFrame()


def fetch_solar_wind(cache_path: Path) -> pd.DataFrame:
    """
    SWPC retired products/solar-wind/*-7-day.json (404). Replacement feeds hold only ~24 h,
    so each poll is merged into a local CSV cache (SWPC's own advice: retain the 1-day file).
    History therefore builds up the longer the script runs (kept for 30 days).
    """
    parts = []
    for name, cols in (("rtsw_wind_1m.json", {"proton_speed": "speed", "proton_density": "density"}),
                       ("rtsw_mag_1m.json", {"bz_gsm": "bz_gsm", "bt": "bt"})):
        d = _fetch_rtsw(name)
        if d.empty or "time_tag" not in d:
            continue
        if "active" in d:                                   # keep the spacecraft NOAA treats as primary
            act = d[d["active"].astype(str).str.lower().isin(["true", "1"])]
            if not act.empty:
                d = act
        d = d.assign(time=pd.to_datetime(d["time_tag"], utc=True, errors="coerce"))
        d = d.dropna(subset=["time"]).set_index("time").sort_index().rename(columns=cols)
        d = d[[c for c in cols.values() if c in d]].apply(pd.to_numeric, errors="coerce")
        parts.append(d.resample("1min").median())           # also merges second-offsets between spacecraft
    new = pd.concat(parts, axis=1) if parts else pd.DataFrame()

    old = pd.DataFrame()
    if cache_path.exists():
        old = pd.read_csv(cache_path, index_col=0)
        old.index = pd.to_datetime(old.index, utc=True, errors="coerce")
        old = old[old.index.notna()]
    sw = pd.concat([old, new]).sort_index()
    if sw.empty:
        return sw
    sw = sw[~sw.index.duplicated(keep="last")].dropna(how="all")
    sw = sw[sw.index >= pd.Timestamp.now(tz="UTC") - timedelta(days=30)]
    sw.to_csv(cache_path)
    return sw


def fetch_kp() -> pd.Series:
    parts = []
    d = to_frame(get_json(f"{SWPC}/json/planetary_k_index_1m.json"))        # ~1-min estimated Kp
    if not d.empty and "time_tag" in d:
        col = "estimated_kp" if "estimated_kp" in d else "kp_index"
        parts.append(pd.Series(pd.to_numeric(d[col], errors="coerce").values,
                               index=pd.to_datetime(d["time_tag"], utc=True, errors="coerce")))
    d = to_frame(get_json(f"{SWPC}/products/noaa-planetary-k-index.json"))  # 3-hourly Kp
    if not d.empty:
        col = next((c for c in d.columns if str(c).lower() in ("kp", "kp_index")), None)
        if col:
            parts.append(pd.Series(pd.to_numeric(d[col], errors="coerce").values,
                                   index=pd.to_datetime(d["time_tag"], utc=True, errors="coerce")))
    if not parts:
        return pd.Series(dtype=float)
    s = pd.concat(parts).dropna()
    s = s[s.index.notna()].sort_index()
    return s[~s.index.duplicated(keep="last")]


def fetch_dst() -> pd.Series:
    d = to_frame(get_json(f"{SWPC}/products/kyoto-dst.json"))
    if d.empty:
        return pd.Series(dtype=float)
    tcol = "time_tag" if "time_tag" in d else d.columns[0]
    vcol = "dst" if "dst" in d else d.columns[1]
    s = pd.Series(pd.to_numeric(d[vcol], errors="coerce").values,
                  index=pd.to_datetime(d[tcol], utc=True, errors="coerce")).dropna()
    s = s[s.index.notna()].sort_index()
    return s[~s.index.duplicated()]


# --------------------------------------------------------------------------- #
# Flare + shock detection
# --------------------------------------------------------------------------- #
FLARE_COLS = ["start", "end", "peak_time", "peak_flux", "cls", "in_progress"]


def detect_flares(xr: pd.DataFrame, threshold=1e-6, merge_gap_min=10, min_dur_min=1) -> pd.DataFrame:
    """Contiguous periods with flux >= threshold (default C1), merged across short gaps."""
    s = xr["flux"]
    above = s >= threshold
    if not above.any():
        return pd.DataFrame(columns=FLARE_COLS)
    grp = (above != above.shift()).cumsum()
    ev = [[seg.index[0], seg.index[-1], seg.idxmax(), seg.max()]
          for _, seg in s[above].groupby(grp[above])]
    merged: list = []
    for e in ev:
        if merged and (e[0] - merged[-1][1]) <= timedelta(minutes=merge_gap_min):
            m = merged[-1]
            m[1] = e[1]
            if e[3] > m[3]:
                m[2], m[3] = e[2], e[3]
        else:
            merged.append(list(e))
    df = pd.DataFrame(merged, columns=["start", "end", "peak_time", "peak_flux"])
    df = df[(df["end"] - df["start"]) >= timedelta(minutes=min_dur_min)].copy()
    df["cls"] = df["peak_flux"].map(flare_class)
    df["in_progress"] = df["end"] >= (s.index[-1] - timedelta(minutes=2))
    return df.reset_index(drop=True)[FLARE_COLS]


def detect_shocks(sw: pd.DataFrame, jump_kms=50.0, merge_min=60) -> list:
    """Sudden solar-wind speed jump at L1 (CME shock candidate). Detection lags ~3 min."""
    if sw.empty or "speed" not in sw:
        return []
    s = sw["speed"].resample("1min").median().interpolate(limit=5)
    after = s.rolling(5, min_periods=3).median()
    before = s.shift(5).rolling(10, min_periods=5).median()
    hits = (after - before) >= jump_kms
    shocks, prev = [], None
    for t in hits[hits].index:
        if prev is None or (t - prev) > timedelta(minutes=merge_min):
            shocks.append(t)
        prev = t
    return shocks


# --------------------------------------------------------------------------- #
# Alerts
# --------------------------------------------------------------------------- #
_email_warned = {"done": False}


def send_email(subject: str, body: str, args) -> bool:
    """
    Email alert via SMTP (Gmail by default). Credentials come from environment variables, never the script:
      FLARE_SMTP_PASSWORD  Gmail *App Password* (16 chars; needs 2-Step Verification on the Google account)
      FLARE_SMTP_USER      sender login (defaults to the --email-to address)
    Email is opt-in: nothing is sent unless --email-to / FLARE_EMAIL_TO is given.
      FLARE_SMTP_HOST / FLARE_SMTP_PORT   default smtp.gmail.com / 587
    """
    to = getattr(args, "email_to", None)
    pw = (os.environ.get("FLARE_SMTP_PASSWORD") or "").replace(" ", "")
    if not to or not pw:
        if not _email_warned["done"]:
            log("Email alerts are OFF. To enable: pass --email-to <address> (or set FLARE_EMAIL_TO) and set FLARE_SMTP_PASSWORD (Gmail App Password).")
            _email_warned["done"] = True
        return False
    user = os.environ.get("FLARE_SMTP_USER") or to
    host = os.environ.get("FLARE_SMTP_HOST", "smtp.gmail.com")
    port = int(os.environ.get("FLARE_SMTP_PORT", "587"))
    msg = EmailMessage()
    msg["Subject"], msg["From"], msg["To"] = subject, user, to
    msg.set_content(body)
    try:
        ctx = ssl.create_default_context()
        if port == 465:
            with smtplib.SMTP_SSL(host, port, context=ctx, timeout=30) as srv:
                srv.login(user, pw)
                srv.send_message(msg)
        else:
            with smtplib.SMTP(host, port, timeout=30) as srv:
                srv.starttls(context=ctx)
                srv.login(user, pw)
                srv.send_message(msg)
        log(f"Email sent to {to}: {subject}")
        return True
    except (smtplib.SMTPException, OSError) as exc:
        log(f"Email FAILED ({type(exc).__name__}): {exc}")
        return False


def process_alerts(flares: pd.DataFrame, args, out_dir: Path) -> None:
    state_path = out_dir / "alert_state.json"
    seed = not state_path.exists()
    state = json.loads(state_path.read_text()) if not seed else {}
    min_flux = class_to_flux(args.alert_class)
    levels = sorted({min_flux, *[x for x in (1e-5, 5e-5, 1e-4, 1e-3) if x > min_flux]})
    webhook = args.webhook or os.environ.get("FLARE_WEBHOOK_URL")

    for f in flares.itertuples():
        lvl = sum(f.peak_flux >= x for x in levels)
        key = f.start.isoformat()
        if lvl == 0:
            continue
        if seed and not f.in_progress:      # don't flood on first run with history
            state[key] = lvl
            continue
        if lvl > state.get(key, 0):
            state[key] = lvl
            tag = "IN PROGRESS (peak so far)" if f.in_progress else "peak"
            msg = (f"SOLAR FLARE {f.cls} | {tag} {f.peak_time:%Y-%m-%d %H:%M}Z | "
                   f"radio-blackout scale {r_scale(f.peak_flux)} (sunlit side) | "
                   f"{kenya_note(f.peak_flux, f.peak_time)} | "
                   f"geomagnetic impact only if an Earth-directed CME was launched (>=~15 h later)")
            print("\a" + "!" * 70 + f"\n{msg}\n" + "!" * 70, flush=True)
            with open(out_dir / "alerts.log", "a") as fh:
                fh.write(f"{datetime.now(timezone.utc).isoformat()} {msg}\n")
            if webhook:
                try:
                    requests.post(webhook, json={"text": msg}, timeout=10)
                except requests.RequestException as exc:
                    log(f"webhook failed: {exc}")
            send_email(f"Solar flare {f.cls} ({r_scale(f.peak_flux)}) - {f.peak_time:%d %b %H:%M} UTC",
                       msg + f"\n\nDashboard file: {(out_dir / 'dashboard.html').resolve()}", args)
    state = dict(sorted(state.items())[-500:])
    state_path.write_text(json.dumps(state))


# --------------------------------------------------------------------------- #
# Kenya / East Africa HF radio-blackout view
# --------------------------------------------------------------------------- #
# Method: NOAA SWPC D-RAP. Highest affected frequency (1 dB absorption, vertical path) at the
# sub-solar point is HAF = 10*log10(flux[W/m2]) + 65 MHz, tapering as cos(zenith)**0.75.
# HF blackouts only occur on the sunlit side because the X-rays must reach the D-region.
HF_SITES = {                       # name: (lat, lon)
    "Nairobi": (-1.286, 36.817), "Mombasa": (-4.043, 39.668), "Kisumu": (-0.092, 34.768),
    "Lodwar": (3.119, 35.597), "Kampala": (0.347, 32.582),
    "Dar es Salaam": (-6.792, 39.208), "Addis Ababa": (9.030, 38.740),
}
KENYA_SITES = ("Nairobi", "Mombasa", "Kisumu", "Lodwar")
EAT = timedelta(hours=3)           # East Africa Time = UTC+3


def sun_elevation(t, lat: float, lon: float) -> float:
    """Solar elevation angle in degrees (low-precision almanac formulas, good to ~0.1 deg)."""
    t = pd.Timestamp(t)
    t = t.tz_convert("UTC").tz_localize(None) if t.tzinfo is not None else t
    n = t.to_julian_date() - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360.0
    g = np.radians((357.528 + 0.9856003 * n) % 360.0)
    lam = np.radians(L + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    eps = np.radians(23.439 - 4e-7 * n)
    dec = np.arcsin(np.sin(eps) * np.sin(lam))
    ra = np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))
    gmst_h = (18.697374558 + 24.06570982441908 * n) % 24.0
    ha = np.radians(15.0 * gmst_h + lon) - ra
    la = np.radians(lat)
    return float(np.degrees(np.arcsin(np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(ha))))


def haf_mhz(flux: float, elev_deg: float) -> float:
    """D-RAP highest affected frequency (MHz) at a site; 0 when the sun is below the horizon."""
    if elev_deg <= 0 or flux <= 0:
        return 0.0
    base = max(10.0 * np.log10(flux) + 65.0, 0.0)
    return round(base * np.sin(np.radians(elev_deg)) ** 0.75, 1)       # cos(zenith) = sin(elevation)


def day_state(elev: float) -> str:
    return "Day" if elev >= 10 else "Low sun" if elev > 0 else "Night"


def drap_condition(flux: float) -> str:
    for lim, name in ((2e-6, "Normal background"), (1e-5, "Elevated X-ray flux"), (5e-5, "Moderate"),
                      (1e-4, "Strong"), (2e-3, "Severe")):
        if flux < lim:
            return name
    return "Extreme"


def kenya_note(flux: float, t) -> str:
    """One-line daylight/HF note for alerts."""
    e = sun_elevation(t, *HF_SITES["Nairobi"])
    if e > 0:
        return f"KENYA: SUNLIT (Nairobi sun {e:.0f} deg) - HF degraded up to ~{haf_mhz(flux, e):.0f} MHz"
    return f"KENYA: NIGHT (Nairobi sun {e:.0f} deg) - little or no HF effect from this flare"


def kenya_now(flux: float) -> dict:
    now = pd.Timestamp.now(tz="UTC")
    e = sun_elevation(now, *HF_SITES["Nairobi"])
    return {"elev": e, "state": day_state(e), "haf": haf_mhz(flux, e),
            "cond": drap_condition(flux), "eat": (now + EAT).strftime("%H:%M")}


def night_spans(t0, t1, lat: float, lon: float, step_min: int = 10) -> list:
    ts = pd.date_range(pd.Timestamp(t0).floor(f"{step_min}min"), t1, freq=f"{step_min}min")
    spans, start = [], None
    for t in ts:
        dark = sun_elevation(t, lat, lon) <= 0
        if dark and start is None:
            start = t
        elif not dark and start is not None:
            spans.append((start, t))
            start = None
    if start is not None:
        spans.append((start, ts[-1]))
    return spans


def kenya_hf_table(flares: pd.DataFrame, args) -> pd.DataFrame:
    """Per-flare daylight status and estimated HF impact over Kenya / East Africa at flare peak."""
    if flares.empty:
        return pd.DataFrame()
    sel = flares[flares["peak_flux"] >= class_to_flux(args.hf_class)]
    if sel.empty:                                           # keep the section useful on quiet weeks
        sel = flares.nlargest(5, "peak_flux")
    rows = []
    for r in sel.sort_values("peak_time", ascending=False).itertuples():
        el = {name: sun_elevation(r.peak_time, *ll) for name, ll in HF_SITES.items()}
        e = el["Nairobi"]
        haf = haf_mhz(r.peak_flux, e)
        impact = (f"HF fade likely, up to ~{haf:.0f} MHz" if e >= 10 else
                  f"Sun low: reduced fade (~{haf:.0f} MHz)" if e > 0 else
                  "Night: no HF effect")
        rows.append({
            "flare_peak_utc": r.peak_time, "peak_EAT": (r.peak_time + EAT).strftime("%H:%M"),
            "class": r.cls, "R_scale": r_scale(r.peak_flux),
            "nairobi_sun_elev_deg": round(e, 1), "nairobi": day_state(e), "nairobi_HAF_MHz": haf,
            "kenya_sites_sunlit": f"{sum(el[k] > 0 for k in KENYA_SITES)}/{len(KENYA_SITES)}",
            "east_africa_sunlit": f"{sum(v > 0 for v in el.values())}/{len(HF_SITES)}",
            "impact": impact,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Latency method B: heuristic timing
# --------------------------------------------------------------------------- #
def heuristic_latency(flares, shocks, kp, dst, args) -> pd.DataFrame:
    big = flares[flares["peak_flux"] >= class_to_flux(args.assoc_class)]
    now = pd.Timestamp.now(tz="UTC")
    rows = []
    for r in big.itertuples():
        t0 = r.peak_time + timedelta(hours=args.min_h)
        t1 = r.peak_time + timedelta(hours=args.max_h)
        shock = next((t for t in shocks if t0 <= t <= t1), pd.NaT)
        kp_t = first_time(kp, lambda w: w >= args.kp, t0, t1) if len(kp) else pd.NaT
        w = dst[(dst.index >= t0) & (dst.index <= t1)] if len(dst) else dst
        dmin_t, dmin = (w.idxmin(), w.min()) if len(w) and w.min() <= args.dst else (pd.NaT, np.nan)
        ref = next((x for x in (shock, kp_t, dmin_t) if pd.notna(x)), pd.NaT)
        competing = 0
        if pd.notna(ref):
            competing = int(((big["peak_time"] >= ref - timedelta(hours=args.max_h)) &
                             (big["peak_time"] <= ref - timedelta(hours=args.min_h)) &
                             (big["peak_time"] != r.peak_time)).sum())
        rows.append({
            "flare_peak_utc": r.peak_time, "class": r.cls,
            "L1_shock_utc": shock, "peak_to_shock_h": hours(r.peak_time, shock),
            "kp_cross_utc": kp_t, f"peak_to_Kp>={args.kp:g}_h": hours(r.peak_time, kp_t),
            "dst_min_nT": dmin, "dst_min_utc": dmin_t, "peak_to_dst_min_h": hours(r.peak_time, dmin_t),
            "status": "complete" if t1 <= now else "window open",
            "competing_flares": competing,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Latency method A: DONKI-linked chain FLR -> CME -> (IPS) -> GST
# --------------------------------------------------------------------------- #
_donki_cache = {"t": 0.0, "df": pd.DataFrame()}


def _etype(i: str) -> str:
    m = re.search(r"-(FLR|CME|GST|IPS|SEP|MPC|RBE|HSS)-", i)
    return m.group(1) if m else "?"


def donki_get(kind: str, params: dict):
    """DONKI public API (CCMC). Paths are case-sensitive: FLR, CME, IPS, GST."""
    d = get_json(f"{DONKI}/{kind}", params=params, retries=2, timeout=(10, 40))
    return d if isinstance(d, list) else None


def donki_latency(args) -> pd.DataFrame:
    if args.no_donki:
        return pd.DataFrame()
    if time.time() - _donki_cache["t"] < 1800:      # refresh at most every 30 min
        return _donki_cache["df"]
    _donki_cache["t"] = time.time()          # also set on failure so a dead link isn't retried every cycle
    end = datetime.now(timezone.utc).date()
    start = end - timedelta(days=min(args.donki_days, 60))
    p = {"startDate": start.isoformat(), "endDate": end.isoformat()}
    flr = donki_get("FLR", p)
    if flr is None:
        log("DONKI unreachable or returned no JSON; skipping Latency A (will retry in 30 min).")
        return _donki_cache["df"]
    cme, ips, gst = (donki_get(k, p) or [] for k in ("CME", "IPS", "GST"))

    adj, when, kp_max = defaultdict(set), {}, {}

    def link(a, items):
        for l in items or []:
            b = l.get("activityID") if isinstance(l, dict) else l
            if b:
                adj[a].add(b)
                adj[b].add(a)

    for f in flr:
        when[f["flrID"]] = pd.to_datetime(f.get("peakTime"), utc=True, errors="coerce")
        link(f["flrID"], f.get("linkedEvents"))
    for c in cme:
        when[c["activityID"]] = pd.to_datetime(c.get("startTime"), utc=True, errors="coerce")
        link(c["activityID"], c.get("linkedEvents"))
    for i in ips:
        when[i["activityID"]] = pd.to_datetime(i.get("eventTime"), utc=True, errors="coerce")
        link(i["activityID"], i.get("linkedEvents"))
    for g in gst:
        when[g["gstID"]] = pd.to_datetime(g.get("startTime"), utc=True, errors="coerce")
        kp_max[g["gstID"]] = max((k.get("kpIndex", 0) for k in g.get("allKpIndex", []) or []), default=np.nan)
        link(g["gstID"], g.get("linkedEvents"))

    min_flux = class_to_flux(args.assoc_class)
    rows = []
    for f in flr:
        fid, cls = f["flrID"], f.get("classType") or ""
        try:
            if class_to_flux(cls) < min_flux:
                continue
        except ValueError:
            continue
        seen, frontier, found = {fid}, [(fid, [fid])], []
        for _ in range(3):                      # FLR -> CME -> (IPS) -> GST
            nxt = []
            for node, path in frontier:
                for nb in adj[node]:
                    if nb in seen:
                        continue
                    seen.add(nb)
                    t = _etype(nb)
                    if t == "GST":
                        found.append(path + [nb])
                    elif t in ("CME", "IPS"):
                        nxt.append((nb, path + [nb]))
            frontier = nxt
        if not found:
            rows.append({"flare_id": fid, "class": cls, "flare_peak_utc": when[fid], "gst_id": None})
            continue
        path = min(found, key=lambda p_: when.get(p_[-1], pd.Timestamp.max.tz_localize("UTC")))
        gid = path[-1]
        cme_id = next((x for x in path if _etype(x) == "CME"), None)
        ips_id = next((x for x in path if _etype(x) == "IPS"), None)
        rows.append({
            "flare_id": fid, "class": cls, "flare_peak_utc": when[fid],
            "cme_id": cme_id, "cme_start_utc": when.get(cme_id),
            "flare_to_cme_h": hours(when[fid], when.get(cme_id)),
            "ips_utc": when.get(ips_id), "flare_to_ips_h": hours(when[fid], when.get(ips_id)),
            "gst_id": gid, "gst_start_utc": when.get(gid), "gst_max_kp": kp_max.get(gid),
            "flare_to_gst_h": hours(when[fid], when.get(gid)),
        })
    df = pd.DataFrame(rows).sort_values("flare_peak_utc", ascending=False).reset_index(drop=True) if rows else pd.DataFrame()
    _donki_cache.update(df=df)
    return df


# --------------------------------------------------------------------------- #
# Dashboard
# --------------------------------------------------------------------------- #
# (min flux W/m2, alert level, headline, effect on HF radio, colour)
LEVELS = [
    (0.0,  "Quiet",       "Quiet Sun: no flare alert",       "No radio impact expected.", "#1e8449"),
    (1e-6, "Low",         "Low activity (C-class flare)",    "No NOAA radio-blackout level reached. HF radio is normally unaffected.", "#b7950b"),
    (1e-5, "R1 Minor",    "R1 Minor radio blackout level",   "Weak or minor degradation of HF radio on the sunlit side; occasional loss of radio contact.", "#e67e22"),
    (5e-5, "R2 Moderate", "R2 Moderate radio blackout level", "Limited HF radio blackout on the sunlit side; loss of radio contact for tens of minutes.", "#d35400"),
    (1e-4, "R3 Strong",   "R3 Strong radio blackout level",  "Wide-area HF radio blackout on the sunlit side; loss of radio contact for about an hour.", "#c0392b"),
    (1e-3, "R4 Severe",   "R4 Severe radio blackout level",  "HF radio blackout on most of the sunlit side for one to two hours.", "#922b21"),
    (2e-3, "R5 Extreme",  "R5 Extreme radio blackout level", "Complete HF radio blackout on the entire sunlit side for several hours.", "#641e16"),
]
LEVEL_CLASS_TXT = ["Below C1 (A or B class)", "C1 to C9.9", "M1 to M4.9", "M5 to M9.9",
                   "X1 to X9.9", "X10 to X19.9", "X20 and above"]

PRETTY_COLS = {
    "flare_peak_utc": "Flare peak (UTC)", "peak_EAT": "Peak time (EAT)", "class": "Class", "R_scale": "Radio blackout level",
    "nairobi_sun_elev_deg": "Nairobi sun elevation (deg)", "nairobi": "Nairobi at peak",
    "nairobi_HAF_MHz": "HF affected up to (MHz)", "kenya_sites_sunlit": "Kenyan sites in daylight",
    "east_africa_sunlit": "East African sites in daylight", "impact": "Expected impact",
    "L1_shock_utc": "Shock seen at L1 (UTC)", "peak_to_shock_h": "Flare to shock (hours)",
    "kp_cross_utc": "Kp threshold reached (UTC)", "dst_min_nT": "Lowest Dst (nT)",
    "dst_min_utc": "Time of lowest Dst (UTC)", "peak_to_dst_min_h": "Flare to lowest Dst (hours)",
    "status": "Status", "competing_flares": "Other candidate flares",
    "flare_id": "NASA flare ID", "cme_id": "Linked CME", "cme_start_utc": "CME start (UTC)",
    "flare_to_cme_h": "Flare to CME (hours)", "ips_utc": "Shock at Earth area (UTC)",
    "flare_to_ips_h": "Flare to shock (hours)", "gst_id": "Linked storm", "gst_start_utc": "Storm start (UTC)",
    "gst_max_kp": "Storm max Kp", "storm": "Geomagnetic storm", "flare_to_gst_h": "Flare to storm (hours)",
}


def alert_level(flux: float):
    i = max(k for k, lv in enumerate(LEVELS) if flux >= lv[0])
    return i, LEVELS[i]


def kp_label(k: float) -> str:
    r = int(round(k))
    if r <= 2:
        return "Quiet"
    if r == 3:
        return "Unsettled"
    if r == 4:
        return "Active"
    return {5: "G1 minor storm", 6: "G2 moderate storm", 7: "G3 strong storm", 8: "G4 severe storm"}.get(r, "G5 extreme storm")


def dst_label(d: float) -> str:
    return ("Quiet" if d > -30 else "Weak storm" if d > -50 else "Moderate storm" if d > -100
            else "Intense storm" if d > -250 else "Super storm")


def kp_color(k: float) -> str:
    r = round(k)
    return ("#1e8449" if r <= 3 else "#b7950b" if r == 4 else "#e67e22" if r == 5 else "#d35400" if r == 6
            else "#c0392b" if r == 7 else "#922b21")


def build_figure(xr, flares, sw, kp, dst, shocks, args):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    now = pd.Timestamp.now(tz="UTC")
    t_start = xr.index.min()
    titles = ("<b>Solar X-ray flux (GOES)</b>: how strong is the Sun's flare activity?",
              "<b>Solar wind speed</b> measured at L1, ~1.5 million km sunward of Earth",
              "<b>Interplanetary magnetic field (Bz)</b>: sustained negative values drive storms",
              "<b>Dst index</b>: strength of the geomagnetic storm at Earth",
              "<b>Kp index</b>: global geomagnetic activity (0 to 9; 5 or more = storm)")
    fig = make_subplots(rows=5, cols=1, shared_xaxes=True, vertical_spacing=0.07,
                        row_heights=[0.27, 0.15, 0.15, 0.17, 0.17], subplot_titles=titles)

    # 1 - X-ray flux
    fig.add_trace(go.Scatter(x=xr.index, y=xr["flux"], mode="lines", name="X-ray flux",
                             line=dict(width=1, color="#1f77b4")), row=1, col=1)
    for lvl, name in ((1e-6, "C1"), (1e-5, "M1"), (1e-4, "X1")):
        fig.add_trace(go.Scatter(x=[t_start, now], y=[lvl, lvl], mode="lines+text", text=["", name],
                                 textposition="top left", line=dict(dash="dot", width=1, color="#999"),
                                 showlegend=False, hoverinfo="skip"), row=1, col=1)
    if len(flares):
        fig.add_trace(go.Scatter(
            x=flares["peak_time"], y=flares["peak_flux"], mode="markers+text",
            text=[c if f >= 5e-6 else "" for c, f in zip(flares["cls"], flares["peak_flux"])], hovertext=flares["cls"],
            textposition="top center", name="Detected flares",
            marker=dict(size=9, color=[FLARE_COLORS[c[0]] for c in flares["cls"]],
                        line=dict(width=1, color="#333"))), row=1, col=1)
    fig.update_yaxes(type="log", range=[-8, -3], row=1, col=1)
    try:                                                    # shade night over Nairobi: no HF blackout there
        for a, b in night_spans(t_start, now + timedelta(hours=1), *HF_SITES["Nairobi"]):
            fig.add_vrect(x0=a.isoformat(), x1=b.isoformat(), fillcolor="rgba(100,110,160,0.14)",
                          line_width=0, layer="below", row=1, col=1)
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name="Night in Nairobi (no HF blackout)",
                                 marker=dict(size=11, symbol="square", color="rgba(100,110,160,0.45)")),
                      row=1, col=1)
    except Exception as exc:
        log(f"night shading skipped: {exc}")

    # 2, 3 - solar wind
    if not sw.empty:
        fig.add_trace(go.Scatter(x=sw.index, y=sw["speed"], mode="lines", name="Solar wind speed",
                                 line=dict(width=1, color="#1e8449")), row=2, col=1)
        if shocks:
            fig.add_trace(go.Scatter(x=shocks, y=[sw["speed"].asof(t) for t in shocks], mode="markers",
                                     name="Shock candidate", marker=dict(symbol="triangle-down", size=12,
                                                                         color="#c0392b")), row=2, col=1)
        fig.add_trace(go.Scatter(x=sw.index, y=sw["bz_gsm"], mode="lines", name="Bz",
                                 line=dict(width=1, color="#8e44ad")), row=3, col=1)
    # 4 - Dst
    if len(dst):
        fig.add_trace(go.Scatter(x=dst.index, y=dst.values, mode="lines", name="Dst",
                                 line=dict(width=1.5, color="#c0392b")), row=4, col=1)
        fig.add_trace(go.Scatter(x=[t_start, now], y=[args.dst, args.dst], mode="lines+text",
                                 text=["", f"Moderate storm ({args.dst:g} nT)"], textposition="top left",
                                 showlegend=False, hoverinfo="skip",
                                 line=dict(dash="dot", color="#c0392b", width=1)), row=4, col=1)
    # 5 - Kp as 3-hourly bars, coloured by storm level
    if len(kp):
        kb = kp.resample("3h").max().dropna()
        fig.add_trace(go.Bar(x=kb.index + timedelta(hours=1.5), y=kb.values, width=3 * 3600 * 1000 * 0.85,
                             name="Kp (3-hourly)", marker=dict(color=[kp_color(k) for k in kb.values]),
                             hovertemplate="Kp %{y:.1f}<extra></extra>"), row=5, col=1)
        fig.add_trace(go.Scatter(x=[t_start, now], y=[args.kp, args.kp], mode="lines+text",
                                 text=["", f"Storm threshold (Kp {args.kp:g})"], textposition="top left",
                                 showlegend=False, hoverinfo="skip",
                                 line=dict(dash="dot", color="#c0392b", width=1)), row=5, col=1)

    # axes: boxed panels, own tick labels and a y-axis label on every panel
    box = dict(showline=True, mirror=True, linecolor="#9aa5b1", gridcolor="#e8ebef", zeroline=False)
    fig.update_xaxes(showticklabels=True, tickformat="%b %d", range=[t_start, now + timedelta(hours=1)], **box)
    fig.update_yaxes(**box)
    for r, t in enumerate(["X-ray flux (W/m\u00b2)", "Speed (km/s)", "Bz (nT)", "Dst (nT)", "Kp index (0-9)"], 1):
        fig.update_yaxes(title_text=t, row=r, col=1)
    fig.update_yaxes(range=[0, 9.5], dtick=3, row=5, col=1)
    fig.update_xaxes(title_text=f"Date and time (UTC), {now.year}", row=5, col=1)
    fig.update_annotations(font=dict(size=13), x=0.0, xanchor="left")
    fig.update_layout(template="plotly_white", plot_bgcolor="#f7f9fb", height=1050, hovermode="x unified",
                      margin=dict(l=70, r=20, t=140, b=60),
                      legend=dict(orientation="h", xanchor="left", x=0, yanchor="bottom", y=1.06))
    return fig


CSS = """
body{background:#fff;color:#1f2933;font-family:system-ui,'Segoe UI',sans-serif;margin:0;padding:12px 20px}
.top{display:flex;justify-content:space-between;align-items:flex-start;gap:14px}
.top button{background:#1f6fb2;color:#fff;border:0;border-radius:6px;padding:7px 14px;font-size:13px;cursor:pointer;white-space:nowrap}
.top button:hover{background:#185a93}
h1{font-size:23px;margin:0 0 10px}
.banner{display:grid;grid-template-columns:repeat(7,minmax(0,1fr));gap:10px;margin-bottom:8px}
@media(max-width:1250px){.banner{grid-template-columns:repeat(4,minmax(0,1fr))}}
@media(max-width:700px){.banner{grid-template-columns:repeat(2,minmax(0,1fr))}}
.card{background:#f3f5f8;border:1px solid #e1e5ea;border-radius:8px;padding:8px 12px}
.big{font-size:22px;font-weight:700;line-height:1.2} .lbl{font-size:12px;color:#5f6b7a}
h2{margin:18px 0 4px;font-size:18px;border-top:2px solid #e1e5ea;padding-top:10px}
.colh{font-size:16px;font-weight:700;margin:0 0 4px}
.note{color:#5f6b7a;font-size:13px;margin:2px 0 6px}
.tbl{border-collapse:collapse;font-size:13px} .tbl th,.tbl td{padding:3px 9px;border-bottom:1px solid #e1e5ea;text-align:left}
.tbl th{color:#5f6b7a;background:#f7f8fa} .scroll{overflow-x:auto}
.alert{border:1px solid #e1e5ea;border-left:10px solid #1e8449;border-radius:8px;padding:8px 16px;margin:6px 0;background:#fbfcfd}
.alert p{margin:4px 0}
.two{display:grid;grid-template-columns:minmax(220px,1fr) 3fr;gap:22px;align-items:start;margin-top:6px}
.two.even{grid-template-columns:1fr 1.5fr}
@media(max-width:1000px){.two,.two.even{grid-template-columns:1fr}}
.chips div{background:#f3f5f8;border:1px solid #e1e5ea;border-radius:8px;padding:4px 12px;margin:0 0 6px;font-size:13px}
.hl td{background:#fff3cd;font-weight:700}
.warn{display:none;background:#fdecea;border:1px solid #f5c6cb;color:#922b21;padding:8px 14px;border-radius:6px;margin-bottom:10px}
.guide{columns:2;column-gap:28px;margin:4px 0 0;padding-left:18px} .guide li{margin:3px 0;font-size:13px;color:#3e4c59;break-inside:avoid}
.site{font-weight:700;color:#1f6fb2}
details{margin:6px 0;font-size:13px;color:#5f6b7a} summary{cursor:pointer}
@page{size:A4 landscape;margin:6mm}
@media print{
 body{padding:0;-webkit-print-color-adjust:exact;print-color-adjust:exact}
 .noprint,details{display:none!important}
 h1{font-size:20px;margin:0 0 6px}
 .banner{grid-template-columns:repeat(7,minmax(0,1fr))} .big{font-size:17px}
 .two{grid-template-columns:minmax(220px,1fr) 3fr} .two.even{grid-template-columns:1fr 1.5fr}
 .tbl{font-size:12px} .tbl th,.tbl td{padding:2px 8px}
 .guide li{font-size:12px;margin:2px 0}
 h2{margin:10px 0 4px;padding-top:6px}
 .pg2,.pg3{break-before:page}
 .pg2{max-height:735px;overflow:hidden}
 .nasa .tbl tr:nth-child(n+7){display:none}
 .pg2 h2{border-top:none;margin-top:0;padding-top:0}
 h2,.colh{break-after:avoid} .alert,.card,.tbl tr,.plotly-graph-div{break-inside:avoid}
}
"""

JS = """
const GEN = new Date("__GEN__"), LAST = new Date("__LAST__"), STALE = __STALE__;
const pad = n => String(n).padStart(2, "0");
function ago(ms){ const m = Math.round(ms / 60000); return m < 1 ? "just now" : m < 90 ? m + " min ago" : (m / 60).toFixed(1) + " h ago"; }
function tick(){
  const n = new Date(), e = new Date(n.getTime() + 3 * 3600 * 1000);
  document.getElementById("clock").textContent = pad(e.getUTCHours()) + ":" + pad(e.getUTCMinutes()) + ":" + pad(e.getUTCSeconds());
  document.getElementById("dataage").textContent = ago(n - LAST);
  document.getElementById("refreshed").textContent = ago(n - GEN);
  document.getElementById("stale").style.display = (n - GEN > STALE) ? "block" : "none";
}
tick(); setInterval(tick, 1000);

// Printing / saving as PDF: Plotly draws at screen width, so fit the charts to the paper first.
const PRINT_W = 1050, PRINT_H = 690;
function fitForPrint(){
  Array.from(document.querySelectorAll(".plotly-graph-div")).forEach(g => {
    try { Plotly.relayout(g, {width: PRINT_W, height: PRINT_H, "margin.t": 120, "margin.b": 50}); } catch (e) {}
    g.style.width = PRINT_W + "px"; g.style.height = PRINT_H + "px";   // shrink the host box too (prevents a blank page)
  });
}
function savePdf(){ fitForPrint(); setTimeout(() => window.print(), 500); }
window.addEventListener("beforeprint", fitForPrint);
window.addEventListener("afterprint", () => setTimeout(() => location.reload(), 300));
"""


def _pretty(df: pd.DataFrame) -> pd.DataFrame:
    ren = {}
    for c in df.columns:
        if c in PRETTY_COLS:
            ren[c] = PRETTY_COLS[c]
        elif str(c).startswith("peak_to_Kp>="):
            ren[c] = f"Flare to Kp {str(c)[len('peak_to_Kp>='):-2]} or more (hours)"
    return df.rename(columns=ren)


def _table_only(df, empty_msg="No events in this period.") -> str:
    if df is None or df.empty:
        return f"<p class='note'>{empty_msg}</p>"
    d = _pretty(df.copy())
    for c in d.columns:
        if str(d[c].dtype).startswith("datetime"):
            d[c] = d[c].dt.strftime("%d %b %H:%M")
    return d.to_html(index=False, border=0, na_rep="\u2013", classes="tbl", float_format="%.1f")


def _donki_view(d: pd.DataFrame, n: int = 8) -> pd.DataFrame:
    if d is None or d.empty:
        return d
    v = d.head(n).copy()
    for c in ("flare_to_cme_h", "flare_to_gst_h", "gst_max_kp"):
        if c not in v:
            v[c] = np.nan
    v["storm"] = np.where(v["gst_id"].notna(), "Storm linked", "No storm linked")
    return v[["flare_peak_utc", "class", "storm", "flare_to_cme_h", "flare_to_gst_h", "gst_max_kp"]]


def _live_view(h: pd.DataFrame, kp_thr: float) -> pd.DataFrame:
    if h is None or h.empty:
        return h
    v = h.copy()
    v["status"] = v["status"].replace({"complete": "Complete", "window open": "Still monitoring"})
    cols = ["flare_peak_utc", "class", "peak_to_shock_h", f"peak_to_Kp>={kp_thr:g}_h", "dst_min_nT",
            "peak_to_dst_min_h", "competing_flares", "status"]
    return v[[c for c in cols if c in v.columns]]


def write_dashboard(xr, flares, sw, kp, dst, shocks, heur, donki, hf, args, out_dir: Path) -> None:
    fig = build_figure(xr, flares, sw, kp, dst, shocks, args)
    now = pd.Timestamp.now(tz="UTC")
    last_t, last_flux = xr.index[-1], xr["flux"].iloc[-1]
    sat = (f"GOES-{int(xr['satellite'].iloc[-1])}" if "satellite" in xr and pd.notna(xr["satellite"].iloc[-1])
           else "GOES (primary)")
    idx, lvl = alert_level(last_flux)
    kn = kenya_now(last_flux)
    n_lit = sum(sun_elevation(now, *ll) > 0 for ll in HF_SITES.values())

    if kn["state"] == "Night":
        kn_txt, kn_sub = "Unaffected", "Night: flares cannot disrupt HF here"
        kn_line = "Kenya is in darkness, so HF radio is not affected by flares right now (radio blackouts only occur on the sunlit side)."
    elif idx <= 1:
        kn_txt, kn_sub = "No disturbance", "Normal daytime conditions"
        kn_line = "Kenya is in daylight, and there is no flare disturbance to HF radio at the moment."
    else:
        kn_txt, kn_sub = f"Degraded to ~{kn['haf']:.0f} MHz", f"{lvl[1]} level"
        kn_line = f"Kenya is in daylight: HF radio may be degraded up to about {kn['haf']:.0f} MHz (Nairobi estimate)."
    kn_line += f" {n_lit} of {len(HF_SITES)} monitored East African sites are in daylight."

    kp_txt = f"{kp.iloc[-1]:.1f} &middot; {kp_label(kp.iloc[-1])}" if len(kp) else "n/a"
    dst_txt = f"{dst.iloc[-1]:.0f} nT &middot; {dst_label(dst.iloc[-1])}" if len(dst) else "n/a"

    win = {"6-hour": "6 hours", "1-day": "24 hours", "3-day": "3 days", "7-day": "7 days"}[args.window]
    if len(flares):
        cnt = flares["cls"].str[0].value_counts()
        big = flares.loc[flares["peak_flux"].idxmax()]
        chips = "".join(f"<div><b>{L}-class flares:</b> {int(cnt.get(L, 0))}</div>" for L in ("X", "M", "C"))
        chips += (f"<div><b>Largest:</b> {big['cls']} at {big['peak_time']:%d %b %H:%M} UTC</div>"
                  f"<div><b>Latest:</b> {flares.iloc[-1]['cls']} at {flares.iloc[-1]['peak_time']:%d %b %H:%M} UTC</div>")
    else:
        chips = "<div>No flares of C1 or above in this period</div>"

    level_rows = "".join(
        f'<tr class="{"hl" if i == idx else ""}"><td>{cls}{" &larr; now" if i == idx else ""}</td>'
        f'<td><span style="color:{lv[4]}">&#9679;</span> {lv[1]}</td><td>{lv[3]}</td></tr>'
        for i, (lv, cls) in enumerate(zip(LEVELS, LEVEL_CLASS_TXT)))
    levels_tbl = ('<table class="tbl"><tr><th>Flare class (peak X-ray strength)</th><th>Alert level</th>'
                  '<th>What it means for HF radio (sunlit side only)</th></tr>' + level_rows + "</table>")

    fig_html = fig.to_html(full_html=False, include_plotlyjs="cdn", config={"displaylogo": False, "responsive": True})
    stale_ms = (args.interval * 3 + 120) * 1000
    js = (JS.replace("__GEN__", now.strftime("%Y-%m-%dT%H:%M:%SZ"))
            .replace("__LAST__", last_t.strftime("%Y-%m-%dT%H:%M:%SZ")).replace("__STALE__", str(stale_ms)))
    title = "Space Weather Monitor: Solar Flares, Geomagnetic Activity and HF Radio Impact over Kenya / East Africa"
    q = "How long after a flare did a geomagnetic storm follow?"

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<meta http-equiv="refresh" content="{args.interval}">
<title>{title}</title><style>{CSS}</style></head><body>
<div class="top"><h1>{title}</h1><button class="noprint" onclick="savePdf()">Save as PDF</button></div>
<div class="warn" id="stale">&#9888; This page has stopped updating. The script must stay running to keep it live: restart it to resume.</div>

<div class="banner">
 <div class="card"><div class="lbl">Sun's flare strength now</div><div class="big">{flare_class(last_flux)}</div><div class="lbl">{lvl[1]}</div></div>
 <div class="card"><div class="lbl">Geomagnetic activity now (Kp)</div><div class="big">{kp_txt}</div><div class="lbl">5 or more = storm</div></div>
 <div class="card"><div class="lbl">Geomagnetic storm index now (Dst)</div><div class="big">{dst_txt}</div><div class="lbl">Below -50 nT = storm</div></div>
 <div class="card"><div class="lbl">Nairobi sun now</div><div class="big">{kn['elev']:.0f}&deg; {kn['state']}</div><div class="lbl">{kn['eat']} EAT</div></div>
 <div class="card"><div class="lbl">HF radio at Nairobi now (estimate)</div><div class="big">{kn_txt}</div><div class="lbl">{kn_sub}</div></div>
 <div class="card"><div class="lbl">Latest satellite reading ({sat})</div><div class="big" id="dataage">...</div><div class="lbl">Solar X-ray data</div></div>
 <div class="card"><div class="lbl">Time now (EAT)</div><div class="big" id="clock">--:--:--</div><div class="lbl">Dashboard refreshed <span id="refreshed">...</span></div></div>
</div>

<section class="pg1">
<h2>Real-time alert status</h2>
<div class="alert" style="border-left-color:{lvl[4]}">
 <div class="lbl">CURRENT FLARE ALERT LEVEL (from the latest X-ray reading)</div>
 <div class="big" style="color:{lvl[4]}">{lvl[2]}</div>
 <p>{lvl[3]}</p>
 <p><b>Kenya / East Africa now:</b> {kn_line}</p>
</div>
<div class="two">
 <div class="chips"><p class="note"><b>Flares in the last {win}</b></p>{chips}</div>
 <div><p class="note"><b>Alert levels by flare class</b> (NOAA radio blackout scale). The highlighted row is the current level.</p>{levels_tbl}</div>
</div>
</section>

<section class="pg2">
<h2>Space weather over the last {win}</h2>
{fig_html}
</section>

<section class="pg3">
<p class="note"><b>How to read the charts</b></p>
<ul class="guide">
 <li><b>X-ray flux:</b> brightness of the Sun in X-rays. Flare classes A, B, C, M, X are each 10 times stronger than the last. M and X flares can black out HF radio.</li>
 <li><b>Solar wind speed:</b> measured about 1.5 million km sunward of Earth. A sudden jump can signal a CME shock arriving, with roughly 15 to 60 minutes of warning.</li>
 <li><b>Bz:</b> north-south direction of the solar wind's magnetic field. Sustained negative (southward) values let solar wind energy into Earth's magnetic field.</li>
 <li><b>Dst:</b> hourly storm-strength index at Earth. Values below about -50 nT mean a moderate storm.</li>
 <li><b>Kp:</b> global geomagnetic activity every 3 hours, from 0 to 9. Kp 5 or more is a storm (G1 to G5).</li>
</ul>

<h2>HF radio blackout impact over Kenya / East Africa at flare peak</h2>
<p class="note">Flares {args.hf_class} or above (the five largest if none). Values are model-based estimates (NOAA D-RAP method). Sites: <span class="site">Nairobi, Mombasa, Kisumu, Lodwar</span> (Kenya); <span class="site">Kampala, Dar es Salaam, Addis Ababa</span>. Grey bands on the X-ray chart = night in Nairobi.</p>
<div class="scroll">{_table_only(hf)}</div>
<details><summary>How the HF impact is estimated</summary>
NOAA's D-RAP method: highest affected frequency = 10&middot;log10(flux) + 65 MHz where the sun is overhead, reduced by cos(zenith)^0.75 elsewhere
(1 dB absorption, vertical path; slanted paths are affected at higher frequencies). Radio blackouts only occur on the sunlit side. The upper atmosphere stays
sunlit a little after ground sunset, so just after sunset is a grey zone.</details>

<div class="two even">
 <div><div class="colh">{q} (NASA DONKI)</div>
  <p class="note">Delays linked by NASA analysts (flare, CME, storm); latest {min(len(donki) if donki is not None else 0, 8)} M-class flares of the last {min(args.donki_days, 60)} days. Most flares do not cause storms.</p>
  <div class="scroll nasa">{_table_only(_donki_view(donki), "NASA DONKI data not available right now.")}</div></div>
 <div><div class="colh">{q} (LIVE DATA)</div>
  <p class="note">Looks {args.min_h:g} to {args.max_h:g} hours after each {args.assoc_class}+ flare for a solar-wind shock, Kp {args.kp:g}+ and Dst {args.dst:g} nT or lower.</p>
  <div class="scroll">{_table_only(_live_view(heur, args.kp))}</div></div>
</div>
<details><summary>About these two tables</summary>
A dash means nothing was found (yet). Kp is reported every 3 hours, so Kp timing can be up to 3 hours early. "Other candidates" above 0 means several flares could
explain the same event, so the match is uncertain. Live shock detection only covers solar-wind data saved by this script (NOAA now serves about 24 hours),
so older flares show no shock there; use the NASA table for those.</details>
</section>
<script>{js}</script>
</body></html>"""
    tmp = out_dir / "dashboard.tmp"
    tmp.write_text(html, encoding="utf-8")
    try:
        tmp.replace(out_dir / "dashboard.html")            # atomic: browser never sees a half-written page
    except OSError:
        (out_dir / "dashboard.html").write_text(html, encoding="utf-8")


# --------------------------------------------------------------------------- #
# Main loop
# --------------------------------------------------------------------------- #
def run_cycle(args, out_dir: Path) -> None:
    xr = fetch_xray(args.window)
    if xr.empty:
        log("No X-ray data this cycle; will retry.")
        return False
    flares = detect_flares(xr, threshold=class_to_flux(args.detect_class))
    process_alerts(flares, args, out_dir)

    sw, kp, dst = fetch_solar_wind(out_dir / "solar_wind_cache.csv"), fetch_kp(), fetch_dst()
    shocks = detect_shocks(sw)
    heur = heuristic_latency(flares, shocks, kp, dst, args)
    hf = kenya_hf_table(flares, args)
    donki = donki_latency(args)

    flares.to_csv(out_dir / "flares_detected.csv", index=False)
    heur.to_csv(out_dir / "latency_heuristic.csv", index=False)
    hf.to_csv(out_dir / "hf_kenya_impact.csv", index=False)
    if not donki.empty:
        donki.to_csv(out_dir / "latency_donki.csv", index=False)
    write_dashboard(xr, flares, sw, kp, dst, shocks, heur, donki, hf, args, out_dir)
    log(f"OK: {len(flares)} flares >= {args.detect_class}, {len(shocks)} L1 shock candidates, "
        f"data age {(pd.Timestamp.now(tz='UTC') - xr.index[-1]).total_seconds() / 60:.1f} min "
        f"-> {out_dir / 'dashboard.html'}")
    return True


def parse_args():
    p = argparse.ArgumentParser(description="GOES flare alert dashboard + geomagnetic latency")
    p.add_argument("--window", default="7-day", choices=["6-hour", "1-day", "3-day", "7-day"],
                   help="X-ray history to pull (7-day recommended for latency analysis)")
    p.add_argument("--once", action="store_true", help="run a single pass and exit (default: keep running and refresh)")
    p.add_argument("--loop", action="store_true", help=argparse.SUPPRESS)   # kept for older commands; live mode is now the default
    p.add_argument("--no-browser", action="store_true", help="do not open the dashboard in the browser automatically")
    p.add_argument("--interval", type=int, default=60, help="poll / page-refresh seconds (>=60)")
    p.add_argument("--detect-class", default="C1.0", help="lowest class tracked as a flare")
    p.add_argument("--alert-class", default="M1.0", help="lowest class that triggers an alert")
    p.add_argument("--assoc-class", default="M1.0", help="lowest class used for impact-latency analysis")
    p.add_argument("--hf-class", default="M1.0", help="lowest class listed in the Kenya HF impact table")
    p.add_argument("--kp", type=float, default=5.0, help="Kp threshold (5 = G1 storm)")
    p.add_argument("--dst", type=float, default=-50.0, help="Dst threshold in nT")
    p.add_argument("--min-h", type=float, default=10.0, help="earliest plausible CME arrival (h)")
    p.add_argument("--max-h", type=float, default=120.0, help="latest CME arrival considered (h)")
    p.add_argument("--donki-days", type=int, default=30, help="DONKI look-back (days)")
    p.add_argument("--no-donki", action="store_true")
    p.add_argument("--webhook", default=None, help="Slack-style webhook (or set FLARE_WEBHOOK_URL)")
    p.add_argument("--email-to", default=os.environ.get("FLARE_EMAIL_TO"),
                   help="send flare alerts to this address (optional; or set FLARE_EMAIL_TO). Off by default")
    p.add_argument("--test-email", action="store_true", help="send a test email and exit")
    p.add_argument("--out-dir", default="flare_dashboard")
    a = p.parse_args()
    a.interval = max(a.interval, 60)
    return a


def main() -> None:
    try:
        import plotly  # noqa: F401
    except ImportError:
        print("Plotly is required for the dashboard. Install it with:  pip install plotly\n"
              "(Anaconda users: conda install plotly), then run the script again.")
        return
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.test_email:
        ok = send_email("Test: solar flare alert email", "If you can read this, email alerts are working.", args)
        log("Test email sent." if ok else "Test email NOT sent - see message above.")
        return
    page = (out_dir / "dashboard.html").resolve()
    opened = False
    if not args.once:
        log(f"Live mode: refreshing every {args.interval} s. Leave this running; stop it with Ctrl+C (or the Stop button in Spyder).")
    try:
        while True:
            try:
                ok = run_cycle(args, out_dir)
                if ok and not opened and not args.no_browser:
                    webbrowser.open(page.as_uri())
                    opened = True
            except Exception as exc:  # keep the monitor alive
                log(f"cycle error: {type(exc).__name__}: {exc}")
            if args.once:
                break
            time.sleep(args.interval)
    except KeyboardInterrupt:
        log("Stopped by user.")


if __name__ == "__main__":
    sys.exit(main())
