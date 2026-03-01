#!/usr/bin/env python3
"""
Spacecraft identification by Doppler curve matching against JPL Horizons.

Given an unknown TDM file and a list of candidate spacecraft (JPL Horizons IDs),
queries each candidate's range-rate, converts to expected Doppler, and finds
the best match by minimising RMS residual after removing a DC offset.

The DC offset represents the difference between the spacecraft's actual transmit
frequency and the SDR centre frequency — typically tens of kHz for lunar missions.
The RMS residual (after removing DC) is the true quality metric: a good match
gives RMS < 200 Hz; a wrong candidate gives RMS in the thousands.

Usage:
  python identify_spacecraft.py --tdm SP5LOT_20260221_140504.tdm \\
      --freq 2260790300 --station 16.6752,52.3699,0.07 \\
      --candidates -155 -85 -152

Known JPL Horizons spacecraft IDs (S-band lunar missions):
  -155  KPLO / Danuri  (Korea, 2022–, ~2260.8 MHz)
  -85   Chandrayaan-2 orbiter  (India, 2019–)
  -152  Chandrayaan-2 (alternate)
  -98   New Horizons  (Pluto flyby, deep space)
  (add more as needed — negative integers are spacecraft in Horizons)
"""

import argparse
import json
import math
import re
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone, timedelta


C_KMS = 299_792.458   # km/s


# ---------------------------------------------------------------------------
# Parse TDM
# ---------------------------------------------------------------------------

def parse_tdm(path):
    """Return list of (datetime UTC, freq_offset_hz) from RECEIVE_FREQ_2 lines."""
    times, freqs = [], []
    with open(path) as f:
        for line in f:
            m = re.match(r'RECEIVE_FREQ_2\s*=\s*(\S+)\s+([\d.eE+-]+)', line.strip())
            if not m:
                continue
            ts = re.sub(r'T(\d{2}:\d{2}:\d{2}):(\d+)', r'T\1.\2', m.group(1))
            mp = re.match(r'(\d{4})-(\d{3})T(\d{2}):(\d{2}):(\d{2})(?:[.,](\d+))?', ts)
            if not mp:
                continue
            year, doy, hh, mm, ss, frac = mp.groups()
            base = datetime(int(year), 1, 1, tzinfo=timezone.utc) + \
                   timedelta(days=int(doy) - 1)
            t = base.replace(hour=int(hh), minute=int(mm), second=0) + \
                timedelta(seconds=int(ss) + (float('0.' + frac) if frac else 0.0))
            times.append(t)
            freqs.append(float(m.group(2)))
    return times, freqs


# ---------------------------------------------------------------------------
# JPL Horizons query
# ---------------------------------------------------------------------------

def query_horizons(spk_id, site_coord, t_start, t_stop, step='1m'):
    """
    Query range-rate (deldot, km/s) and elevation for a spacecraft from JPL Horizons.
    Returns list of (datetime UTC, deldot_km_s, elevation_deg).
    """
    lon, lat, alt = site_coord.split(',')
    params = {
        'format':      'json',
        'COMMAND':     f"'{spk_id}'",
        'OBJ_DATA':    "'NO'",
        'MAKE_EPHEM':  "'YES'",
        'TABLE_TYPE':  "'OBSERVER'",
        'CENTER':      "'coord@399'",
        'COORD_TYPE':  "'GEODETIC'",
        'SITE_COORD':  f"'{lon},{lat},{alt}'",
        'START_TIME':  f"'{t_start}'",
        'STOP_TIME':   f"'{t_stop}'",
        'STEP_SIZE':   f"'{step}'",
        'QUANTITIES':  "'4,20'",
        'CAL_FORMAT':  "'CAL'",
        'TIME_DIGITS': "'MINUTES'",
    }
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=30) as r:
        data = json.loads(r.read().decode())

    result = data.get('result', '')
    if 'No ephemeris' in result or 'Cannot' in result or 'ERROR' in result[:200]:
        return None

    rows = []
    in_data = False
    for line in result.splitlines():
        if '$$SOE' in line:
            in_data = True
            continue
        if '$$EOE' in line:
            break
        if not in_data:
            continue
        # Format: date  flag  azimuth  elevation  range_AU  deldot_km_s
        m = re.match(
            r'\s*(\d{4}-\w{3}-\d{2}\s+\d{2}:\d{2})'
            r'\s+\S+\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s+([-\d.]+)',
            line
        )
        if m:
            dt = datetime.strptime(m.group(1), '%Y-%b-%d %H:%M').replace(tzinfo=timezone.utc)
            elev   = float(m.group(3))   # degrees above horizon
            deldot = float(m.group(4))   # km/s, positive = receding
            rows.append((dt, deldot, elev))
    return rows if rows else None


# ---------------------------------------------------------------------------
# Match one candidate
# ---------------------------------------------------------------------------

def match_candidate(spk_id, active_tdm, center_freq, site_coord, t_start, t_stop):
    """
    Query Horizons for spk_id, compute expected Doppler, compare with TDM.
    Returns dict with match statistics, or None if query fails.
    """
    print(f"  Querying Horizons for ID {spk_id} ...", end=' ', flush=True)
    hor = query_horizons(spk_id, site_coord, t_start, t_stop)
    if not hor:
        print("no data")
        return None
    print(f"{len(hor)} rows  (elev {hor[0][2]:.1f}°–{hor[-1][2]:.1f}°)")

    # Convert deldot -> Doppler offset
    # deldot > 0 (receding) -> received freq < center -> negative offset
    hor_dop = [(t, -deldot * center_freq / C_KMS, elev)
               for t, deldot, elev in hor]

    # Match TDM points to nearest Horizons point (within 90 s)
    pairs = []
    for t_tdm, f_tdm in active_tdm:
        best_dt, best_hor = None, None
        for t_h, d_h, _ in hor_dop:
            dt = abs((t_tdm - t_h).total_seconds())
            if dt <= 90 and (best_dt is None or dt < best_dt):
                best_dt = dt
                best_hor = d_h
        if best_hor is not None:
            pairs.append((t_tdm, f_tdm, best_hor, f_tdm - best_hor))

    if not pairs:
        print(f"    -> no time matches found")
        return None

    diffs = [p[3] for p in pairs]
    dc_offset = sum(diffs) / len(diffs)
    residuals = [d - dc_offset for d in diffs]
    rms = math.sqrt(sum(r**2 for r in residuals) / len(residuals))

    return {
        'id':        spk_id,
        'n_pairs':   len(pairs),
        'dc_offset': dc_offset,
        'rms':       rms,
        'pairs':     pairs,
        'hor_dop':   [(t, d) for t, d, _ in hor_dop],
        'elev_range': (hor[0][2], hor[-1][2]),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tdm',        required=True, help='TDM file to identify')
    ap.add_argument('--freq',       type=float, required=True, help='SDR centre frequency [Hz]')
    ap.add_argument('--station',    required=True, help='lon,lat,alt_km  (WGS84)')
    ap.add_argument('--candidates', nargs='+', required=True,
                    help='JPL Horizons spacecraft IDs to test (e.g. -155 -85)')
    ap.add_argument('--plot',       action='store_true', help='Save comparison plot')
    args = ap.parse_args()

    # Load TDM — keep only active (non-zero) measurements
    print(f"\nLoading TDM: {args.tdm}")
    times, freqs = parse_tdm(args.tdm)
    active = [(t, f) for t, f in zip(times, freqs) if abs(f) > 10]
    print(f"  Total measurements : {len(times)}")
    print(f"  Active (|offset|>10 Hz) : {len(active)}")
    if not active:
        print("ERROR: no active measurements found.")
        sys.exit(1)

    t_start_str = active[0][0].strftime('%Y-%m-%d %H:%M')
    t_stop_str  = active[-1][0].strftime('%Y-%m-%d %H:%M')
    print(f"  Active window : {t_start_str} – {t_stop_str} UTC")

    # Query each candidate
    print(f"\nTesting {len(args.candidates)} candidate(s):")
    results = []
    for cand in args.candidates:
        r = match_candidate(cand, active, args.freq, args.station,
                            t_start_str, t_stop_str)
        if r:
            results.append(r)
            print(f"    DC offset = {r['dc_offset']:+.1f} Hz  |  "
                  f"RMS residual = {r['rms']:.1f} Hz  |  n = {r['n_pairs']}")

    if not results:
        print("\nNo candidates matched.")
        sys.exit(1)

    # Rank by RMS residual
    results.sort(key=lambda x: x['rms'])
    best = results[0]

    print(f"\n{'='*60}")
    print(f"  BEST MATCH: ID {best['id']}")
    print(f"  Elevation range : {best['elev_range'][0]:.1f}° – {best['elev_range'][1]:.1f}°")
    print(f"  DC offset (transmit freq – SDR centre) : {best['dc_offset']:+.1f} Hz")
    print(f"  Estimated transmit frequency : "
          f"{(args.freq + best['dc_offset'])/1e6:.6f} MHz")
    print(f"  RMS residual after DC removal : {best['rms']:.1f} Hz")
    print(f"  Matched pairs : {best['n_pairs']}")
    print(f"{'='*60}")

    if len(results) > 1:
        print(f"\nAll candidates ranked by RMS:")
        for r in results:
            marker = " <-- best" if r is best else ""
            print(f"  ID {r['id']:6s}  RMS={r['rms']:7.1f} Hz  DC={r['dc_offset']:+.1f} Hz{marker}")

    # Optional plot
    if args.plot:
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates

        fig, ax = plt.subplots(figsize=(13, 5))
        fig.suptitle(f'Spacecraft identification — best match: ID {best["id"]}\n'
                     f'RMS residual: {best["rms"]:.1f} Hz  |  '
                     f'DC offset: {best["dc_offset"]:+.1f} Hz  |  '
                     f'Est. TX freq: {(args.freq + best["dc_offset"])/1e6:.6f} MHz',
                     fontsize=12, fontweight='bold')

        # Horizons prediction (shifted by DC)
        hor_t = [t for t, _ in best['hor_dop']]
        hor_f = [d + best['dc_offset'] for _, d in best['hor_dop']]
        ax.plot(hor_t, hor_f, color='steelblue', linewidth=2.0,
                label=f'JPL Horizons ID {best["id"]} (DC-shifted)')
        # TDM active measurements
        ax.scatter([t for t, f in active], [f for _, f in active],
                   s=4, color='darkorange', alpha=0.6, label='TDM (measured)')

        ax.set_xlabel('UTC')
        ax.set_ylabel(f'Doppler offset from {args.freq/1e6:.4f} MHz [Hz]')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0, 60, 10)))
        ax.tick_params(axis='x', rotation=25)
        ax.grid(True, alpha=0.35)
        ax.legend(fontsize=10)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:+.0f}'))

        ax.text(0.02, 0.05,
                f'Best match: ID {best["id"]}  |  RMS: {best["rms"]:.1f} Hz  |  '
                f'DC: {best["dc_offset"]:+.1f} Hz  |  n={best["n_pairs"]}',
                transform=ax.transAxes, fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.tight_layout()
        out = 'spacecraft_identification.png'
        plt.savefig(out, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved: {out}")


if __name__ == '__main__':
    main()
