#!/usr/bin/env python3
"""
Spacecraft identification by Doppler curve matching against JPL Horizons.

Given a TDM file from an unknown recording, automatically queries all known
spacecraft from JPL Horizons, computes expected Doppler, and identifies the
best match by minimising RMS residual after removing a DC offset.

Centre frequency and bandwidth are read from the TDM file automatically.
The only input required from the user is the observer location.

Usage:
  python identify_spacecraft.py \\
      --tdm unknown.tdm \\
      --station 16.6752,52.3699,0.07

  python identify_spacecraft.py \\
      --tdm unknown.tdm \\
      --station 16.6752,52.3699,0.07 \\
      --plot
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
# Built-in candidate list
# Negative integers are spacecraft in JPL Horizons.
# Horizons returns "no data" for IDs without valid ephemerides — safe to include
# historical or uncertain entries; they are skipped automatically.
# ---------------------------------------------------------------------------

DEFAULT_SPACECRAFT = {
    # Lunar orbiters
    '-155': 'KPLO / Danuri (Korea, 2022–)',
    '-85':  'Chandrayaan-2 orbiter (India, 2019–)',
    '-152': 'Chandrayaan-3 propulsion module (India, 2023–)',
    '-12':  'LADEE (NASA, deorbited 2014)',
    '-167': 'LADEE (alt ID)',
    '-177': 'GRAIL-A / Ebb (NASA, 2012)',
    '-181': 'GRAIL-B / Flow (NASA, 2012)',
    '-191': 'LCROSS (NASA, 2009)',
    '-198': 'CAPSTONE (NASA, 2022–)',
    # Mars / deep space (S-band)
    '-74':  'MRO — Mars Reconnaissance Orbiter (NASA, 2006–)',
    '-116': 'Mars Odyssey (NASA, 2001–)',
    '-130': 'MAVEN (NASA, 2014–)',
    '-98':  'New Horizons (NASA, 2006–)',
    # Other
    '-18':  'WIND (NASA, 1994–)',
    '-55':  'Ulysses (ESA/NASA, ended 2009)',
    '-68':  'Dawn (NASA, ended 2018)',
}


# ---------------------------------------------------------------------------
# Parse TDM
# ---------------------------------------------------------------------------

def _parse_doy_timestamp(s):
    s = re.sub(r'T(\d{2}:\d{2}:\d{2}):(\d+)', r'T\1.\2', s.strip())
    mp = re.match(r'(\d{4})-(\d{3})T(\d{2}):(\d{2}):(\d{2})(?:[.,](\d+))?', s)
    if not mp:
        return None
    year, doy, hh, mm, ss, frac = mp.groups()
    base = datetime(int(year), 1, 1, tzinfo=timezone.utc) + \
           timedelta(days=int(doy) - 1)
    return base.replace(hour=int(hh), minute=int(mm), second=0) + \
           timedelta(seconds=int(ss) + (float('0.' + frac) if frac else 0.0))


def parse_tdm(path):
    """
    Parse a CCSDS TDM file.
    Returns times, freqs, freq_offset [Hz], sample_rate [Hz].
    """
    times, freqs = [], []
    freq_offset = None
    sample_rate = None

    with open(path) as f:
        for line in f:
            line = line.strip()
            m = re.match(r'FREQ_OFFSET\s*=\s*([\d.eE+-]+)', line)
            if m:
                freq_offset = float(m.group(1))
            m = re.search(r'_(\d+)_fc\.sigmf', line)   # iq_to_tdm filename convention
            if m:
                sample_rate = float(m.group(1))
            m = re.match(r'RECEIVE_FREQ_2\s*=\s*(\S+)\s+([\d.eE+-]+)', line)
            if m:
                t = _parse_doy_timestamp(m.group(1))
                if t:
                    times.append(t)
                    freqs.append(float(m.group(2)))

    return times, freqs, freq_offset, sample_rate


# ---------------------------------------------------------------------------
# JPL Horizons query
# ---------------------------------------------------------------------------

def query_horizons(spk_id, site_coord, t_start, t_stop):
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
        'STEP_SIZE':   "'1m'",
        'QUANTITIES':  "'4,20'",
        'CAL_FORMAT':  "'CAL'",
        'TIME_DIGITS': "'MINUTES'",
    }
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read().decode())
    except Exception:
        return None

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
        m = re.match(
            r'\s*(\d{4}-\w{3}-\d{2}\s+\d{2}:\d{2})'
            r'\s+\S+\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s+([-\d.]+)',
            line
        )
        if m:
            dt     = datetime.strptime(m.group(1), '%Y-%b-%d %H:%M').replace(tzinfo=timezone.utc)
            elev   = float(m.group(3))
            deldot = float(m.group(4))
            rows.append((dt, deldot, elev))
    return rows if rows else None


# ---------------------------------------------------------------------------
# Match one candidate
# ---------------------------------------------------------------------------

def match_candidate(spk_id, name, active_tdm, center_freq, bandwidth_hz,
                    site_coord, t_start, t_stop):
    label = f"{spk_id}  {name}"
    print(f"  {label} ... ", end='', flush=True)

    hor = query_horizons(spk_id, site_coord, t_start, t_stop)
    if not hor:
        print("no ephemeris")
        return None

    max_elev = max(e for _, _, e in hor)
    if max_elev < 0:
        print(f"below horizon (max {max_elev:.1f}°)")
        return None

    print(f"elev {min(e for _,_,e in hor):.1f}°–{max_elev:.1f}°", end='  ')

    hor_dop = [(t, -deldot * center_freq / C_KMS) for t, deldot, _ in hor]

    pairs = []
    for t_tdm, f_tdm in active_tdm:
        best_dt, best_hor = None, None
        for t_h, d_h in hor_dop:
            dt = abs((t_tdm - t_h).total_seconds())
            if dt <= 90 and (best_dt is None or dt < best_dt):
                best_dt = dt
                best_hor = d_h
        if best_hor is not None:
            pairs.append((t_tdm, f_tdm, best_hor, f_tdm - best_hor))

    if not pairs:
        print("no time overlap")
        return None

    diffs     = [p[3] for p in pairs]
    dc_offset = sum(diffs) / len(diffs)
    residuals = [d - dc_offset for d in diffs]
    rms       = math.sqrt(sum(r**2 for r in residuals) / len(residuals))

    if bandwidth_hz and abs(dc_offset) > bandwidth_hz / 2:
        print(f"out of band (DC={dc_offset:+.0f} Hz, bw=±{bandwidth_hz/2:.0f} Hz)")
        return None

    print(f"DC={dc_offset:+.0f} Hz  RMS={rms:.1f} Hz")

    return {
        'id':         spk_id,
        'name':       name,
        'n_pairs':    len(pairs),
        'dc_offset':  dc_offset,
        'rms':        rms,
        'pairs':      pairs,
        'hor_dop':    hor_dop,
        'elev_range': (min(e for _,_,e in hor), max_elev),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tdm',     required=True,
                    help='TDM file to identify')
    ap.add_argument('--station', required=True,
                    help='Observer location: lon,lat,alt_km  (WGS84)')
    ap.add_argument('--plot',    action='store_true',
                    help='Save Doppler comparison plot to spacecraft_identification.png')
    args = ap.parse_args()

    print(f"\nLoading TDM: {args.tdm}")
    times, freqs, freq_offset, sample_rate = parse_tdm(args.tdm)

    if freq_offset is None:
        print("ERROR: FREQ_OFFSET not found in TDM.")
        sys.exit(1)

    active = [(t, f) for t, f in zip(times, freqs) if abs(f) > 10]
    print(f"  Centre frequency       : {freq_offset/1e6:.6f} MHz")
    if sample_rate:
        print(f"  Bandwidth              : ±{sample_rate/2:.0f} Hz")
    print(f"  Total measurements     : {len(times)}")
    print(f"  Active (|offset|>10 Hz): {len(active)}")
    if not active:
        print("ERROR: no active measurements found.")
        sys.exit(1)

    t_start_str = active[0][0].strftime('%Y-%m-%d %H:%M')
    t_stop_str  = active[-1][0].strftime('%Y-%m-%d %H:%M')
    print(f"  Active window          : {t_start_str} – {t_stop_str} UTC")

    print(f"\nTesting {len(DEFAULT_SPACECRAFT)} known spacecraft:")
    results = []
    for cid, name in DEFAULT_SPACECRAFT.items():
        r = match_candidate(cid, name, active, freq_offset, sample_rate,
                            args.station, t_start_str, t_stop_str)
        if r:
            results.append(r)

    if not results:
        print("\nNo match found. The spacecraft may not be in the built-in list.")
        sys.exit(1)

    results.sort(key=lambda x: x['rms'])
    best = results[0]

    print(f"\n{'='*60}")
    print(f"  BEST MATCH: {best['id']}  {best['name']}")
    print(f"  Elevation              : "
          f"{best['elev_range'][0]:.1f}° – {best['elev_range'][1]:.1f}°")
    print(f"  DC offset (TX – SDR)   : {best['dc_offset']:+.1f} Hz")
    print(f"  Transmit frequency     : "
          f"{(freq_offset + best['dc_offset'])/1e6:.6f} MHz")
    print(f"  RMS residual           : {best['rms']:.1f} Hz")
    print(f"  Matched points         : {best['n_pairs']}")
    print(f"{'='*60}")

    if len(results) > 1:
        print(f"\nAll matched candidates ranked by RMS:")
        for r in results:
            marker = "  <-- best" if r is best else ""
            print(f"  {r['id']:6s}  {r['name'][:42]:42s}  "
                  f"RMS={r['rms']:7.1f} Hz{marker}")

    if args.plot:
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates

        fig, ax = plt.subplots(figsize=(13, 5))
        fig.suptitle(
            f"Spacecraft identification — {best['id']}  {best['name']}\n"
            f"RMS: {best['rms']:.1f} Hz  |  "
            f"DC: {best['dc_offset']:+.1f} Hz  |  "
            f"TX freq: {(freq_offset + best['dc_offset'])/1e6:.6f} MHz",
            fontsize=11, fontweight='bold')

        hor_t = [t for t, _ in best['hor_dop']]
        hor_f = [d + best['dc_offset'] for _, d in best['hor_dop']]
        ax.plot(hor_t, hor_f, color='steelblue', linewidth=2.0,
                label=f"JPL Horizons {best['id']} (DC-shifted)")
        ax.scatter([t for t, f in active], [f for _, f in active],
                   s=4, color='darkorange', alpha=0.6, label='TDM (measured)')

        ax.set_xlabel('UTC')
        ax.set_ylabel(f"Doppler offset from {freq_offset/1e6:.4f} MHz [Hz]")
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0, 60, 10)))
        ax.tick_params(axis='x', rotation=25)
        ax.grid(True, alpha=0.35)
        ax.legend(fontsize=10)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:+.0f}'))

        plt.tight_layout()
        out = 'spacecraft_identification.png'
        plt.savefig(out, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved: {out}")


if __name__ == '__main__':
    main()
