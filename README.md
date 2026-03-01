# spacecraft-doppler-id

Identify an unknown spacecraft by matching its measured Doppler curve against
JPL Horizons predictions for candidate spacecraft.

**Input:** a CCSDS TDM file with one-way Doppler measurements (e.g. from
[iq-to-tdm](https://github.com/SP5LOT/iq-to-tdm))
**Output:** ranked list of candidates by RMS residual + estimated transmit frequency

---

## How it works

1. Load `RECEIVE_FREQ_2` measurements from the TDM file
2. For each candidate spacecraft, query JPL Horizons for range-rate (`deldot`, km/s)
3. Convert range-rate → expected Doppler: `Δf = −deldot × f_centre / c`
4. Compute the mean difference between measured and predicted Doppler (DC offset)
5. Remove the DC offset and compute the RMS residual
6. Rank candidates by RMS — the lowest RMS is the best match

### DC offset vs RMS residual

| Metric | What it means |
|--------|--------------|
| **DC offset** | Difference between the spacecraft's actual transmit frequency and the SDR centre frequency. Typically tens of kHz for lunar missions. *Not* an error — it is the information. |
| **RMS residual** | How well the *shape* of the Doppler curve matches. **This is the identification metric.** |

A correct identification gives RMS of order **100 Hz**. A wrong candidate gives
RMS of thousands of Hz (the curve shapes simply do not match).

Example from KPLO/Danuri validation (2026-02-21, SP5LOT):

| Candidate | RMS | DC offset | Verdict |
|-----------|-----|-----------|---------|
| **-155 KPLO** | **120 Hz** | +32 536 Hz | ✓ match |
| -152 Chandrayaan-2 | 1 837 Hz | — | ✗ |
| -85 | 5 072 Hz | — | ✗ |

The DC offset of +32 kHz reveals that KPLO's actual downlink frequency is
approximately **2260.822 MHz**, not the 2260.790 MHz the SDR was tuned to.

---

## Requirements

```
Python 3.8+
matplotlib  (only for --plot)
```

Internet access to `ssd.jpl.nasa.gov` (JPL Horizons API).

---

## Usage

```bash
python identify_spacecraft.py \
    --tdm  SP5LOT_20260221_140504.tdm \
    --freq 2260790300 \
    --station 16.6752,52.3699,0.07 \
    --candidates -155 -85 -152
```

| Argument | Description |
|----------|-------------|
| `--tdm` | TDM file to identify (CCSDS KVN, `RECEIVE_FREQ_2`) |
| `--freq` | SDR centre frequency in Hz |
| `--station` | Observer coordinates: `lon,lat,alt_km` (WGS84) |
| `--candidates` | JPL Horizons spacecraft IDs to test (negative integers) |
| `--plot` | Save a Doppler comparison plot to `spacecraft_identification.png` |

### Station coordinates

The station coordinates must match the observer location in the TDM.
Use the same values as in [iq-to-tdm](https://github.com/SP5LOT/iq-to-tdm).

---

## Known JPL Horizons IDs

| ID | Spacecraft | Notes |
|----|-----------|-------|
| -155 | KPLO / Danuri | Korea, lunar orbit, 2022– |
| -85 | Chandrayaan-2 orbiter | India, lunar orbit, 2019– |
| -152 | Chandrayaan-2 (alt) | |
| -98 | New Horizons | Deep space |

Horizons uses negative integers for spacecraft. If you do not know the ID,
try a range of values — Horizons returns `no data` for unknown IDs.

The full list of spacecraft with ephemerides can be searched at
https://ssd.jpl.nasa.gov/horizons/

---

## Case study: identifying an unknown recording

On 2026-02-21, station SP5LOT recorded an IQ file starting at 14:05 UTC on
2260.790 MHz. The signal was active from 14:20 to 15:04 UTC with Doppler
drifting from +34 267 Hz to +28 822 Hz. The spacecraft was unknown.

```bash
# Step 1 — convert IQ to TDM
python iq_to_tdm.py \
    --input gqrx_20260221_140504_2260790300_125000_fc.sigmf-meta \
    --station SP5LOT --output SP5LOT_20260221_140504.tdm

# Step 2 — identify
python identify_spacecraft.py \
    --tdm SP5LOT_20260221_140504.tdm \
    --freq 2260790300 \
    --station 16.6752,52.3699,0.07 \
    --candidates -155 -85 -152 \
    --plot
```

Output:

```
Loading TDM: SP5LOT_20260221_140504.tdm
  Total measurements : 4393
  Active (|offset|>10 Hz) : 2592
  Active window : 2026-02-21 14:20 – 2026-02-21 15:04 UTC

Testing 3 candidate(s):
  Querying Horizons for ID -155 ... 45 rows  (elev 50.2°–48.5°)
    DC offset = +32535.9 Hz  |  RMS residual = 120.0 Hz  |  n = 2592
  Querying Horizons for ID -85  ... 45 rows  (elev 50.2°–48.5°)
    DC offset = +28471.2 Hz  |  RMS residual = 5071.6 Hz  |  n = 2592
  Querying Horizons for ID -152 ... 45 rows  (elev 50.2°–48.5°)
    DC offset = +28455.5 Hz  |  RMS residual = 1837.3 Hz  |  n = 2592

============================================================
  BEST MATCH: ID -155
  Elevation range : 50.2° – 48.5°
  DC offset (transmit freq – SDR centre) : +32535.9 Hz
  Estimated transmit frequency : 2260.822835 MHz
  RMS residual after DC removal : 120.0 Hz
  Matched pairs : 2592
============================================================
```

**Result: KPLO/Danuri** — an earlier orbital pass from the same day, visible
at 50° elevation from SP5LOT throughout the observation window.

The method was validated against a second confirmed KPLO recording from the
same day (15:47–17:05 UTC), which gave RMS = 96.7 Hz and DC = +32 000 Hz —
consistent with the identification above.

---

## Limitations

- Requires the spacecraft to have a **stable unmodulated carrier** (or a
  detectable residual carrier) during the observation window.
- Horizons step size is 1 minute; TDM measurements are matched to the nearest
  Horizons epoch within 90 s.
- Only spacecraft with JPL Horizons ephemerides can be tested.

---

## Related

- [iq-to-tdm](https://github.com/SP5LOT/iq-to-tdm) — convert SDR IQ recordings to CCSDS TDM.
  Produces the TDM files used as input to this tool.

## License

MIT
