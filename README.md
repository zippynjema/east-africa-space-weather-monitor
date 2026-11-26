<<<<<<< HEAD
=======
# Space Weather Monitor: Solar Flares, Geomagnetic Activity and HF Radio Impact over Kenya / East Africa

A live dashboard, written in Python, that tracks solar flares and geomagnetic conditions and estimates their effect on HF radio over Kenya and East Africa. It answers three questions:

1. Is the Sun producing flares right now, and how strong are they?
2. Is East Africa in daylight when a flare peaks? (HF radio blackouts only affect the sunlit side of Earth.)
3. When a flare was followed by a geomagnetic storm, how long did the effect take to arrive?

It supports, and does not replace, official warnings from the [NOAA Space Weather Prediction Center](https://www.swpc.noaa.gov/).

## Features

- **Real-time flare alerts** from GOES X-ray data, graded on NOAA's radio blackout scale (R1 to R5), on screen, in a log file, and optionally by email or a Slack-style webhook.
- **HF radio impact over Kenya / East Africa**: daylight check and estimated highest affected frequency at the moment each flare peaks, for Nairobi, Mombasa, Kisumu, Lodwar, Kampala, Dar es Salaam and Addis Ababa.
- **Geomagnetic conditions**: solar wind speed and Bz at L1, plus the Dst and Kp indices with plain-language labels.
- **Flare-to-storm delay**, estimated two ways: from NASA analyst-linked events and from live data.
- **Live dashboard page** that reloads every minute, shows how old the data is, and warns if the script stops. A **Save as PDF** button produces a three-page landscape A4 report.

## Screenshots

![Real-time alert status](docs/alert-status.png)
![Seven days of space weather data](docs/charts.png)
![HF radio impact and flare-to-storm delays](docs/hf-impact-and-storm-delays.png)

## Quick start

Requires Python 3.9 or newer and an internet connection.

```bash
git clone https://github.com/<zippynjema>/<east-africa-space-weather-monitor>.git
cd <east-africa-space-weather-monitor>

python -m venv .venv
.venv\Scripts\activate          # Windows  (macOS/Linux: source .venv/bin/activate)
pip install -r requirements.txt

python goes_flare_alert_dashboard.py
```

The script opens `flare_dashboard/dashboard.html` in your browser and keeps refreshing it every 60 seconds. **Leave the script running**: the page only updates while it runs (stop it with Ctrl+C, or the Stop button in Spyder). The dashboard loads Plotly from a CDN, so the browser needs internet access.

### Spyder users

Run the script in its own console. To pass options, use `--args`:

```
%runfile goes_flare_alert_dashboard.py --args "--window 3-day --no-browser"
```

## Usage

| Command | What it does |
|---|---|
| `python goes_flare_alert_dashboard.py` | Live mode (default): refreshes every 60 s and opens the browser |
| `python goes_flare_alert_dashboard.py --once` | One pass, then exit |
| `python goes_flare_alert_dashboard.py --no-browser` | Do not open the browser automatically |
| `python goes_flare_alert_dashboard.py --no-donki` | Skip NASA DONKI (use if it is unreachable from your network) |
| `python goes_flare_alert_dashboard.py --test-email` | Send a test email and exit |

Common options:

| Option | Default | Meaning |
|---|---|---|
| `--window` | `7-day` | X-ray history to pull: `6-hour`, `1-day`, `3-day`, `7-day` |
| `--interval` | `60` | Refresh interval in seconds (minimum 60) |
| `--detect-class` | `C1.0` | Lowest class tracked as a flare |
| `--alert-class` | `M1.0` | Lowest class that triggers an alert |
| `--hf-class` | `M1.0` | Lowest class listed in the Kenya HF impact table |
| `--assoc-class` | `M1.0` | Lowest class used for the flare-to-storm analysis |
| `--kp` / `--dst` | `5` / `-50` | Kp and Dst (nT) thresholds that count as a storm |
| `--min-h` / `--max-h` | `10` / `120` | Hours after a flare in which a storm may be matched to it |
| `--donki-days` | `30` | NASA DONKI look-back in days (maximum 60) |
| `--out-dir` | `flare_dashboard` | Folder for generated files |

Run `python goes_flare_alert_dashboard.py --help` for the full list.

## Configuration (optional)

Email and webhook alerts are **off by default**. To enable them, copy `.env.example` to `.env` and fill in your values. The script reads `.env` from the folder it sits in, and `.env` is excluded from Git, so your secrets stay local.

| Variable | Purpose |
|---|---|
| `FLARE_EMAIL_TO` | Address that receives flare alerts (email stays off if this or the password is empty) |
| `FLARE_SMTP_USER` | Sender login (defaults to the recipient address) |
| `FLARE_SMTP_PASSWORD` | SMTP password. For Gmail use an **App Password** (needs 2-Step Verification) |
| `FLARE_SMTP_HOST` / `FLARE_SMTP_PORT` | Mail server, default `smtp.gmail.com` / `587` |
| `FLARE_WEBHOOK_URL` | Slack-style webhook for alerts |

On the first run, flares that already happened are recorded silently, so you only get alerts for new ones.

## Output files

Everything is written to `flare_dashboard/` (ignored by Git):

| File | Contents |
|---|---|
| `dashboard.html` | The live dashboard page |
| `flares_detected.csv` | Every detected flare with class and times |
| `hf_kenya_impact.csv` | Daylight check and HF impact estimate per flare |
| `latency_heuristic.csv` | Flare-to-storm timing estimated from live data |
| `latency_donki.csv` | Flare-to-storm chain linked by NASA analysts |
| `solar_wind_cache.csv` | Solar wind history saved by this script (kept for 30 days) |
| `alerts.log`, `alert_state.json` | Alert history and the record of flares already alerted |

## How it works

### Data sources

| Data | Source |
|---|---|
| X-ray flux (GOES primary, 0.1 to 0.8 nm) | [NOAA SWPC](https://services.swpc.noaa.gov/json/goes/) |
| Solar wind speed and Bz at L1 | NOAA SWPC real-time solar wind (`json/rtsw`) |
| Kp index | NOAA SWPC **estimated** Kp (3-hourly and 1-minute) |
| Dst index | WDC for Geomagnetism, Kyoto **quick-look** (real-time) Dst, redistributed by NOAA SWPC |
| Flares, CMEs, shocks, geomagnetic storms | [NASA CCMC DONKI](https://ccmc.gsfc.nasa.gov/tools/DONKI/) |

### Alert levels

| Flare class | Level | Effect on HF radio (sunlit side) |
|---|---|---|
| Below C1 | Quiet | None expected |
| C1 to C9.9 | Low | Normally unaffected |
| M1 to M4.9 | R1 Minor | Weak or minor degradation; occasional loss of contact |
| M5 to M9.9 | R2 Moderate | Limited blackout; loss of contact for tens of minutes |
| X1 to X9.9 | R3 Strong | Wide-area blackout; loss of contact for about an hour |
| X10 to X19.9 | R4 Severe | Blackout on most of the sunlit side for 1 to 2 hours |
| X20 and above | R5 Extreme | Complete blackout on the sunlit side for several hours |

### HF impact estimate

The estimate follows NOAA's [D-Region Absorption Prediction (D-RAP)](https://update-swpc.woc.noaa.gov/content/global-d-region-absorption-prediction-documentation) method. At the point where the sun is overhead, the highest affected frequency (HAF) is

```
HAF (MHz) = 10 x log10(X-ray flux in W/m2) + 65
```

and away from that point it is reduced by `cos(zenith angle)^0.75`. It is zero when the sun is below the horizon. The sun's elevation is computed for each site at the flare's peak time. The result is a 1 dB, vertical-path absorption estimate: slanted, long-distance paths are affected at higher frequencies, and the upper atmosphere stays sunlit a little after ground sunset, so conditions just after sunset are a grey zone.

### Flare-to-storm delay

- **Method A, NASA DONKI:** follows the analyst-linked chain flare, CME, interplanetary shock, geomagnetic storm, and reports the delay at each step. This is the reliable attribution.
- **Method B, live data:** for each flare, searches `--min-h` to `--max-h` hours later for a solar wind shock, Kp at or above the threshold, and Dst at or below the threshold. It can match the wrong event, so it also reports how many other flares could explain the same one.

Most flares do not cause geomagnetic storms. Storms usually need an Earth-directed CME, which takes roughly 15 hours to 4 days to arrive.

## Limitations

- HF figures are model-based estimates, not measurements, and use conditions at the flare peak.
- Kp is NOAA's near-real-time estimate, not the official index from [GFZ Potsdam](https://www.gfz.de/en/section/geomagnetism/data-products-services/geomagnetic-kp-index). Dst is a provisional quick-look value that can change.
- NOAA's solar wind feed only holds about 24 hours, so the script builds its own history from the time it first runs. Earlier days, and times when the script was not running, appear as gaps.
- Kp is reported every 3 hours, which limits the timing precision of the live-data flare-to-storm estimate.
- Data feeds change from time to time (NOAA retired its 7-day solar wind files and NASA moved the DONKI API in 2026), so the script may need maintenance.
- This is not an official warning service.

## Deploying on a server

The page only updates while the script runs, so a public deployment needs an always-on host running the script (for example with `--no-browser`) and a web server that serves the `flare_dashboard/` folder. Before publishing the page, read the data terms below, and add a footer that names the data sources and states that the HF figures are estimates.

## Data sources and acknowledgements

- NOAA Space Weather Prediction Center: GOES X-ray, solar wind, Kp and Dst data. https://www.swpc.noaa.gov/
- NASA Community Coordinated Modeling Center: DONKI database. https://ccmc.gsfc.nasa.gov/tools/DONKI/
- World Data Center for Geomagnetism, Kyoto: Dst index. https://wdc.kugi.kyoto-u.ac.jp/ . Real-time (quick-look) data is provisional. Kyoto asks users to acknowledge the data suppliers and WDC Kyoto, to contact them before using real-time data in publications or presentations, and not to use the indices commercially. See their [data service page](https://wdc.kugi.kyoto-u.ac.jp//wdc/Sec3.html).
- NOAA D-Region Absorption Prediction documentation: https://update-swpc.woc.noaa.gov/content/global-d-region-absorption-prediction-documentation

## Licence

This project is licensed under the `Apache License` - see the [LICENSE](LICENSE) file for details.

## Author

Zipporah Njeri, Kenya Space Agency. Feedback and collaboration are welcome: [www.linkedin.com/in/zipporahnjeri]

## Disclaimer

This software is provided for information and research purposes only, without warranty of any kind. Use official NOAA SWPC products for operational space weather warnings.
>>>>>>> 294bf1c (Add GOES flare alert dashboard and project files)
