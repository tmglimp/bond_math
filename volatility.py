# volatility.py
import requests
import numpy as np
import pandas as pd
from datetime import date, timedelta
import config

today = date.today()
yesterday = today - timedelta(days=1)
hundred_one_before_yday = yesterday - timedelta(days=101)  # t-101

config.fred_key
baseurl = "https://api.stlouisfed.org/fred/series/observations"
url_ext = ("series_id={series_id}"
    "&api_key={fred_key}"
    "&observation_start={start}"
    "&observation_end={end}"
    "&file_type=json")

def fetch_yields_df(series_id: str) -> pd.DataFrame:
    full_url = f"{baseurl}?{url_ext.format(series_id=series_id, fred_key=config.fred_key, start=hundred_one_before_yday.isoformat(), end=yesterday.isoformat())}"
    resp = requests.get(full_url, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    obs = data.get("observations", [])
    rows = [
        (obs_i["date"], float(obs_i["value"]))
        for obs_i in obs
        if obs_i.get("value") not in (".", None, "")
    ]
    df = pd.DataFrame(rows, columns=["date", "yield_pct"])
    if not df.empty:
        df["date"] = pd.to_datetime(df["date"], errors="coerce")
        df = df.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)
    return df

def yield_log_vol(df: pd.DataFrame, series_id: str) -> tuple[pd.DataFrame, float]:
    yields = df.copy()
    yields["series"] = series_id
    if yields.empty:
        yields["log_yield"] = pd.Series(dtype="float64")
        return yields, np.nan

    y = yields["yield_pct"].where(yields["yield_pct"] > 0)
    yields["log_yield"] = round((np.log(y)), 6)
    log_vals = yields["log_yield"].dropna().to_numpy()
    print(log_vals)
    n = len(log_vals)
    m = n-2
    print(n)
    if n > 1:
        std_log = float(np.std(log_vals, ddof=(m)))
        adj = round((std_log * np.sqrt(n)),6)
    else:
        adj = np.nan
        yields["adj_log_yield"] = adj

    return yields, adj

def derive_vol(series_list=None):

    if series_list is None:
        series_list = ["THREEFY1", "THREEFY2", "THREEFY3","THREEFY4", "THREEFY5", "THREEFY6", "THREEFY7", "DGS1","DGS2","DGS3","DGS5","DGS7"]

    print(f"Date range queried: {hundred_one_before_yday} → {yesterday}")
    for sid in series_list:
        try:
            df = fetch_yields_df(sid)
            df, logvol_adj = yield_log_vol(df, sid)
            setattr(config, f"{sid}_ln_y100_std", round((logvol_adj), 6))
            adj_val = getattr(config, f"{sid}_ln_y100_std")
            print(f"{sid}: modified std={adj_val}")
        except Exception as e:
            print(f"[ERROR] {sid}: {e}")
            setattr(config, f"{sid}_ln_y100_std", np.nan)

if __name__ == "__main__":
    derive_vol()
