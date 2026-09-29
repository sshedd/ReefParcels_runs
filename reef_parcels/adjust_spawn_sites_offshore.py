import glob
import math
import numpy as np
import pandas as pd
import xarray as xr

SITE_FILE = "spawn_files/my_spawn_sites.csv"
CURRENT_FILE = sorted(glob.glob("current_files/1km_MITgcm/*.cdf"))[0]
OUT_FILE = "data/spawn_sites_offshore_qc.csv"

STEP_KM = 0.20
MAX_SHIFT_KM = 5.0


def offset_point(lon, lat, bearing_deg, distance_km):
    """Move distance_km along bearing_deg."""
    R = 6371.0
    d = distance_km / R

    b = math.radians(bearing_deg)
    lat1 = math.radians(lat)
    lon1 = math.radians(lon)

    lat2 = math.asin(
        math.sin(lat1) * math.cos(d)
        + math.cos(lat1) * math.sin(d) * math.cos(b)
    )

    lon2 = lon1 + math.atan2(
        math.sin(b) * math.sin(d) * math.cos(lat1),
        math.cos(d) - math.sin(lat1) * math.sin(lat2),
    )

    return math.degrees(lon2) % 360.0, math.degrees(lat2)


# ------------------------------------------------------------
# Read original empirical sampling sites
# ------------------------------------------------------------

sites = pd.read_csv(SITE_FILE)

sites["original_lat"] = sites["lat"].astype(float)
sites["original_lon"] = sites["lon"].astype(float)

sites["model_lon"] = np.where(
    sites["original_lon"] < 0,
    sites["original_lon"] + 360.0,
    sites["original_lon"],
)


# ------------------------------------------------------------
# Read actual 1-km hydrodynamic wet/dry grid
# ------------------------------------------------------------

ds = xr.open_dataset(CURRENT_FILE, decode_times=False)

lon_t = ds["Longitude_t"].values
lat_t = ds["Latitude_t"].values

surface = ds["temp"].isel(Time=0, Depth_t=0).values
wet = np.isfinite(surface)

ds.close()


def nearest_grid_index(lon, lat):
    ix = int(np.argmin(np.abs(lon_t - lon)))
    iy = int(np.argmin(np.abs(lat_t - lat)))
    return iy, ix


def is_wet(lon, lat):
    iy, ix = nearest_grid_index(lon, lat)
    return bool(wet[iy, ix])


# ------------------------------------------------------------
# Find wet-side boundary of the current grid
# ------------------------------------------------------------

boundary = np.zeros_like(wet, dtype=bool)

boundary[1:, :] |= wet[1:, :] & ~wet[:-1, :]
boundary[:-1, :] |= wet[:-1, :] & ~wet[1:, :]
boundary[:, 1:] |= wet[:, 1:] & ~wet[:, :-1]
boundary[:, :-1] |= wet[:, :-1] & ~wet[:, 1:]

by, bx = np.where(boundary)

boundary_lon = lon_t[bx]
boundary_lat = lat_t[by]

print("Hydrodynamic wet-boundary cells:", len(boundary_lon))


# ------------------------------------------------------------
# Relocate sites
# ------------------------------------------------------------

results = []

for _, row in sites.iterrows():

    lon = float(row["model_lon"])
    lat = float(row["original_lat"])

    release_lon = np.nan
    release_lat = np.nan
    shift = np.nan
    bearing = np.nan
    status = "FAILED"

    # No relocation needed
    if is_wet(lon, lat):

        release_lon = lon
        release_lat = lat
        shift = 0.0
        status = "already_wet"

    else:

        # Distance from site to every wet-side boundary cell
        dx = (
            (boundary_lon - lon)
            * 111.32
            * np.cos(np.radians(lat))
        )

        dy = (
            (boundary_lat - lat)
            * 110.57
        )

        distance = np.sqrt(dx**2 + dy**2)

        # Nearest wet-side model coastline
        i = np.argmin(distance)

        # Direction from empirical site toward that coastline
        bearing = (
            np.degrees(
                np.arctan2(dx[i], dy[i])
            ) + 360
        ) % 360

        # Walk outward along that direction until currents are reached
        for dist in np.arange(
            STEP_KM,
            MAX_SHIFT_KM + STEP_KM,
            STEP_KM,
        ):

            test_lon, test_lat = offset_point(
                lon,
                lat,
                bearing,
                dist,
            )

            if is_wet(test_lon, test_lat):

                # Move one additional step into the wet region
                safe_dist = min(
                    dist + STEP_KM,
                    MAX_SHIFT_KM,
                )

                safe_lon, safe_lat = offset_point(
                    lon,
                    lat,
                    bearing,
                    safe_dist,
                )

                if is_wet(safe_lon, safe_lat):
                    release_lon = safe_lon
                    release_lat = safe_lat
                    shift = safe_dist
                else:
                    release_lon = test_lon
                    release_lat = test_lat
                    shift = dist

                status = "shifted"
                break

    result = row.to_dict()

    result.update(
        {
            "release_lat": release_lat,
            "release_lon": release_lon,
            "shift_km": shift,
            "offshore_bearing": bearing,
            "current_valid": (
                is_wet(release_lon, release_lat)
                if np.isfinite(release_lon)
                and np.isfinite(release_lat)
                else False
            ),
            "relocation_status": status,
        }
    )

    results.append(result)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

out = pd.DataFrame(results)
out.to_csv(OUT_FILE, index=False)


# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

print("\nOFFSHORE RELOCATION QC")
print("----------------------")
print("Sites:", len(out))
print(
    "Already wet:",
    (out["relocation_status"] == "already_wet").sum()
)
print(
    "Shifted:",
    (out["relocation_status"] == "shifted").sum()
)
print(
    "FAILED:",
    (out["relocation_status"] == "FAILED").sum()
)
print(
    "Invalid final current locations:",
    (~out["current_valid"]).sum()
)

moved = out[out["relocation_status"] == "shifted"]

print("\nShift distance (km):")
print(moved["shift_km"].describe())

print("\nLargest shifts:")

cols = [
    "site",
    "original_lat",
    "original_lon",
    "release_lat",
    "release_lon",
    "shift_km",
    "offshore_bearing",
    "relocation_status",
]

print(
    out.sort_values(
        "shift_km",
        ascending=False
    )[cols].head(20).to_string(index=False)
)

print("\nSaved:", OUT_FILE)
