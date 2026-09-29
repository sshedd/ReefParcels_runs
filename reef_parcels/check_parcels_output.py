# check_parcels_output.py

"""
Simple quality-check script for a Parcels larval dispersal model output.

Input:
    output.zarr

This script checks:
1. Whether the Zarr file opens correctly
2. What variables are inside
3. Whether longitude, latitude, and time look reasonable
4. Whether particles moved
5. Whether there are missing values
6. Basic trajectory plots
7. Start vs. final particle locations
"""

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# ------------------------------------------------------------
# 1. Load the model output
# ------------------------------------------------------------

zarr_file = "output.zarr"

print("\nOpening Zarr file...")
ds = xr.open_zarr(zarr_file)

print("\nFile opened successfully.")
print("\nDataset summary:")
print(ds)


# ------------------------------------------------------------
# 2. Check what variables are present
# ------------------------------------------------------------

print("\nVariables in the dataset:")
for var in ds.data_vars:
    print(f" - {var}")

print("\nDimensions:")
print(ds.dims)


# ------------------------------------------------------------
# 3. Check for required variables
# ------------------------------------------------------------

required_vars = ["lon", "lat", "time"]

print("\nChecking required variables...")

for var in required_vars:
    if var not in ds:
        raise ValueError(f"Missing required variable: {var}")
    else:
        print(f"Found variable: {var}")


# ------------------------------------------------------------
# 4. Basic missing-value checks
# ------------------------------------------------------------

print("\nChecking for missing values...")

for var in ["lon", "lat"]:
    n_missing = ds[var].isnull().sum().compute().item()
    print(f"{var}: {n_missing} missing values")


# ------------------------------------------------------------
# 5. Check coordinate ranges
# ------------------------------------------------------------

print("\nChecking longitude and latitude ranges...")

lon_min = float(ds["lon"].min().compute())
lon_max = float(ds["lon"].max().compute())
lat_min = float(ds["lat"].min().compute())
lat_max = float(ds["lat"].max().compute())

print(f"Longitude range: {lon_min:.4f} to {lon_max:.4f}")
print(f"Latitude range:  {lat_min:.4f} to {lat_max:.4f}")

if lon_min < -180 or lon_max > 360:
    print("Warning: longitude values may be outside expected range.")

if lat_min < -90 or lat_max > 90:
    print("Warning: latitude values are outside expected range.")


# ------------------------------------------------------------
# 6. Check start and final locations
# ------------------------------------------------------------

print("\nExtracting start and final particle positions...")

start = ds.isel(obs=0)
end = ds.isel(obs=-1)

start_lon = start["lon"].values
start_lat = start["lat"].values
end_lon = end["lon"].values
end_lat = end["lat"].values


# ------------------------------------------------------------
# 7. Estimate particle displacement
# ------------------------------------------------------------

print("\nEstimating simple displacement...")

displacement_degrees = np.sqrt(
    (end_lon - start_lon) ** 2 +
    (end_lat - start_lat) ** 2
)

mean_disp = np.nanmean(displacement_degrees)
max_disp = np.nanmax(displacement_degrees)
n_no_move = np.sum(displacement_degrees == 0)

print(f"Mean displacement, in degrees: {mean_disp:.6f}")
print(f"Max displacement, in degrees:  {max_disp:.6f}")
print(f"Particles with zero displacement: {n_no_move}")

if mean_disp == 0:
    print("Warning: particles may not have moved.")

print(ds["lon"].isel(trajectory=0).values[:20])
print(ds["lat"].isel(trajectory=0).values[:20])

# ------------------------------------------------------------
# 8. Plot a subset of trajectories
# ------------------------------------------------------------

print("\nPlotting example trajectories...")

n_particles_to_plot = min(50, ds.sizes["trajectory"])

plt.figure(figsize=(8, 6))

for i in range(n_particles_to_plot):
    plt.plot(
        ds["lon"].isel(trajectory=i),
        ds["lat"].isel(trajectory=i),
        alpha=0.4
    )

plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Example Larval Trajectories")
plt.grid(True)
plt.tight_layout()
plt.savefig("trajectory_check.png", dpi=300)
plt.close()

print("Saved trajectory plot as trajectory_check.png")


# ------------------------------------------------------------
# 9. Plot release vs final locations
# ------------------------------------------------------------

print("\nPlotting release and final locations...")

plt.figure(figsize=(8, 6))

plt.scatter(start_lon, start_lat, s=15, label="Release locations", alpha=0.7)
plt.scatter(end_lon, end_lat, s=15, label="Final locations", alpha=0.7)

plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Release vs Final Particle Locations")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("start_vs_final_locations.png", dpi=300)
plt.close()

print("Saved start/final plot as start_vs_final_locations.png")


# ------------------------------------------------------------
# 10. Export a simple summary table
# ------------------------------------------------------------

print("\nSaving final particle summary table...")

summary = pd.DataFrame({
    "trajectory": ds["trajectory"].values,
    "start_lon": start_lon,
    "start_lat": start_lat,
    "final_lon": end_lon,
    "final_lat": end_lat,
    "displacement_degrees": displacement_degrees
})

summary.to_csv("particle_summary.csv", index=False)

print("Saved particle summary as particle_summary.csv")


# ------------------------------------------------------------
# 11. Final message
# ------------------------------------------------------------

print("\nDone.")
print("Check these output files:")
print(" - trajectory_check.png")
print(" - start_vs_final_locations.png")
print(" - particle_summary.csv")
