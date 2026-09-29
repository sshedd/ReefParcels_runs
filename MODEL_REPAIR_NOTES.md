# ReefParcels Model Repair Notes

## Purpose

This repository contains the working repair and adaptation of the ReefParcels
larval dispersal model for a Pocillopora acuta connectivity analysis using
empirical sampling sites around Oahu.

The starting model was inherited from an existing ReefParcels workflow and had
been modified by multiple users. The immediate goal of this repair is to:

1. obtain a stable and reproducible model run,
2. use the available 1-km MITgcm current data,
3. release larvae from the empirical sampling sites,
4. adapt reproductive settings for P. acuta,
5. preserve the existing ReefParcels biological and settlement framework where
   it has not yet been intentionally changed, and
6. document all departures from the inherited model.

The current primary working script is:

    reef_parcels/parcels_pocillopora_complex_2023_REPAIR.py

The repair was based primarily on:

    reef_parcels/parcels_pocillopora_complex_2023_TRY2.py

The TRY2 version was retained rather than overwritten so that the repair process
can be reconstructed.

---

# 1. Git/repository repair

The `reef_parcels` directory was originally represented in the parent Git
repository as a gitlink/submodule entry rather than as ordinary tracked source
files.

The referenced gitlink commit was not available in the current repository, and
the directory itself did not contain an independent `.git` repository.

The stale gitlink was therefore removed from the Git index with:

    git rm --cached reef_parcels

The actual `reef_parcels` directory was preserved on disk and its Python source
files are now being tracked directly by this repository.

Large environmental datasets, model outputs, Zarr stores, Python cache files,
and repair backup files are excluded using `.gitignore`.

---

# 2. Ocean current data

## Available current data

Two MITgcm datasets were found in the inherited project.

### 1-km dataset

Location:

    current_files/1km_MITgcm/

Number of CDF files:

    852

Time coverage:

    2011-04-01 12:00
    through
    2013-07-30 12:00

Approximate spatial coverage:

    Longitude: 198.2000 to 205.9857 degrees E
    Latitude:   17.0000 to 22.1901 degrees N

The files contain staggered-grid velocity components:

    u(Time, Depth_t, Latitude_t, Longitude_u)
    v(Time, Depth_t, Latitude_v, Longitude_t)
    w(Time, Depth_w, Latitude_t, Longitude_t)

### 4-km dataset

Location:

    current_files/4km_MITgcm_unused/

Number of CDF files:

    2107

Approximate size:

    697 GB

Time coverage:

    2011-03-19 through 2016-12-23

Approximate spatial coverage:

    Longitude: 175.0400 to 210.0341 degrees E
    Latitude:   14.9200 to 35.0404 degrees N

All 852 dates represented in the 1-km dataset were found to have matching
calendar dates in the 4-km dataset.

However, the two simulations use different raw model time coordinates/time
origins. Combining them as nested Parcels velocity fields produced
TimeExtrapolationError failures.

Matching the files by filename was also determined to be invalid because the
filenames represent timestep identifiers from different simulations.

## Current decision

The repaired model currently uses the 1-km velocity fields only.

The 4-km velocity fallback is disabled.

This is intentional and should not be reversed without correctly reconciling the
time coordinates of the two MITgcm simulations.

All empirical release sites are geographically inside the 1-km domain.

A future production version may restore the 4-km fields if particles need to be
followed after leaving the 1-km domain.

---

# 3. Empirical sampling/release sites

Original empirical sites:

    spawn_files/my_spawn_sites.csv

Number of sites:

    193

Columns:

    lat
    lon
    site

Approximate coordinate range:

    Longitude: -158.246505 to -157.676
    Latitude:   21.2535 to 21.699

The original ReefParcels habitat-based spawning-site generator is not being used
to determine release locations for this analysis.

Instead, particles are released from the empirical sampling sites.

---

# 4. Offshore relocation of release sites

Some empirical sampling locations occur on land or too close to the coastline
to interpolate the MITgcm velocity fields reliably.

A preprocessing script was created:

    reef_parcels/adjust_spawn_sites_offshore.py

The script uses the 1-km MITgcm wet/dry grid.

For each empirical site:

1. determine whether the site is already in a wet model cell;
2. if dry, locate the nearest wet-side ocean boundary cell;
3. calculate the direction from the site toward that wet cell;
4. move the release position in 0.2-km increments toward the ocean;
5. stop when a wet model location is reached;
6. move an additional 0.2 km in the same direction if that position remains wet.

The original site coordinates and site IDs are retained. Separate release
coordinates are generated so the empirical sampling locations are not
overwritten.

Output:

    data/spawn_sites_offshore_qc.csv

Results:

    Total sites:       193
    Already wet:        29
    Relocated:         164
    Failed:              0

Relocation distance:

    Mean:    ~1.03 km
    Median:   1.00 km
    Maximum:  3.20 km

All 193 resulting release positions were successfully tested against actual
Parcels U, V, and W interpolation at 5 m depth.

---

# 5. Particle initialization from empirical sites

The repaired model reads:

    data/spawn_sites_offshore_qc.csv

and uses:

    site
    release_lat
    release_lon

Particles are initialized at:

    depth = 5.0 m

The release code preserves the empirical `site` value as the particle `siteID`.

For each release timestep, the number of particles initialized at each site is
controlled by:

    n_eggs

In this implementation, `n_eggs` should be interpreted as a COMPUTATIONAL
PARTICLE MULTIPLIER:

    particles released per site per release timestep

It is not interpreted as literal biological egg or planula production.

This interpretation is consistent with the original ReefParcels workflow, in
which reproductive rates were calibrated to provide sufficient particles for
connectivity analysis rather than being treated as measured biological
fecundity.

Current setting:

    n_eggs = 1

This is intentionally the minimum particle multiplier while model performance
and release scheduling are being validated.

---

# 6. P. acuta reproductive configuration

The current model is being adapted for Pocillopora acuta.

Current settings:

    n_eggs = 1

    typeSpawner = 1

    monthSpawner = False
    spawnMonths = "NA"

    timedSpawner = True
    timeSpawnBegin = 03:00
    timeSpawnEnd = 06:00

    peakSpawner = False
    peakSpawnMonths = "NA"

    moonSpawner = False

## Interpretation

`typeSpawner = 1` uses the ReefParcels brooder pathway.

`monthSpawner = False` represents year-round reproduction.

`moonSpawner = False` means no lunar-phase restriction is currently applied.

The value stored in `spawnPhase` is ignored when `moonSpawner = False`.

## Literature basis and modeling simplification

Lam et al. (2023), Diversity 15(2):218, was used to inform the reproductive
timing configuration.

The study supports:

- reproduction throughout the year;
- predominantly nocturnal planula release;
- variable reproductive activity among colonies and lunar months;
- no simple lunar-phase restriction suitable for directly imposing as a fixed
  ReefParcels release window.

The current 03:00-06:00 release window is a MODEL SIMPLIFICATION based on the
late-night concentration of release shown in the study's hourly release data.

It should NOT be described as an explicit 03:00-06:00 biological release window
reported by the authors.

The model currently does not randomly select a fixed number of reproductive
days per lunar month.

---

# 7. Spawning-event code modification

File modified:

    reef_parcels/lanaau_build_spawn_2023.py

The inherited `getSpawningEvents()` logic did not contain a branch for the
combination:

    moonSpawner = False
    monthSpawner = False
    peakSpawner = False
    timedSpawner = True

This combination is required for the current P. acuta configuration:
year-round spawning with a restricted time-of-day window but no lunar or
seasonal restriction.

The following logic was added:

    # no moon, no month, no peak, time of day
    # P. acuta - year-round brooded larval release
    if not moonSpawner and not monthSpawner and not peakSpawner and timedSpawner:
        spawn_starts = [x + timeSpawnBegin for x in spawningdays]

        if timeSpawnBegin > timeSpawnEnd:
            # release window crosses midnight
            spawn_stops = [
                x + datetime.timedelta(hours=24) + timeSpawnEnd
                for x in spawningdays
            ]
        else:
            spawn_stops = [x + timeSpawnEnd for x in spawningdays]

With the current 03:00-06:00 window, the event generator produces alternating:

    3-hour spawning segment
    21-hour non-spawning segment

for each day in the available forcing period.

The full forcing period currently produces:

    1702 event segments
    851 spawning days

---

# 8. Advection NaN bug and repair

File modified:

    reef_parcels/lanaau_functions_complex.py

A major particle-loss bug was identified in:

    AdvectionRK4_3D_Lanaau

The inherited kernel counted nonzero U, V, and W velocity samples during the RK4
calculation and divided the accumulated velocity by those counts.

If all four RK4 samples for one velocity component were exactly zero, the
corresponding count was zero.

This caused a division by zero:

    0 / 0

which produced NaN coordinates.

The particle could then disappear from valid output without being classified as
a normal biological or boundary-related death.

The kernel was changed so that a coordinate is updated only when the number of
valid/nonzero samples for that component is greater than zero:

    if n_valid_u > 0:
        particle.lon += (
            (u1 + 2*u2 + 2*u3 + u4)
            / n_valid_u * particle.dt
        )

    if n_valid_v > 0:
        particle.lat += (
            (v1 + 2*v2 + 2*v3 + v4)
            / n_valid_v * particle.dt
        )

    if n_valid_w > 0:
        particle.depth += (
            (w1 + 2*w2 + 2*w3 + w4)
            / n_valid_w * particle.dt
        )

If all sampled velocities for a component are zero, displacement in that
component is therefore zero instead of NaN.

This repair eliminated the observed internal-lifespan NaN trajectories in the
subsequent tests.

---

# 9. Parcels compatibility / boundary diagnostics

The current environment uses:

    Parcels 3.1.4
    Python 3.11

The repaired workflow uses the current Parcels `StatusCode` interface rather
than relying on the older `ErrorCode` API where applicable.

The model currently includes explicit out-of-bounds diagnostics.

Out-of-bounds particles are classified using:

    isDead = 1
    howDead = 1

and their final coordinates are recorded before deletion.

Current `howDead` interpretation used during QC:

    0 = alive / no recorded death
    1 = out of bounds
    2 = beached
    3 = random mortality

---

# 10. Current kernel chain

The repaired complex model currently executes:

    InitDepth
    GrowFishComplex
    KillFish
    SettleFish2
    AdvectionRK4_3D_Lanaau
    RWDiffusion_Lanaau
    Behavior
    Buoyancy
    HandleBoundaries
    KeepInOcean
    CheckOutOfBounds

Parcels currently emits warnings because some inherited kernels directly modify
particle coordinates.

Those kernels have NOT yet been broadly rewritten because the immediate priority
is preserving the inherited model behavior while establishing a stable working
baseline.

---

# 11. Validation tests

## Five-day diagnostic test

A previous repaired test produced:

    Trajectories:       7,913
    Observations:         120
    Unique site IDs:      193
    Never valid:            0

Final-status QC:

    Alive:              7,087
    Dead:                 826

Deaths:

    Out of bounds:          0
    Beached:              590
    Random mortality:     236
    Maximum PLD:            0

Settled stage-2 particles:

    3,420

Critically:

    trajectories with NaN longitude/latitude
    inside their recorded lifespan = 0

This test was used to verify the advection divide-by-zero repair.

---

## P. acuta 48-hour validation test

A 48-hour test was run using the current P. acuta 03:00-06:00 release schedule.

Results:

    Trajectories:       5,983
    Observations:          48
    Unique sites:         193
    Never valid:            0

Internal NaN trajectories:

    0

Mortality:

    Alive:              5,822
    Dead:                 161

Death causes:

    howDead = 0:        5,822
    howDead = 2:           70
    howDead = 3:           91

Out-of-bounds deaths:

    0

This indicates that the repaired model can currently initialize all 193 sites
and execute the core model without the previous NaN/disappearing-particle
failure.

---

# 12. Particle-count implications

With:

    193 sites
    n_eggs = 1
    12-minute release interval
    03:00-06:00 daily release window
    851 spawning days

the current configuration is estimated to generate approximately:

    2.63 million particle trajectories

over the full forcing period.

This is a computational release count, not an estimate of actual P. acuta
larval production.

The number of simultaneously active particles may become large because cohorts
can overlap for many days. Runtime therefore should not be extrapolated only
from the short validation runs.

---

# 13. Current performance benchmark

Before attempting the complete 851-day simulation, a longer benchmark is being
performed.

The execution loop was temporarily changed from:

    events[:4]

to:

    events[:60]

Because each day currently consists of two event segments, 60 event segments
represent approximately 30 model days.

The purpose of this benchmark is to assess:

- runtime as particle cohorts accumulate;
- memory behavior;
- output size;
- particle stability;
- NaN occurrence;
- out-of-bounds behavior; and
- whether a full production simulation is computationally practical.

Results should be added here when the benchmark finishes.

---

# 14. Important unresolved issues

## 4-km current fallback

The 4-km velocity fields are currently disabled because their time coordinates
are incompatible with the 1-km simulation when directly nested in Parcels.

This remains unresolved.

## Legacy settlement/habitat fields

Although empirical release sites now replace the inherited spawning-site
generation, the model still contains inherited P. meandrina settlement-related
inputs, including files such as:

    p_meandrina_settle.nc
    p_mea_distance_final_2023.nc
    p_mea_angle_final_2023.nc
    p_meandrina_spawn.nc

These have NOT yet been fully evaluated or replaced for P. acuta.

The release-site change should therefore not be interpreted as a complete
species conversion of every biological/habitat component in the model.

## Other biological parameters

Several complex-model biological parameters are still inherited from the
existing ReefParcels configuration.

These include parameters related to:

- pelagic larval duration;
- mortality;
- growth;
- buoyancy;
- swimming behavior;
- settlement competency; and
- maximum adult depth.

They should be reviewed individually before describing the final model as a
fully P. acuta-specific biological parameterization.

## PLD/hatch-time handling

The current complex model adds `hatchTime` to the PLD constants before they are
stored in the Parcels FieldSet.

Because the model uses the brooder pathway for P. acuta, the interaction between
brooder stage initialization, hatchTime, minimum PLD, and settlement competency
should be verified before the production run.

No change has yet been made to this behavior.

## Full production runtime

The complete 851-day run has not yet been benchmarked successfully with the
current repaired configuration.

Do not launch the complete production simulation until the longer benchmark has
been evaluated.

---

# 15. Files currently central to the repaired workflow

Primary model:

    reef_parcels/parcels_pocillopora_complex_2023_REPAIR.py

Baseline used for repair:

    reef_parcels/parcels_pocillopora_complex_2023_TRY2.py

Modified spawning logic:

    reef_parcels/lanaau_build_spawn_2023.py

Modified complex kernels:

    reef_parcels/lanaau_functions_complex.py

Offshore-site preprocessing:

    reef_parcels/adjust_spawn_sites_offshore.py

Original empirical sites:

    spawn_files/my_spawn_sites.csv

Processed release sites:

    data/spawn_sites_offshore_qc.csv

---

# 16. Next steps

Before a production run:

1. evaluate the 30-day benchmark;
2. record benchmark runtime, trajectory count, deaths, settlement, NaNs, and
   output size;
3. verify settlement behavior and settlement output variables;
4. review inherited biological parameters for P. acuta;
5. resolve or explicitly accept the remaining P. meandrina habitat/settlement
   dependencies;
6. decide whether the 1-km domain alone is sufficient or whether the 4-km
   fallback must be repaired;
7. estimate full-run runtime and storage requirements from the longer benchmark;
8. preserve the tested production configuration in Git before launching the
   full run.

---

# 17. Current status

At the time of this documentation:

- all 193 empirical sites have valid offshore-adjusted release positions;
- 1-km MITgcm velocity interpolation works at all release sites;
- the 4-km velocity fallback is intentionally disabled;
- the RK4 NaN/division-by-zero failure has been repaired;
- year-round timed P. acuta spawning events are functioning;
- `n_eggs = 1` is being used as a computational particle multiplier;
- the 48-hour P. acuta test completed with zero internal-lifespan NaN
  trajectories and zero out-of-bounds deaths;
- a longer benchmark is the next validation step;
- the model should still be considered a repaired/testing version rather than a
  finalized production P. acuta model.
