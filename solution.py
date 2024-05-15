import sys
import os, csv
from datetime import datetime, timezone, timedelta
import pandas as pd
import numpy as np
import navpy
import pandas
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

# sys.path.insert(0, "/gnss-analysis-main/gnssutils/")
from ephemeris_manager import EphemerisManager

parent_directory = os.path.split(os.getcwd())[0]
ephemeris_data_directory = os.path.join(parent_directory, "data")
sys.path.insert(0, parent_directory)

# Get path to sample file in data directory, which is located in the parent directory of this notebook
input_filepath = os.path.join("data/", "gnss_log_2024_04_13_19_53_33.txt")

with open(input_filepath) as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row[0][0] == "#":
            if "Fix" in row[0]:
                android_fixes = [row[1:]]
            elif "Raw" in row[0]:
                measurements = [row[1:]]
        else:
            if row[0] == "Fix":
                android_fixes.append(row[1:])
            elif row[0] == "Raw":
                measurements.append(row[1:])

android_fixes = pd.DataFrame(android_fixes[1:], columns=android_fixes[0])
measurements = pd.DataFrame(measurements[1:], columns=measurements[0])

# Format satellite IDs
measurements.loc[measurements["Svid"].str.len() == 1, "Svid"] = (
    "0" + measurements["Svid"]
)
measurements.loc[measurements["ConstellationType"] == "1", "Constellation"] = "G"
measurements.loc[measurements["ConstellationType"] == "3", "Constellation"] = "R"
measurements["SvName"] = measurements["Constellation"] + measurements["Svid"]

# Remove all non-GPS measurements
measurements = measurements.loc[measurements["Constellation"] == "G"]

# Convert columns to numeric representation
measurements["Cn0DbHz"] = pd.to_numeric(measurements["Cn0DbHz"])
measurements["TimeNanos"] = pd.to_numeric(measurements["TimeNanos"])
measurements["FullBiasNanos"] = pd.to_numeric(measurements["FullBiasNanos"])
measurements["ReceivedSvTimeNanos"] = pd.to_numeric(measurements["ReceivedSvTimeNanos"])
measurements["PseudorangeRateMetersPerSecond"] = pd.to_numeric(
    measurements["PseudorangeRateMetersPerSecond"]
)
measurements["ReceivedSvTimeUncertaintyNanos"] = pd.to_numeric(
    measurements["ReceivedSvTimeUncertaintyNanos"]
)

# A few measurement values are not provided by all phones
# We'll check for them and initialize them with zeros if missing
if "BiasNanos" in measurements.columns:
    measurements["BiasNanos"] = pd.to_numeric(measurements["BiasNanos"])
else:
    measurements["BiasNanos"] = 0
if "TimeOffsetNanos" in measurements.columns:
    measurements["TimeOffsetNanos"] = pd.to_numeric(measurements["TimeOffsetNanos"])
else:
    measurements["TimeOffsetNanos"] = 0

# print(measurements.columns)

measurements["GpsTimeNanos"] = measurements["TimeNanos"] - (
    measurements["FullBiasNanos"] - measurements["BiasNanos"]
)
gpsepoch = datetime(1980, 1, 6, 0, 0, 0)
measurements["UnixTime"] = pd.to_datetime(
    measurements["GpsTimeNanos"], utc=True, origin=gpsepoch
)
measurements["UnixTime"] = measurements["UnixTime"]

# Split data into measurement epochs
measurements["Epoch"] = 0
measurements.loc[
    measurements["UnixTime"] - measurements["UnixTime"].shift()
    > timedelta(milliseconds=200),
    "Epoch",
] = 1
measurements["Epoch"] = measurements["Epoch"].cumsum()


WEEKSEC = 604800
LIGHTSPEED = 2.99792458e8

# This should account for rollovers since it uses a week number specific to each measurement

measurements["tRxGnssNanos"] = (
    measurements["TimeNanos"]
    + measurements["TimeOffsetNanos"]
    - (measurements["FullBiasNanos"].iloc[0] + measurements["BiasNanos"].iloc[0])
)
measurements["GpsWeekNumber"] = np.floor(1e-9 * measurements["tRxGnssNanos"] / WEEKSEC)
measurements["tRxSeconds"] = (
    1e-9 * measurements["tRxGnssNanos"] - WEEKSEC * measurements["GpsWeekNumber"]
)
measurements["tTxSeconds"] = 1e-9 * (
    measurements["ReceivedSvTimeNanos"] + measurements["TimeOffsetNanos"]
)
# Calculate pseudorange in seconds
measurements["prSeconds"] = measurements["tRxSeconds"] - measurements["tTxSeconds"]

# Conver to meters
measurements["PrM"] = LIGHTSPEED * measurements["prSeconds"]
measurements["PrSigmaM"] = (
    LIGHTSPEED * 1e-9 * measurements["ReceivedSvTimeUncertaintyNanos"]
)


input_filepath = os.path.join("data")
manager = EphemerisManager(input_filepath)


def calculate_satellite_position(ephemeris, transmit_time):
    mu = 3.986005e14
    OmegaDot_e = 7.2921151467e-5
    F = -4.442807633e-10
    sv_position = pd.DataFrame()
    sv_position["sv"] = ephemeris.index
    sv_position.set_index("sv", inplace=True)
    sv_position["GPS_time"] = transmit_time - ephemeris["t_oe"]
    A = ephemeris["sqrtA"].pow(2)
    n_0 = np.sqrt(mu / A.pow(3))
    n = n_0 + ephemeris["deltaN"]
    M_k = ephemeris["M_0"] + n * sv_position["GPS_time"]
    E_k = M_k
    err = pd.Series(data=[1] * len(sv_position.index))
    i = 0
    while err.abs().min() > 1e-8 and i < 10:
        new_vals = M_k + ephemeris["e"] * np.sin(E_k)
        err = new_vals - E_k
        E_k = new_vals
        i += 1

    sinE_k = np.sin(E_k)
    cosE_k = np.cos(E_k)
    delT_r = F * ephemeris["e"].pow(ephemeris["sqrtA"]) * sinE_k
    delT_oc = transmit_time - ephemeris["t_oc"]
    sv_position["delT_sv"] = (
        ephemeris["SVclockBias"]
        + ephemeris["SVclockDrift"] * delT_oc
        + ephemeris["SVclockDriftRate"] * delT_oc.pow(2)
    )

    v_k = np.arctan2(
        np.sqrt(1 - ephemeris["e"].pow(2)) * sinE_k, (cosE_k - ephemeris["e"])
    )

    Phi_k = v_k + ephemeris["omega"]

    sin2Phi_k = np.sin(2 * Phi_k)
    cos2Phi_k = np.cos(2 * Phi_k)

    du_k = ephemeris["C_us"] * sin2Phi_k + ephemeris["C_uc"] * cos2Phi_k
    dr_k = ephemeris["C_rs"] * sin2Phi_k + ephemeris["C_rc"] * cos2Phi_k
    di_k = ephemeris["C_is"] * sin2Phi_k + ephemeris["C_ic"] * cos2Phi_k

    u_k = Phi_k + du_k

    r_k = A * (1 - ephemeris["e"] * np.cos(E_k)) + dr_k

    i_k = ephemeris["i_0"] + di_k + ephemeris["IDOT"] * sv_position["GPS_time"]

    x_k_prime = r_k * np.cos(u_k)
    y_k_prime = r_k * np.sin(u_k)

    Omega_k = (
        ephemeris["Omega_0"]
        + (ephemeris["OmegaDot"] - OmegaDot_e) * sv_position["GPS_time"]
        - OmegaDot_e * ephemeris["t_oe"]
    )

    sv_position["x_k"] = x_k_prime * np.cos(Omega_k) - y_k_prime * np.cos(i_k) * np.sin(
        Omega_k
    )
    sv_position["y_k"] = x_k_prime * np.sin(Omega_k) + y_k_prime * np.cos(i_k) * np.cos(
        Omega_k
    )
    sv_position["z_k"] = y_k_prime * np.sin(i_k)
    return sv_position


epoch = 0
num_sats = 0
sv_position_epoch = pd.DataFrame()
measurements.insert(0, "SatX", 0)
measurements.insert(0, "SatY", 0)
measurements.insert(0, "SatZ", 0)

while num_sats < 9:
    one_epoch = measurements.loc[(measurements["Epoch"] == epoch)].drop_duplicates(
        subset="SvName"
    )
    timestamp = one_epoch.iloc[0]["UnixTime"].to_pydatetime(warn=False)
    one_epoch.set_index("SvName", inplace=True)
    num_sats = len(one_epoch.index)
    sats = one_epoch.index.unique().tolist()
    ephemeris = manager.get_ephemeris(timestamp, sats)
    sv_position_epoch_temp = calculate_satellite_position(
        ephemeris, one_epoch["tTxSeconds"]
    )
    if epoch == 0:
        sv_position_epoch = sv_position_epoch_temp
    else:
        frames = [sv_position_epoch, sv_position_epoch_temp]
        sv_position_epoch = pd.concat(frames)
    for index, row in sv_position_epoch_temp.iterrows():
        # print(index[1:], row["x_k"], row["y_k"], row["z_k"])
        measurements.loc[
            ((measurements.Svid == index[1:]) & (measurements.Epoch == epoch)),
            "SatX",
        ] = row["x_k"]
        measurements.loc[
            ((measurements.Svid == index[1:]) & (measurements.Epoch == epoch)),
            "SatY",
        ] = row["y_k"]
        measurements.loc[
            ((measurements.Svid == index[1:]) & (measurements.Epoch == epoch)),
            "SatZ",
        ] = row["z_k"]
        measurements.loc[
            ((measurements.Svid == index[1:]) & (measurements.Epoch == epoch)),
            "GPS_time",
        ] = row["GPS_time"]
    epoch += 1


# print(timestamp)
# print(one_epoch[["UnixTime", "tTxSeconds", "GpsWeekNumber"]])


# Run the function and check out the results:
# print(one_epoch["tTxSeconds"])
sv_position = calculate_satellite_position(ephemeris, one_epoch["tTxSeconds"])
measurements_temp = measurements
gps_time = measurements["GPS_time"]
new_data = pd.DataFrame().assign(
    GPS_time=gps_time,
    SatPRN=measurements["Svid"],
    Pseudo_range=measurements["PseudorangeRateMetersPerSecond"],
    CN0=measurements["Cn0DbHz"],
    SatX=measurements["SatX"],
    SatY=measurements["SatY"],
    SatZ=measurements["SatZ"],
)
measurements.to_csv("measurements.csv")
new_data.to_csv("ST_POS.csv")
# print(sv_position)

# initial guesses of receiver clock bias and position
b0 = 0
x0 = np.array([0, 0, 0])
xs = sv_position[["x_k", "y_k", "z_k"]].to_numpy()
# print(sv_position)
# Apply satellite clock bias to correct the measured pseudorange values
pr = one_epoch["PrM"] + LIGHTSPEED * sv_position["delT_sv"]
pr = pr.to_numpy()


def least_squares(xs, measured_pseudorange, x0, b0):
    dx = 100 * np.ones(3)
    b = b0
    # set up the G matrix with the right dimensions. We will later replace the first 3 columns
    # note that b here is the clock bias in meters equivalent, so the actual clock bias is b/LIGHTSPEED
    G = np.ones((measured_pseudorange.size, 4))
    iterations = 0
    while np.linalg.norm(dx) > 1e-3:
        r = np.linalg.norm(xs - x0, axis=1)
        phat = r + b0
        deltaP = measured_pseudorange - phat
        G[:, 0:3] = -(xs - x0) / r[:, None]
        sol = np.linalg.inv(np.transpose(G) @ G) @ np.transpose(G) @ deltaP
        dx = sol[0:3]
        db = sol[3]
        x0 = x0 + dx
        b0 = b0 + db
    norm_dp = np.linalg.norm(deltaP)
    return x0, b0, norm_dp


x, b, dp = least_squares(xs, pr, x0, b0)
print("rough pos ", navpy.ecef2lla(x))
print(b / LIGHTSPEED)
print(dp)


ecef_list = []
ecef_list_by_time = []
for epoch in measurements["Epoch"].unique():
    one_epoch = measurements.loc[
        (measurements["Epoch"] == epoch) & (measurements["prSeconds"] < 0.1)
    ]
    one_epoch = one_epoch.drop_duplicates(subset="SvName").set_index("SvName")
    if len(one_epoch.index) > 4:
        timestamp = one_epoch.iloc[0]["UnixTime"].to_pydatetime(warn=False)
        sats = one_epoch.index.unique().tolist()
        ephemeris = manager.get_ephemeris(timestamp, sats)
        sv_position = calculate_satellite_position(ephemeris, one_epoch["tTxSeconds"])

        xs = sv_position[["x_k", "y_k", "z_k"]].to_numpy()
        pr = one_epoch["PrM"] + LIGHTSPEED * sv_position["delT_sv"]
        pr = pr.to_numpy()

        x, b, dp = least_squares(xs, pr, x, b)
        ecef_list.append(x)
        ecef_list_by_time.append((x, sv_position["GPS_time"].min()))

# print(ecef_list_by_time)

lla = [navpy.ecef2lla(coord) for (coord, time) in ecef_list_by_time]
ecef = [coord for (coord, time) in ecef_list_by_time]
locations = []

for ecef_coord, time in ecef_list_by_time:
    lla_coord = navpy.ecef2lla(ecef_coord)
    ecef_x = ecef_coord[0]
    ecef_y = ecef_coord[1]
    ecef_z = ecef_coord[2]
    lla_lat = lla_coord[0]
    lla_lon = lla_coord[1]
    lla_alt = lla_coord[2]
    row = (time, ecef_x, ecef_y, ecef_z, lla_lat, lla_lon, lla_alt)
    locations.append(row)

print(locations)
locations_df = pd.DataFrame(
    locations, columns=["GPS_time", "Pos.X", "Pos.Y", "Pos,Z", "Lat", "Lon", "Alt"]
)


final_df = pd.merge(new_data, locations_df, on="GPS_time")
final_df.to_csv("final_df.csv")

import simplekml

kml = simplekml.Kml()
for coord in lla:
    lat, lon, alt = coord
    kml.newpoint(name="", coords=[(lat, lon, alt)])
kml.save("coords.kml")
