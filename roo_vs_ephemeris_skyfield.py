from math import asin, atan2, floor
import georinex as gr
from pathlib import Path
import numpy as np
import pandas as pd
import glob
import os
import numerical_methods as num_meth
import requests
from datetime import datetime, timedelta
from alive_progress import alive_bar
from skyfield.api import N, W, wgs84, load, Distance, utc, Angle
from skyfield.toposlib import ITRSPosition
import pytz

# CONFIG
# Place telescope data in <telescope_data_folder>/QZSS1 or QZSS3
# Program will automatically detect the correct satellite (provided the filename contains the name of the correct satellite)
# Program will also automatically download and choose the correct ephemeris and save it to <ephemeris_folder> 
# Might cause problems if the telescope's data includes 2 different days
# The UTC time, interpolated RA/DEC, telescope RA/DEC and the difference will be saved in csv files under <output_dir>

ephemeris_folder = "qzss_ephemeris/" # Choose location of QZSS ephemeris files
telescope_data_folder = "1_telescope_data/" # Choose where to put the telescope obs csvs
output_dir = "2_compared_data/" # Choose where to save final comparison data
load_sp3 = Path("./load_sp3") # Won't work without this for some reason (for gr.loadsp3() function)

ephemeris_type = "final" # Choose between rapid or final

roo = [-37.680589141, 145.061634327, 155.083] # Location of the ROO

####################################################################

QZS1 = "J01"
QZS2 = "J02"
QZS3 = "J07"
QZS4 = "J03"
PRN5 = "G05"
PRN18 = "G18"

ephemeris_final_url = "https://sys.qzss.go.jp/archives/final-sp3/"
ephemeris_rapid_url = "https://sys.qzss.go.jp/archives/rapid-sp3/"

planets = load('de421.bsp')
earth = planets['earth']

obs_location = earth + wgs84.latlon(roo[0], roo[1], roo[2])

def choose_ephemeris(ephemeris_type, ephemeris_folder, time):
    gps_fepoch = datetime(1980, 1, 6)
    time_since_epoch = time - gps_fepoch
    gps_week = floor(time_since_epoch.total_seconds()/86400 / 7)
    gps_day = floor(time_since_epoch.total_seconds()/86400) % 7
    gps_weekday = str(gps_week) + str(gps_day)
    download_link = None
    file_path = None
    file_name = None

    file_name = None
    year = str(time.year) + "/"
    if ephemeris_type == "rapid":
        file_name = "qzr" + gps_weekday + ".sp3"
        download_link = ephemeris_rapid_url + year + file_name
        file_path = Path(ephemeris_folder + file_name)

    elif ephemeris_type == "final":
        file_name = "qzf" + gps_weekday + ".sp3"
        download_link = ephemeris_final_url + year + file_name
        file_path = Path(ephemeris_folder + file_name)
    else:
        print("Please choose a valid ephemeris!")

    if file_path.is_file():
        print(file_name + " already exists, using")
    else:
        print("Downloading final ephemeris " + download_link + "...")
        ephemeris_file = requests.get(download_link, allow_redirects=True)
        open(ephemeris_folder + file_name, 'wb').write(ephemeris_file.content)
    
    return file_path

def interpolate(x, y, z, satellite_time_utc, telescope_datetime):
    satellite_time_utc_timestamp = np.array([pd.to_datetime([time]).astype(np.int64) / 10**6 for time in satellite_time_utc]).flatten()
    telescope_datetime_timestamp = np.array([pd.to_datetime([str(time)]).astype(np.int64) / 10**6 for time in telescope_datetime]).flatten()
    x_interp = []
    y_interp = []
    z_interp = []
    for time in telescope_datetime_timestamp:
        x_interp.append(num_meth.interp_lagrange(satellite_time_utc_timestamp, x, time, 9)[0])
        y_interp.append(num_meth.interp_lagrange(satellite_time_utc_timestamp, y, time, 9)[0])
        z_interp.append(num_meth.interp_lagrange(satellite_time_utc_timestamp, z, time, 9)[0])
    return x_interp, y_interp, z_interp

def transform(satellite_itrs_x, satellite_itrs_y, satellite_itrs_z, telescope_datetime, obs_location):
    satellite = Distance(km=[satellite_itrs_x, satellite_itrs_y, satellite_itrs_z])
    satellite_itrs = ITRSPosition(satellite) + earth
    timezone = pytz.utc
    obstime_utc = timezone.localize(telescope_datetime)
    ts = load.timescale()
    t = ts.utc(obstime_utc)

    topo_ra, topo_dec, distance = obs_location.at(t).observe(satellite_itrs).apparent().radec()
    
    return topo_ra, topo_dec
    
def extract_ephemeris(satellite, sp3):
    satellite_itrs_x = sp3.sel(sv=satellite).position.sel(ECEF='x').to_numpy()
    satellite_itrs_y = sp3.sel(sv=satellite).position.sel(ECEF='y').to_numpy()
    satellite_itrs_z = sp3.sel(sv=satellite).position.sel(ECEF='z').to_numpy()
    satellite_time_gps = sp3.sel(sv=satellite).position.time.to_numpy()
    satellite_time_utc = []
    for time in satellite_time_gps:
        satellite_time_utc.append(time - np.timedelta64(18000, 'ms'))
    return satellite_itrs_x, satellite_itrs_y, satellite_itrs_z, satellite_time_utc
    
def interpolate_transform(satellite, telescope_datetime, ephemeris_folder):
    path_to_sp3 = choose_ephemeris(ephemeris_type, ephemeris_folder, telescope_datetime[0])
    sp3 = gr.load_sp3(path_to_sp3, load_sp3)
    print(sp3)
    
    satellite_itrs_x, satellite_itrs_y, satellite_itrs_z, satellite_time_utc = extract_ephemeris(satellite, sp3)
        
    x_interp, y_interp, z_interp = interpolate(satellite_itrs_x, satellite_itrs_y, satellite_itrs_z, satellite_time_utc, telescope_datetime)
    
    ephemeris_ra = []
    ephemeris_dec = []
    print("Transforming...")
    with alive_bar(len(x_interp)) as bar:
        for i in range(len(x_interp)):
            ra, dec = transform(x_interp[i], y_interp[i], z_interp[i], telescope_datetime[i], obs_location)
            ephemeris_ra.append(ra)
            ephemeris_dec.append(dec)
            bar()
    return ephemeris_ra, ephemeris_dec

def main(ephemeris_folder, telescope_data_folder, output_dir, offset):
    telescope_obs_files = glob.glob(telescope_data_folder + "*.csv")

    for file in telescope_obs_files:
        #Read telescope observation file
        file_name = os.path.basename(file)
        print("Reading telescope observation file " + file_name + "...")
        telescope_obs_file = np.genfromtxt(file, delimiter = ",", dtype = None, encoding='utf-8')

        #Parse observation file data
        telescope_datetime = np.array([datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%f') + offset for time in telescope_obs_file[1:, 0]])
        telescope_ra = telescope_obs_file[1:, 3].astype(float)
        telescope_dec = telescope_obs_file[1:, 4].astype(float)

        if ("qzs1" in file):
            satellite = QZS1
            ephemeris_ra, ephemeris_dec = interpolate_transform(satellite, telescope_datetime, ephemeris_folder)
        elif ("qzs2" in file):
            satellite = QZS2
            ephemeris_ra, ephemeris_dec = interpolate_transform(satellite, telescope_datetime, ephemeris_folder)
        elif ("qzs3" in file):
            satellite = QZS3
            ephemeris_ra, ephemeris_dec = interpolate_transform(satellite, telescope_datetime, ephemeris_folder)
        elif ("qzs4" in file):
            satellite = QZS4
            ephemeris_ra, ephemeris_dec = interpolate_transform(satellite, telescope_datetime, ephemeris_folder)

        telescope_datetime_timestamp = np.array([(datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%f') + offset).isoformat() for time in telescope_obs_file[1:, 0]])
        ephemeris_ra = [x.arcseconds()/3600 for x in ephemeris_ra]
        ephemeris_dec = [x.arcseconds()/3600 for x in ephemeris_dec]
        ra_difference = np.subtract(ephemeris_ra, telescope_ra) * 3600
        dec_difference = np.subtract(ephemeris_dec, telescope_dec) * 3600
        heading = ["Timestamp", "Telescope RA", "Telescope DEC", "Ephemeris RA", "Ephemeris DEC", "RA Difference", "DEC Difference"]
        final_data = np.stack((telescope_datetime_timestamp, telescope_ra, telescope_dec, ephemeris_ra, ephemeris_dec, ra_difference, dec_difference), axis=1)
        formatted_data = np.vstack((heading, final_data))
        if offset < timedelta(seconds = 0):
            output_file = output_dir + str(int(offset.microseconds / 1000 - 1000)) + "ms_compared_" + file_name
        elif offset == timedelta(seconds = 1):
            output_file = output_dir + str(int(offset.seconds * 1000)) + "ms_compared_" + file_name
        else:
            output_file = output_dir + str(int(offset.microseconds / 1000)) + "ms_compared_" + file_name
        if not(os.path.isdir(output_dir)):
            os.mkdir(output_dir)
        np.savetxt(output_file, formatted_data, delimiter=',', fmt='%s')
        print("Done :D Final data saved in " + output_file + "\n")

offset_datetime = timedelta(milliseconds=int(0))
main(ephemeris_folder, telescope_data_folder, output_dir, offset = offset_datetime)

