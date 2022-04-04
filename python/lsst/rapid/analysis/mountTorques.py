# This file is part of rapid_analysis.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import asyncio
import lsst.rapid.analysis.butlerUtils as bu
from .utils import dayObsIntToString
from astro_metadata_translator import ObservationInfo

try:
    from lsst_efd_client import merge_packed_time_series as mpts
except ImportError:
    pass

GOOD_IMAGE_TYPES = ['OBJECT', 'SKYEXP', 'ENGTEST', 'SCIENCE']


def _getEfdData(client, dataSeries, startTime, endTime):
    """A synchronous warpper for geting the data from the EFD.

    This exists so that the top level functions don't all have to be async def.
    """
    loop = asyncio.get_event_loop()
    return loop.run_until_complete(client.select_time_series(dataSeries, ['*'], startTime.utc, endTime.utc))


def plotMountTracking(dataId, butler, client, figure, saveFilename, logger):
    """Queries EFD for given exposure and checks if there were tracking errors.

    Parameters
    ----------
    dataId : `dict`
        The dataId for quich to plot the mount torques.
    butler : `lsst.daf.butler.Butler`
        The butler to use to retrieve the image metadata.
    client : `lsst_efd_client.Client`
        The EFD client to retrieve the mount torques.
    figure : `matplotlib.Figure`
        A matplotlib figure to re-use. Necessary to pass this in to prevent an
        ever-growing figure count and the ensuing memory leak.
    saveFilename : `str`
        Full path and filename to save the plot to.
    logger : `logging.Logger`
        The logger.

    Returns
    -------
    plotted : `bool`
        True if the dataId was plotted, False if it was skipped.
    """
    # lsst-efd-client is not a required import at the top here, but is
    # implicitly required as a client is passed into this function so is not
    # rechecked here.

    start = time.time()

    expRecord = bu.getExpRecordFromDataId(butler, dataId)
    dayString = dayObsIntToString(expRecord.day_obs)
    seqNumString = str(expRecord.seq_num)
    dataIdString = f"{dayString} - seqNum {seqNumString}"

    imgType = expRecord.observation_type.upper()
    tStart = expRecord.timespan.begin.tai.to_value("isot")
    tEnd = expRecord.timespan.end.tai.to_value("isot")
    elevation = 90 - expRecord.zenith_angle
    exptime = expRecord.exposure_time

    # TODO: DM-33859 remove this once it can be got from the expRecord
    md = butler.get('raw.metadata', dataId, detector=0)
    obsInfo = ObservationInfo(md)
    azimuth = obsInfo.altaz_begin.az.value
    logger.debug(f"dataId={dataIdString}, imgType={imgType}, Times={tStart}, {tEnd}")

    if imgType not in GOOD_IMAGE_TYPES:
        logger.info(f'Skipping image type {imgType} for {dataIdString}')
        return False
    if exptime < 1.99:
        logger.info('Skipping sub 2s expsoure')
        return False

    end = time.time()
    elapsed = end-start
    logger.debug(f"Elapsed time for butler query = {elapsed}")

    start = time.time()
    # Time base in the EFD is still a big mess.  Although these times are in
    # UTC, it is necessary to tell the code they are in TAI. Then it is
    # necessary to tell the merge_packed_time_series to use UTC.
    # After doing all of this, there is still a 2 second offset,
    # which is discussed in JIRA ticket DM-29243, but not understood.

    t_start = Time(tStart, scale='tai')
    t_end = Time(tEnd, scale='tai')
    logger.debug(f"Tstart = {t_start.isot}, Tend = {t_end.isot}")

    mount_position = _getEfdData(client, "lsst.sal.ATMCS.mount_AzEl_Encoders", t_start, t_end)
    nasmyth_position = _getEfdData(client, "lsst.sal.ATMCS.mount_Nasmyth_Encoders", t_start, t_end)
    torques = _getEfdData(client, "lsst.sal.ATMCS.measuredTorque", t_start, t_end)
    logger.debug("Length of time series", len(mount_position))

    az = mpts(mount_position, 'azimuthCalculatedAngle', stride=1)
    el = mpts(mount_position, 'elevationCalculatedAngle', stride=1)
    rot = mpts(nasmyth_position, 'nasmyth2CalculatedAngle', stride=1)
    az_torque_1 = mpts(torques, 'azimuthMotor1Torque', stride=1)
    az_torque_2 = mpts(torques, 'azimuthMotor2Torque', stride=1)
    el_torque = mpts(torques, 'elevationMotorTorque', stride=1)
    rot_torque = mpts(torques, 'nasmyth2MotorTorque', stride=1)

    end = time.time()
    elapsed = end-start
    logger.debug(f"Elapsed time to get the data = {elapsed}")
    start = time.time()

    # Calculate the tracking errors
    az_vals = np.array(az.values[:, 0])
    el_vals = np.array(el.values[:, 0])
    rot_vals = np.array(rot.values[:, 0])
    times = np.array(az.values[:, 1])
    logger.debug("Length of packed time series", len(az_vals))

    # Fit with a linear
    az_fit = np.polyfit(times, az_vals, 1)
    el_fit = np.polyfit(times, el_vals, 1)
    rot_fit = np.polyfit(times, rot_vals, 1)
    az_model = az_fit[0] * times + az_fit[1]
    el_model = el_fit[0] * times + el_fit[1]
    rot_model = rot_fit[0] * times + rot_fit[1]

    # Errors in arcseconds
    az_error = (az_vals - az_model) * 3600
    el_error = (el_vals - el_model) * 3600
    rot_error = (rot_vals - rot_model) * 3600

    # Calculate RMS
    az_rms = np.sqrt(np.mean(az_error * az_error))
    el_rms = np.sqrt(np.mean(el_error * el_error))
    rot_rms = np.sqrt(np.mean(rot_error * rot_error))

    end = time.time()
    elapsed = end-start
    logger.debug(f"Elapsed time for error calculations = {elapsed}")
    start = time.time()

    # Plotting
    figure.clear()
    title = f"Mount Tracking {dataIdString}, Azimuth = {azimuth:.1f}, Elevation = {elevation:.1f}"
    plt.suptitle(title, fontsize=18)
    # Azimuth axis
    plt.subplot(3, 3, 1)
    ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
    ax1.set_title("Azimuth axis", fontsize=16)
    ax1.axvline(az.index[0], color="red", linestyle="--")
    ax1.set_xticks([])
    ax1.set_ylabel("Degrees")
    plt.subplot(3, 3, 4)
    plt.plot(times, az_error, color='red')
    plt.title(f"Azimuth RMS error = {az_rms:.2f} arcseconds")
    plt.ylim(-10.0, 10.0)
    plt.xticks([])
    plt.ylabel("Arcseconds")
    plt.subplot(3, 3, 7)
    ax7 = az_torque_1['azimuthMotor1Torque'].plot(legend=True, color='blue')
    ax7 = az_torque_2['azimuthMotor2Torque'].plot(legend=True, color='green')
    ax7.axvline(az.index[0], color="red", linestyle="--")
    ax7.set_ylabel("Torque (motor current in amps)")

    # Elevation axis
    plt.subplot(3, 3, 2)
    ax2 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
    ax2.set_title("Elevation axis", fontsize=16)
    ax2.axvline(az.index[0], color="red", linestyle="--")
    ax2.set_xticks([])
    plt.subplot(3, 3, 5)
    plt.plot(times, el_error, color='green')
    plt.title(f"Elevation RMS error = {el_rms:.2f} arcseconds")
    plt.ylim(-10.0, 10.0)
    plt.xticks([])
    plt.subplot(3, 3, 8)
    ax8 = el_torque['elevationMotorTorque'].plot(legend=True, color='blue')
    ax8.axvline(az.index[0], color="red", linestyle="--")
    ax8.set_ylabel("Torque (motor current in amps)")

    # Nasmyth2 rotator axis
    plt.subplot(3, 3, 3)
    ax3 = rot['nasmyth2CalculatedAngle'].plot(legend=True, color='blue')
    ax3.set_title("Nasmyth2 axis", fontsize=16)
    ax3.axvline(az.index[0], color="red", linestyle="--")
    ax3.set_xticks([])
    plt.subplot(3, 3, 6)
    plt.plot(times, rot_error, color='blue')
    plt.title(f"Nasmyth RMS error = {rot_rms:.2f} arcseconds")
    plt.ylim(-100.0, 100.0)
    plt.subplot(3, 3, 9)
    ax9 = rot_torque['nasmyth2MotorTorque'].plot(legend=True, color='blue')
    ax9.axvline(az.index[0], color="red", linestyle="--")
    ax9.set_ylabel("Torque (motor current in amps)")
    plt.savefig(saveFilename)

    end = time.time()
    elapsed = end-start
    logger.debug(f"Elapsed time for plots = {elapsed}")

    return True
