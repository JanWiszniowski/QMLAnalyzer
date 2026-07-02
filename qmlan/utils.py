"""
Commonly used utils for seismic data processing be the seismic processing in Python packages
--------------------------------------------------------------------------------------------

..
    :copyright:
        Jan Wiszniowski (jwisz@igf.edu.pl)
    :license:
        GNU Lesser General Public License, Version 3
        (https://www.gnu.org/copyleft/lesser.html)
    :version 0.0.1:
        2024-11-07

"""

import math
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.event.base import WaveformStreamID, Comment
from obspy.core.event.magnitude import Magnitude
import obspy.geodetics.base as geo


MinUTCDateTime = UTCDateTime(-4.0e9)
MaxUTCDateTime = UTCDateTime(4.0e9)
_units_dictionary = {"'VEL'": "m/s", "'ACC'": "m/s^2", "'DISP'": "m"}


def get_net_sta(name):
    """
    Function get_net_sta extracts network and station_name codes as strings

    :param name: The trace name. It can be the string or the WaveformStreamID object.
        The text in the string is in the form 'NN.SSS.LL.CCC', where NN is the network code,
        SSS is the station_name code, LL is the location code, and CCC is the channel code.
    :type name: str or ObsPy.WaveformStreamID

    :return: The tuple of the network code the station_name code.
    :rtype: tuple(str, str)

    """
    if isinstance(name, str):
        if '.' not in name:
            return '', name
        names = name.split('.', 3)
        return names[0], names[1]
    elif isinstance(name, WaveformStreamID):
        return name.network_code, name.station_code
    else:
        return '', ''


def get_station_name(name):
    """
    Function get_station_name extracts the station_name name as a string

    :param name: The trace name. It can be the string or the WaveformStreamID object.
        The text in the string is in the form 'NN.SSS.LL.CCC', where NN is the network code,
        SSS is the station_name code, LL is the location code, and CCC is the channel code.
    :type name: str or ObsPy.WaveformStreamID

    :return: The string in the form 'NN.STA', where NN is the network code and SSS is the station_name code.
    :rtype: str

    """
    net, sta = get_net_sta(name)
    return "{}.{}".format(net, sta)


def get_station_id(name):
    """
    Function get_station_id extracts the station_name name as a WaveformStreamID object

    :param name: The trace name. It can be the string or the ObsPy WaveformStreamID object.
        The text in the string is in the form 'NN.SSS.LL.CCC', where NN is the network code,
        SSS is the station_name code, LL is the location code, and CCC is the channel code.
    :type name: str or ObsPy.WaveformStreamID

    :return: The waveform stream object containing only the network code and the station_name code.
    :rtype: ObsPy.WaveformStreamID

    """
    net, sta = get_net_sta(name)
    return WaveformStreamID(network_code=net, station_code=sta)


def get_origin(event):
    """
    Function get_origin extracts the origin from the event.
    If preferred_origin_id of the event is set it return the preferred origin.
    Otherwise, it returns the first origin from the list.
    The function is intended to extract the event origin unconditionally and non-interactively.
    Therefore, if preferred_origin_id is not set and there are multiple origins, the returned origin may be random

    :param event: The seismic event object
    :type event: ObsPy.Event

    :return: The origin (event location) object or None if none origin is defined for the event.
    :rtype: ObsPy.Origin

    """
    if not event.origins:
        return None
    if event.preferred_origin_id is not None:
        return event.preferred_origin_id.get_referred_object()
    return event.origins[0]


def get_focal_mechanism(event, inversion_type=None):
    """
    Function get_focal_mechanism extracts the focal mechanism from the event.
    If preferred_focal_mechanism_id of the event is set it return the preferred focal mechanism.
    Otherwise, it returns the first focal mechanism from the list.
    The function is intended to extract the focal mechanism unconditionally and non-interactively.
    Therefore, if preferred_focal_mechanism_id is not set and there are multiple focal mechanisms,
    the returned focal mechanism may be random.

    :param event: The seismic event object
    :type event: ObsPy.Event
    :param inversion_type: The name of tensor inversion type.
        It must belong to the QuakeML MTInversionType category:
        ``'general'``, ``'zero trace'``, ``'double couple'``, or None.
    :type inversion_type: (str)

    :return: The focal mechanism object or None if none focal_mechanism is defined for the event
        or the focal_mechanism with the defined inversion type does not exist.
    :rtype: ObsPy.FocalMechanism

    """
    if not event.focal_mechanisms:
        return None
    if event.preferred_focal_mechanism_id is not None:
        if inversion_type is None:
            return event.preferred_focal_mechanism_id.get_referred_object()
        fm = event.preferred_focal_mechanism_id.get_referred_object()
        if fm.moment_tensor and fm.moment_tensor.inversion_type == inversion_type:
            return event.preferred_focal_mechanism_id.get_referred_object()
    if inversion_type is None:
        return event.focal_mechanisms[0]
    else:
        for fm in event.focal_mechanisms:
            if fm.moment_tensor and fm.moment_tensor.inversion_type == inversion_type:
                return fm
    return None


def get_hypocentral_distance(origin, station_inventory):
    """
    Function get_hypocentral_distance computes the local hypocentral distance in meters
    from origin coordinates to station_name coordinates.
    The calculations do not take into account the curvature of the earth.

    :param origin: The ObsPy Origin object
    :type origin: ObsPy.Origin
    :param station_inventory: The station inventory object
    :type station_inventory: ObsPy.Inventory

    :return: The hypocentral distance in meters and epicentral distance in degrees
    :rtype: tuple(float, float)

    """
    delta = geo.locations2degrees(origin.latitude, origin.longitude,
                                  station_inventory.latitude, station_inventory.longitude)
    epi_distance = geo.degrees2kilometers(delta) * 1000.0
    hypo_distance = math.sqrt(epi_distance ** 2 + (origin.depth + station_inventory.elevation) ** 2)
    return hypo_distance, delta


def this_magnitude(magnitude_object, **kwargs):
    magnitude_type = kwargs.get('magnitude_type')
    method_id = kwargs.get('method_id')
    if magnitude_type is not None and magnitude_object.magnitude_type != magnitude_type:
        return False
    if method_id is not None and magnitude_object.method_id != method_id:
        return False
    return True


def this_sta_magnitude(magnitude_object, **kwargs):
    magnitude_type = kwargs.get('magnitude_type')
    method_id = kwargs.get('method_id')
    if magnitude_type is not None and magnitude_object.station_magnitude_type != magnitude_type:
        return False
    if method_id is None:
        return True
    if magnitude_object.method_id.text == method_id:
        return True
    if method_id == 'none':
        return magnitude_object.method_id is None
    return False


def get_magnitude(event, **kwarg):
    """
    Function get_magnitude extracts the magnitude of the event.
    If you want to extract a specific magnitude type you can define it as magnitude_type,
    e.g. ``get_magnitude(event, magnitude_type='Mw')``, otherwise, any magnitude will be extracted.
    If the preferred_magnitude_id of the event is set it returns the preferred magnitude.
    Otherwise, it returns the first magnitude from the list.
    The function is intended to extract the magnitude unconditionally and non-interactively.
    Therefore, if preferred_magnitude_id is not set and there are multiple magnitudes,
    the returned magnitude may be random.

    If event magnitude does not exist, but station magnitudes exist, the new magnitude is computed
    as the mean value of all station magnitudes or all station magnitudes of given type.

    :param event: The seismic event object
    :type event: ObsPy.Event
    :param magnitude_type:  (optional)
        Describes the type of magnitude. Proposed values are:
        * local magnitude (``'ML'``),
        * moment magnitude (``'Mw'``),
        * energy (``'Energy'``),
        * etc.
    :type magnitude_type:  str
    :return: The magnitude object or None if the function cannot find or create the magnitude.
        If only station_name magnitudes exist, the new ObsPy Magnitude object is created,
        but it is not appended to the event
    :rtype: ObsPy.Magnitude

    """
    if event.preferred_magnitude_id is not None:
        magnitude_object = event.preferred_magnitude_id.get_referred_object()
        if this_magnitude(magnitude_object, **kwarg):
            return magnitude_object
    if event.magnitudes:
        for m in event.magnitudes:
            if this_magnitude(m, **kwarg):
                return m
    if event.station_magnitudes:
        no_magnitudes = 0
        magnitude = 0.0
        for m in event.station_magnitudes:
            if this_sta_magnitude(m, **kwarg):
                magnitude += m.mag
                no_magnitudes += 1
        if no_magnitudes:
            magnitude /= no_magnitudes
            magnitude_object = Magnitude(mag=magnitude, **kwarg)
            magnitude_object.comments.append(Comment(text=f"Mean of {no_magnitudes} station_name magnitudes"))
            return magnitude_object
    return None


def time_ceil(time, step):
    """
    Returns the time rounded-up to the specified accuracy.

    :param time: The time object
    :type time: ObsPy.UTCDateTime
    :param step: The accuracy units in seconds
    :type step: float

    :return: The new rounded-up time object
    :rtype: ObsPy.UTCDateTime

    Example::

        >> from obspy.core.utcdatetime import UTCDateTime
        >> from core.utils import time_ceil
        >> time = UTCDateTime(2024, 1, 3, 8, 28, 33, 245678)
        >> time_ceil(time,1.0)
        >> UTCDateTime(2024, 1, 3, 8, 28, 34)
        >> time_ceil(time,60.0)
        >> UTCDateTime(2024, 1, 3, 8, 29)
        >> time_ceil(time,0.1)
        UTCDateTime(2024, 1, 3, 8, 28, 33, 300000)
        >> time_ceil(time,0.01)
        >> UTCDateTime(2024, 1, 3, 8, 28, 33, 250000)
        >> time_ceil(time,0.001)
        >> UTCDateTime(2024, 1, 3, 8, 28, 33, 246000)

    """
    time_int = math.ceil(time.timestamp / step)
    return UTCDateTime(time_int * step)


def time_floor(time, step):
    r"""
    Returns the time rounded-down to the specified accuracy.

    :param time: The time object
    :type time: ObsPy.UTCDateTime
    :param step: The accuracy units in seconds
    :type step: float

    :return: The new rounded-down time object
    :rtype: ObsPy.UTCDateTime

    Example::

        >> from obspy.core.utcdatetime import UTCDateTime
        >> from utils import time_floor
        >> time = UTCDateTime(2024, 1, 3, 8, 28, 33, 245678)
        >> time_floor(time,0.001)
        UTCDateTime(2024, 1, 3, 8, 28, 33, 245000)
        >> time_floor(time,0.01)
        UTCDateTime(2024, 1, 3, 8, 28, 33, 240000)
        >> time_floor(time,0.1)
        UTCDateTime(2024, 1, 3, 8, 28, 33, 200000)
        >> time_floor(time,1.0)
        UTCDateTime(2024, 1, 3, 8, 28, 33)
        >> time_floor(time,60.0)
        UTCDateTime(2024, 1, 3, 8, 28)

    """
    time_int = math.floor(time.timestamp / step)
    return UTCDateTime(time_int * step)


def time_ceil_dist(time, step):
    """
    Returns seconds from the time to the time rounded up to the specified accuracy.

    :param time: The time object
    :type time: ObsPy.UTCDateTime
    :param step: The accuracy units in seconds
    :type step: float

    :return: The period in seconds to the rounded-up time
    :rtype: float

    Example::

        >> from obspy.core.utcdatetime import UTCDateTime
        >> time = UTCDateTime(2024, 1, 3, 8, 28, 33, 245678)
        >> time_ceil_dist(time,0.1)
        0.054322
        >> time_ceil_dist(time,1.0)
        0.754322

    """
    return time_ceil(time, step) - time


def time_floor_dist(time, step):
    """
    Returns seconds from the time to the time rounded up to the specified accuracy.

    :param time: The time object
    :type time: ObsPy.UTCDateTime
    :param step: The accuracy units in seconds
    :type step: float

    :return: The period in seconds to the rounded-down time
    :rtype: float

    Example::

        >> from obspy.core.utcdatetime import UTCDateTime
        >> time = UTCDateTime(2024, 1, 3, 8, 28, 33, 245678)
        >> time_floor_dist(time,0.1)
        0.045678
        >> time_floor_dist(time,1.0)
        0.245678

    """
    return time - time_floor(time, step)


def get_units(trace):
    """
    Return the signal units of the trace

    :param trace: The trace object
    :type trace: ObsPy.Trace

    :return: The string with units: 'm/s', 'm/s^2', or 'm', if the response was removed,
        when in the processing_parameters is the remove_response process defined,
        or 'counts' otherwise
    :rtype: str
    """
    if 'processing' in trace.meta:
        for line in trace.meta.processing:
            pos1 = line.find(': ')
            pos2 = line.find('(')
            processing_name = line[pos1 + 2:pos2]
            processing_parameters = dict((a.strip(), b.strip())
                                         for a, b in (element.split('=')
                                                      for element in line[pos2 + 2:-1].split('::')))
            if processing_name == 'remove_response':
                return _units_dictionary[processing_parameters['output']]
    return 'counts'


class ProcessTrace(object):
    """
    The base class of the trace processing. Implementations of objects of classes derived from the ProcessTrace
    do some processing on traces defined in the derived classes initialization

    :param trace:
        The processed trace
    :type trace: ObsPy.Trace
    :param begin_time:
        It limits the period, where a process is performed.
        If begin_time is not defined or it is earlier than the beginning of the trace,
        the process is performed from the beginning of the trace
    :type begin_time: ObsPy.UTCDateTime
    :param end_time:
        It limits the period, where a process is performed.
        If end_time is not defined or it is later than the end of the trace,
        the process is performed to the end of the trace
    :type end_time: ObsPy.UTCDateTime

    """
    def __init__(self, trace, begin_time=None, end_time=None):
        """
        The initialization sets up the data subset for processing and calculate the beginning and end time
        of selected data subset

        :param trace:
            The processed trace
        :param begin_time:
            It limits the period, where a process is performed.
            If begin_time is not defined or it is earlier than the beginning of the trace
            the process is performed from the beginning of the trace
        :param end_time:
            It limits the period, where a process is performed.
            If end_time is not defined or it is later than the end of the trace
            the process is performed to the end of the trace
        """
        begin_idx = 0
        end_idx = len(trace.data)
        self.data = trace.data[begin_idx:end_idx]
        if begin_time:
            idx = math.floor((begin_time - trace.meta.starttime) / trace.meta.delta)
            if idx >= end_idx:
                self.data = None
                return
            elif idx > 0:
                begin_idx = idx
        if end_time:
            idx = math.ceil((end_time - trace.meta.starttime) / trace.meta.delta) + 1
            if idx < 0:
                self.data = None
                return
            elif idx < end_idx:
                end_idx = idx
        self.data = trace.data[begin_idx:end_idx]
        self.start_time = trace.meta.starttime + trace.meta.delta * begin_idx
        self.end_time = trace.meta.starttime + trace.meta.delta * end_idx


class IndexTrace(object):
    """
    Class for operating directly on time limited part of trace data

    :param trace:
        The processed trace
    :type trace: ObsPy.Trace
    :param begin_time:
        It limits the period, where a process is performed.
        If begin_time is not defined or it is earlier than the beginning of the trace
        the process is performed from the beginning of the trace
    :type begin_time: ObsPy.UTCDateTime
    :param end_time:
        It limits the period, where a process is performed.
        If end_time is not defined or it is later than the end of the trace,
        the process is performed to the end of the trace
    :type end_time: ObsPy.UTCDateTime

    **Class variables**

    :start_time: The time of the first data sample index
    :end_time: The time of the next sample after the last data sample.
        It differs from the ObsPy trace end_time, which points to the last sample of the trace
    :begin_idx: The first data sample index
    :end_idx: The last data sample index + 1

    Example::

        >> from utils import IndexTrace
        >> from obspy.core.utcdatetime import UTCDateTime
        >> t1 = UTCDateTime(2024, 1, 3, 8, 28, 00)
        >> t2 = UTCDateTime(2024, 1, 3, 8, 29, 00)
        >> st = read('test.msd')
        >> indexes = IndexTrace(st[1], begin_time=t1, end_time=t2
        >> for idx in range(indexes.begin_idx, indexes.end_idx):
        ... pass

    """
    def __init__(self, trace, begin_time=None, end_time=None):
        """
        The initialization of the IndexTrace class bases on the trace and time limits

        :param trace:
            The processed trace
        :type trace: ObsPy.Trace
        :param begin_time:
            It limits the period, where a process is performed.
            If begin_time is not defined or it is earlier than the beginning of the trace
            the process is performed from the beginning of the trace
        :type begin_time: ObsPy.UTCDateTime
        :param end_time:
            It limits the period, where a process is performed.
            If end_time is not defined or it is later than the end of the trace,
            the process is performed to the end of the trace
        :type end_time: ObsPy.UTCDateTime

        """
        begin_idx = 0
        end_idx = len(trace.data)
        if begin_time:
            idx = math.floor((begin_time - trace.meta.starttime) / trace.meta.delta)
            if idx >= end_idx:
                self.begin_idx = None
                return
            elif idx > 0:
                begin_idx = idx
        if end_time:
            idx = math.ceil((end_time - trace.meta.starttime) / trace.meta.delta) + 1
            if idx < 0:
                self.begin_idx = None
                return
            elif idx < end_idx:
                end_idx = idx
        self.start_time = trace.meta.starttime + trace.meta.delta * begin_idx
        self.end_time = trace.meta.starttime + trace.meta.delta * end_idx
        self.begin_idx = begin_idx
        self.end_idx = end_idx


class ExtremeTraceValues(ProcessTrace):
    """
    Class that assess the extreme trace values: maximum, minimum, and absolute maximum value

    :param trace:
        The processed trace
    :type trace: ObsPy.Trace
    :param begin_time:
        It limits the period, where a process is performed.
        If begin_time is not defined or it is earlier than the beginning of the trace
        the process is performed from the beginning of the trace
    :type begin_time: ObsPy.UTCDateTime
    :param end_time:
        It limits the period, where a process is performed.
        If end_time is not defined or it is later than the end of the trace,
        the process is performed to the end of the trace
    :type end_time: ObsPy.UTCDateTime

    **Class variables**

    :data: The optionally cut to time limits data. The data are not a new array but subarray of the Trace data
    :start_time: The time of the first data sample.
    :end_time: The time of the next sample after the last data sample.
        It differs from the ObsPy trace end_time, which points to the last sample of the trace
    :max_value: Maximum data value
    :max_value: Minimum data value
    :abs_max: Absolute maximum value.
        abs_max = max(abs(min_value), abs(max_value))

    """
    def __init__(self, trace, begin_time=None, end_time=None):
        """
        The initialization of the ExtremeTraceValues class bases on the trace and time limits

        :param trace:
            The processed trace
        :type trace: ObsPy.Trace
        :param begin_time:
            It limits the period, where a process is performed.
            If begin_time is not defined or it is earlier than the beginning of the trace
            the process is performed from the beginning of the trace
        :type begin_time: ObsPy.UTCDateTime
        :param end_time:
            It limits the period, where a process is performed.
            If end_time is not defined or it is later than the end of the trace,
            the process is performed to the end of the trace
        :type end_time: ObsPy.UTCDateTime
        """
        ProcessTrace.__init__(self, trace, begin_time, end_time)
        self.min_value = min(self.data)
        self.max_value = max(self.data)
        self.abs_max = max(abs(self.min_value), abs(self.max_value))
