import numpy as np
from astropy import units as u
from . import Const

# convert 'year:days of year:seconds of day' to years
def year_day_second(yds):
    year,days,seconds = np.array(yds.split(':'),dtype=np.int)
    if year < 80:
        year = year+100
    return year+(days+seconds/86400)/365.25

def geode2geocen(geode_lat):
    f0 = Const.f0
    tan_geocen_lat = np.tan(geode_lat) * (1-f0)**2
    geocen_lat = np.arctan(tan_geocen_lat).to(u.deg)
    return geocen_lat

def geocen2geode(geocen_lat):
    f0 = Const.f0
    tan_geode_lat = np.tan(geocen_lat) / (1-f0)**2
    geode_lat = np.arctan(tan_geode_lat).to(u.deg)
    return geode_lat    







