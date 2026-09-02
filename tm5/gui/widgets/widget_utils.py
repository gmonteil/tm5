#!/usr/bin/env python

import panel as pn
from pathlib import Path
from pandas import DataFrame
import xarray as xr


def experiment_desc( exp : str ) -> str:
    desc = "!!! description missing !!!"
    match exp:
        case 'default':
            desc = f"standard/default emission scenario"
        case 'edgarflat':
            desc = f"Similar to the default case, but using a flat " \
                f"temporal profile for EDGAR anthropogenic emissions."
        case 'regional':
            desc = f"Similar to the default case, but emissions from " \
                f"wetlands, mineral-soils, and anthropogenic sources " \
                f"over the European domain are taken from dedicated datasets " \
                f"generated in AVENGERS WP2."
        case 'regional_no-agri' | 'regional_anthro-no-agri' | 'regional-anthro-no-agri':
            desc = f"Emissions similar to the regional case, " \
                f"but without emissions from the agriculture sector " \
                f"over the European domain."
        case 'regional_no-fossil' | 'regional_anthro-no-fossil' | 'regional-anthro-no-fossil':
            desc = f"Emissions similar to the regional case, " \
                f"but without emissions from the fossil sector " \
                f"over the European domain."
        case 'regional_no-waste' | 'regional_anthro-no-waste' | 'regional-anthro-no-waste':
            desc = f"Emissions similar to the regional case, " \
                f"but without emissions from the waste sector " \
                f"over the European domain."
        case 'regional_no-anthro-france' | 'regional_anthro-no-france' | 'regional-anthro-no-france':
            desc = f"Emissions similar to the regional case, " \
                f"but without anthropogenic emissions over France."
        case 'regional_no-anthro-netherlands' | 'regional_anthro-no-netherlands' | 'regional-anthro-no-netherlands':
            desc = f"Emissions similar to the regional case, " \
                f"but without anthropogenic emissions over " \
                f"the Netherlands."
        case 'half-oh':
            desc = f"Emissions similar to the default case, " \
                f"but using halved CAMS OH concentrations " \
                f"(which are entering the TM5 chemistry)."
        case 'no-germany':
            desc = "Emissions similar to the default case, " \
                f"but without emissions over domain around Germany " \
                f"(6E-15E,47N-55N)."
        case 'no-gns':
            desc = "Emissions similar to the default case, " \
                f"but without emissions over the innermost zoom domain " \
                f"(0E-18E,42N-58N) covering Germany, Netherlands, and Switzerland."
        case 'no-northamerica':
            desc = "Emissions similar to the default case, " \
                f"but without emissions over Northern America " \
                f"(165W-55W,25N-80N)."
        #-- MVO-20260529:ad-hoc catch for the file Zois had placed onto the exploredata platform!
        case 'mytest-emissions':
            desc = "!!!NOT ACTIVE YET!!! for the future it is foreseen that users can upload " \
                f"their own emission fields."
    return desc


def load_observations_metadata(fname: Path) -> DataFrame:
    """
    Complementary function to "load_observations_data": this one loads a bunch of metadata, for each site:
    - site_name
    - site_code
    - country
    - latitude
    - longitude
    - elevation
    - doi
    - filename
    """
    vars_select = ['time', 'value', 'altitude', 'latitude', 'longitude', 'elevation', 'intake_height']
    ds = xr.open_dataset(fname, decode_timedelta=False)[vars_select]
    return DataFrame({
        'site_name': ds.attrs['site_name'],
        'site_code': ds.attrs['site_code'],
        'country': ds.attrs['site_country'],
        'latitude': ds.attrs['site_latitude'],
        'longitude': ds.attrs['site_longitude'],
        'elevation': ds.attrs['site_elevation'],
        'doi': ds.attrs['obspack_identifier_link'],
        'filename': fname
    }, index=[ds.attrs['site_name']])


def plot_site_info(sites: DataFrame, station: str | None):
    site = sites.loc[station]

    text = pn.pane.Markdown(f"""
    ### {site.site_name}

    - latitude: {site.latitude}
    - longitude: {site.longitude}
    - elevation: {site.elevation}
    - DOI: {site.doi}
    """)

    return pn.Column(
        text,
        sites.hvplot.points(
            x='longitude', y='latitude', geo=True, coastline=True, xlim=(-180, 180), ylim=(-90, 90),
            frame_width=300, hover_cols=['site_name']
        ) *
        sites.loc[[station]].hvplot.points(
            x='longitude', y='latitude', geo=True, coastline=True, xlim=(-180, 180), ylim=(-90, 90), frame_width=300, color='r'
        )
    )
