# %%
import pandas as pd
import numpy as np
import geopandas
import plotly.express as px
from dash import Dash, html, dcc, Input, Output, no_update, ctx

# %%
df = pd.read_pickle("GEO_dataframe.pkl")
gdf = geopandas.GeoDataFrame(
    df, geometry=geopandas.points_from_xy(df.lon, df.lat), crs="EPSG:4979")

df_ecmwf = pd.read_pickle("GEO_dataframe_ECMWF.pkl")
gdf_ecmwf = geopandas.GeoDataFrame(
    df_ecmwf, geometry=geopandas.points_from_xy(df.lon, df.lat), crs="EPSG:4979")


# %%
valid_stations = {'VZMX': ['NO_TS_MO_6201083',
              'NO_TS_MO_Fjaltring',
              'NO_TS_MO_Fanoebugt',
              'NO_TS_MO_6600366',
              'NO_TS_MO_6600379'],
             'VTM02': ['NO_TS_MO_6201083',
              'NO_TS_MO_6600367',
              'NO_TS_MO_Fjaltring',
              'NO_TS_MO_Fanoebugt',
              'NO_TS_MO_6600366',
              'NO_TS_MO_6600379'],
             'VDIR': ['NO_TS_MO_Fjaltring', 'NO_TS_MO_Fanoebugt'],
             'VAVH': ['NO_TS_MO_6201083',
              'NO_TS_MO_Fjaltring',
              'NO_TS_MO_Fanoebugt',
              'NO_TS_MO_6600379'],
             'VTPK': ['BO_TS_MO_BothnianSea',
              'BO_TS_TG_Tallinnamadal',
              'BO_TS_MO_Uto',
              'NO_TS_MO_6201083',
              'BO_TS_MO_HelsinkiSuomenlinna',
              'NO_TS_MO_6202110',
              'BO_TS_MO_Knollsgrund',
              'NO_TS_MO_6600367',
              'NO_TS_MO_Fjaltring',
              'BO_TS_TG_Soru',
              'BO_TS_MO_BrofjordenWR',
              'BO_TS_MO_FinngrundetWR',
              'BO_TS_TG_Kuivastu',
              'BO_TS_MO_VaderoarnaWR',
              'NO_TS_MO_Fanoebugt',
              'BO_TS_MO_HelsinkiBuoy',
              'NO_TS_MO_6600366',
              'BO_TS_MO_NorthernBaltic',
              'NO_TS_MO_1000044',
              'NO_TS_MO_6600379',
              'BO_TS_MO_HuvudskarOst'],
             'VPED': ['BO_TS_MO_BothnianSea',
              'BO_TS_MO_Uto',
              'NO_TS_MO_6201083',
              'BO_TS_MO_HelsinkiSuomenlinna',
              'NO_TS_MO_6202110',
              'BO_TS_MO_Knollsgrund',
              'NO_TS_MO_6600367',
              'BO_TS_MO_BrofjordenWR',
              'BO_TS_MO_FinngrundetWR',
              'BO_TS_MO_VaderoarnaWR',
              'BO_TS_MO_HelsinkiBuoy',
              'NO_TS_MO_6600366',
              'BO_TS_MO_NorthernBaltic',
              'NO_TS_MO_1000044',
              'NO_TS_MO_6600379'],
             'VMDR': ['NO_TS_MO_6201083',
              'NO_TS_MO_6600379',
              'BO_TS_MO_HuvudskarOst'],
             'VHM0': ['BO_TS_MO_BothnianSea',
              'BO_TS_TG_Tallinnamadal',
              'BO_TS_MO_Uto',
              'NO_TS_MO_6201083',
              'BO_TS_MO_HelsinkiSuomenlinna',
              'NO_TS_MO_6202110',
              'BO_TS_MO_Knollsgrund',
              'NO_TS_MO_6600367',
              'NO_TS_MO_Fjaltring',
              'BO_TS_TG_Soru',
              'BO_TS_MO_BrofjordenWR',
              'BO_TS_MO_FinngrundetWR',
              'BO_TS_TG_Kuivastu',
              'BO_TS_MO_VaderoarnaWR',
              'NO_TS_MO_Fanoebugt',
              'BO_TS_MO_HelsinkiBuoy',
              'NO_TS_MO_6600366',
              'BO_TS_MO_NorthernBaltic',
              'NO_TS_MO_1000044',
              'NO_TS_MO_6600379',
              'BO_TS_MO_HuvudskarOst'],
             'VEMH': ['BO_TS_TG_Tallinnamadal',
              'BO_TS_MO_Knollsgrund',
              'BO_TS_TG_Soru',
              'BO_TS_MO_BrofjordenWR',
              'BO_TS_MO_FinngrundetWR',
              'BO_TS_TG_Kuivastu',
              'BO_TS_MO_VaderoarnaWR',
              'BO_TS_MO_HuvudskarOst']}

# %%
copernicus_to_ecmwf = {'VZMX': {'long_name': 'Maximum zero crossing wave height (Hmax)',
  'unit': 'm',
  'CF standard_name': 'sea_surface_wave_maximum_height',
  'Copernicus Wave': 'VCMX',
  'ECMWF Wave': 'hmax'},
 'VTZM': {'long_name': 'Period of the highest wave (Thmax)',
  'unit': 's',
  'CF standard_name': 'sea_surface_wave_period_of_highest_wave',
  'Copernicus Wave': '',
  'ECMWF Wave': 'tmax'},
 'VTM02': {'long_name': 'Spectral moments (0,2) wave period (Tm02)',
  'unit': 's',
  'CF standard_name': 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
  'Copernicus Wave': 'VTM02',
  'ECMWF Wave': 'mp2'},
 'VT110': {'long_name': 'Average period highest 1/10 wave (T1/10)',
  'unit': 's',
  'CF standard_name': 'sea_surface_wave_mean_period_of_highest_tenth',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'VHZA': {'long_name': 'Average zero crossing wave height (Hzm)',
  'unit': 'm',
  'CF standard_name': 'sea_surface_wave_mean_height',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'VH110': {'long_name': 'Average height highest 1/10 wave (H1/10)',
  'unit': 'm',
  'CF standard_name': 'sea_surface_wave_mean_height_of_highest_tenth',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'VDIR': {'long_name': 'Wave direction rel. true north',
  'unit': 'degree',
  'CF standard_name': 'sea_surface_wave_from_direction',
  'Copernicus Wave': 'VMDR',
  'ECMWF Wave': 'mwd'},
 'VAVT': {'long_name': 'Average period highest 1/3 wave (T1/3)',
  'unit': 's',
  'CF standard_name': 'sea_surface_wave_significant_period',
  'Copernicus Wave': '',
  'ECMWF Wave': 'mwp'},
 'VAVH': {'long_name': 'Average height highest 1/3 wave (H1/3)',
  'unit': 'm',
  'CF standard_name': 'sea_surface_wave_significant_height',
  'Copernicus Wave': 'VHM0',
  'ECMWF Wave': 'swh'},
 'SVEL': {'long_name': 'Sound velocity',
  'unit': 'm s-1',
  'CF standard_name': 'speed_of_sound_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'SINC': {'long_name': 'Shortwave/solar incoming radiation',
  'unit': 'W m-2',
  'CF standard_name': 'surface_downwelling_shortwave_flux_in_air',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'RVFL': {'long_name': 'River flow rate',
  'unit': 'm3 s-1',
  'CF standard_name': 'water_volume_transport_into_sea_water_from_rivers',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'RDVA': {'long_name': 'Radial sea water velocity away from instrument',
  'unit': 'm s-1',
  'CF standard_name': 'radial_sea_water_velocity_away_from _instrument',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'OSAT': {'long_name': 'Oxygen saturation',
  'unit': '%',
  'CF standard_name': 'fractional_saturation_of_oxygen_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'NSCT': {'long_name': 'South-north current component',
  'unit': 'm s-1',
  'CF standard_name': 'northward_sea_water_velocity',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'LGHT': {'long_name': 'Immerged incoming photosynthetic active radiation',
  'unit': 'µmol m-2 s-1',
  'CF standard_name': 'downwelling_photosynthetic_photon_flux_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'EWCT': {'long_name': 'West-east current component',
  'unit': 'm s-1',
  'CF standard_name': 'eastward_sea_water_velocity',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'DRVA': {'long_name': 'Direction of radial vector away from instrument',
  'unit': 'degree_true',
  'CF standard_name': 'direction_of_radial_vector_away_from_instrument',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'CNDC': {'long_name': 'Total chlorophyll',
  'unit': 'mg m-3',
  'CF standard_name': 'mass_concentration_of_chlorophyll_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'CHLT': {'long_name': 'Total chlorophyll',
  'unit': 'mg m-3',
  'CF standard_name': 'mass_concentration_of_chlorophyll_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'CDOM': {'long_name': 'Cdom',
  'unit': 1e-09,
  'CF standard_name': 'concentration_of_colored_dissolved_organic_matter_in_sea_water_expressed_as_equivalent_mass_fraction_of_quinine_sulfate_dihydrate',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'WSPD': {'long_name': 'Horizontal wind speed',
  'unit': 'm s-1',
  'CF standard_name': 'wind_speed',
  'Copernicus Wave': '',
  'ECMWF Wave': 'wind'},
 'WDIR': {'long_name': 'Wind from direction relative true north',
  'unit': 'degree',
  'CF standard_name': 'wind_from_direction',
  'Copernicus Wave': '',
  'ECMWF Wave': 'dwi'},
 'VTZA': {'long_name': 'Average zero crossing wave period (Tz)',
  'unit': 's',
  'CF standard_name': 'sea_surface_wave_mean_period',
  'Copernicus Wave': '',
  'ECMWF Wave': 'mp2'},
 'VTPK': {'long_name': 'Wave period at spectral peak / peak period (Tp)',
  'unit': 's',
  'CF standard_name': 'sea_surface_wave_period_at_variance_spectral_density_maximum',
  'Copernicus Wave': 'VTPK',
  'ECMWF Wave': 'pp1d'},
 'VPSP': {'long_name': 'Wave directional spreading at spectral peak',
  'unit': 'degree',
  'CF standard_name': 'sea_surface_wave_directional_spread_at_variance_spectral_density_maximum',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'VPED': {'long_name': 'Wave principal direction at spectral peak',
  'unit': 'degree',
  'CF standard_name': 'sea_surface_wave_from_direction_at_variance_spectral_density_maximum',
  'Copernicus Wave': 'VPED',
  'ECMWF Wave': 'mwd'},
 'VMDR': {'long_name': 'Mean wave direction from (Mdir)',
  'unit': 'degree',
  'CF standard_name': 'sea_surface_wave_from_direction',
  'Copernicus Wave': 'VMDR',
  'ECMWF Wave': 'mwd'},
 'VHM0': {'long_name': 'Spectral significant wave height (Hm0)',
  'unit': 'm',
  'CF standard_name': 'sea_surface_wave_significant_height',
  'Copernicus Wave': 'VHM0',
  'ECMWF Wave': 'swh'},
 'VEMH': {'long_name': 'Estimated maximum wave height',
  'unit': 'm',
  'CF standard_name': 'sea_surface_wave_maximum_height',
  'Copernicus Wave': 'VCMX',
  'ECMWF Wave': 'hmax'},
 'TUR4': {'long_name': 'Turbidity',
  'unit': 1,
  'CF standard_name': 'sea_water_turbidity',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'TEMP': {'long_name': 'Sea temperature',
  'unit': 'degrees_C',
  'CF standard_name': 'sea_water_temperature',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'SLEV': {'long_name': 'Water surface height above a specific datum',
  'unit': 'm',
  'CF standard_name': 'water_surface_height_above_reference_datum',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'RELH': {'long_name': 'Relative humidity',
  'unit': '%',
  'CF standard_name': 'relative_humidity',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'PSAL': {'long_name': 'Practical salinity',
  'unit': '0.001',
  'CF standard_name': 'sea_water_practical_salinity',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'PRES': {'long_name': 'Sea pressure',
  'unit': 'dbar',
  'CF standard_name': 'sea_water_pressure',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'HCSP': {'long_name': 'Horizontal current speed',
  'unit': 'm s-1',
  'CF standard_name': 'sea_water_speed',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'HCDT': {'long_name': 'Current to direction relative true north',
  'unit': 'degree',
  'CF standard_name': 'direction_of_sea_water_velocity',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'GSPD': {'long_name': 'Gust wind speed',
  'unit': 'm s-1',
  'CF standard_name': 'wind_speed_of_gust',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'FLU2': {'long_name': 'Chlorophyll-a fluorescence',
  'unit': 'mg m-3',
  'CF standard_name': 'mass_concentration_of_chlorophyll_a_fluorescence_in_sea_water (1)',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'DRYT': {'long_name': 'Air temperature in dry bulb',
  'unit': 'degrees_C',
  'CF standard_name': 'air_temperature',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 'DOXY': {'long_name': 'Dissolved oxygen',
  'unit': 'mmol m-3',
  'CF standard_name': 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'DOX2': {'long_name': 'Dissolved oxygen',
  'unit': 'µmol kg-1',
  'CF standard_name': 'moles_of_oxygen_per_unit_mass_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'DOX1': {'long_name': 'Dissolved oxygen',
  'unit': 'ml l-1',
  'CF standard_name': 'volume_fraction_of_oxygen_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'CPHL': {'long_name': 'Chlorophyll-a',
  'unit': 'mg m-3',
  'CF standard_name': 'mass_concentration_of_chlorophyll_a_in_sea_water',
  'Copernicus Wave': '',
  'ECMWF Wave': ''},
 'ATMP': {'long_name': 'Atmospheric pressure at altitude',
  'unit': 'hPa',
  'CF standard_name': 'air_pressure',
  'Copernicus Wave': '',
  'ECMWF Wave': '_'},
 '': {'long_name': '',
  'unit': '',
  'CF standard_name': '',
  'Copernicus Wave': 'VSDY',
  'ECMWF Wave': 'vst'}}

# %%
parameters = ['VZMX','VTM02','VDIR','VAVH','VTPK','VPED','VMDR','VHM0','VEMH']

# %%
options = []
for par in parameters:
    d = {"label": f'{copernicus_to_ecmwf[par]['long_name']} ({par})',
         "value": par}
    options.append(d)
    


# ── Dash app layout -------------------------------------------------
app = Dash(__name__)

#   | List | Map | Graph |
app.layout = html.Div(
    style={"display": "flex", "gap": "1rem"},
    children=[
        dcc.Dropdown(
            id="param-dropdown",
            options=options,
            value='default',
            clearable=False,
            style={"width": "400px"}
        ),
        dcc.Graph(id="map", style={"width": "800px"}),
        # html.Div(
        #     style={"flex": "1 1 50%"},
        #     children=[
        #         dcc.Graph(id="station-graph", style={"height": "400px"})
        #     ],
        # )
        html.Div(id='graph-container')
        
    ],
)

# ── Callback: Select parameter -> draw points on map -----------------------------
@app.callback(
    Output("map", "figure"),
    Input("param-dropdown", "value"),
    prevent_initial_call=True
)
def update_map(param):
    if param is None:
        return px.scatter_map(lat=[], lon=[],zoom=3.5,center={"lat": 60, "lon": 20}).update_layout(margin=dict(l=0, r=0, t=0, b=0))

    station_lst = [name[9:] for name in valid_stations[param]]
    gdf_filtr = gdf[gdf['names'].isin(station_lst)]
    map_fig = px.scatter_map(
        gdf_filtr,
        lat=gdf_filtr.geometry.y,
        lon=gdf_filtr.geometry.x,
        hover_name='names',
        zoom=3.5,
        center={"lat": 60, "lon": 20}
    )
    map_fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    
    return map_fig


# ── Callback: Select point -> draw grap -----------------------------
@app.callback(
    # Output("station-graph", "figure"),
    Output("graph-container", "children"),

    Input("map", "clickData"),
    Input("param-dropdown", "value"),
    prevent_initial_call=True
)
def update_graph(click,param):
    triggered = ctx.triggered_id
    if triggered == "param-dropdown":
        return [dcc.Graph(figure={}, style={"height": "400px"})]
    if triggered == "map":
        if click is None or param is None:
            return no_update
        station_id = click["points"][0]["hovertext"]
        row = df.loc[df["names"] == station_id].iloc[0]
        row_ecmwf = df_ecmwf.loc[df_ecmwf['names'] == station_id].iloc[0]
        
        ef = row_ecmwf['data_ecmwf'][param].values
        cf = row['data_copernicus'][param].values

        if np.isnan(ef.sum()):

            fig = px.line([np.nan]*12, title=f'No data for {row['names']} station')
            graphs = [dcc.Graph(figure=fig, style={"height": "400px"})]
            return graphs
        elif np.isnan(cf.sum()):
            ts_ecmwf = pd.DataFrame({
                "lead_time": range(0,len(ef)),
                param: ef,
                'model': 'ECMWF'
            })
            
            fig = px.line(ts_ecmwf, x="lead_time", y=param,
                        title=f"MAE Time‑series for {row['names']} ({copernicus_to_ecmwf[param]['unit']})",
                        line_dash = 'model',
                        color='model')
            
            graphs = [dcc.Graph(figure=fig, style={"height": "400px"})]
            return graphs
        else:
            
            ts_ecmwf = pd.DataFrame({
                "lead_time": range(0,len(ef)),
                param: ef,
                'model': 'ECMWF'
            })
            fig_ecmwf = px.line(ts_ecmwf, x="lead_time", y=param,
                        title=f"MAE Time‑series for {row['names']} ({copernicus_to_ecmwf[param]['unit']})",
                        line_dash = 'model',
                        color='model')
            
            ef = ef[:12]
            cf = cf[:12]
            
            score_val = 1 - ef/cf
            score = pd.DataFrame({
                "lead_time": range(0,12),
                "score (MAE)": score_val,
                "type":'Data' 
            })
            zero = pd.DataFrame({
                "lead_time": range(0,12),
                "score (MAE)": [0]*12,
                "type":'Zero' 
            })
            conc = pd.concat([score,zero])
            fig_score = px.line(conc, x="lead_time", y="score (MAE)",
                        title=f"{param} skill score for {station_id}",
                        line_dash="type",
                        color="type")
            
            ts_cop = pd.DataFrame({
                "lead_time": range(0,len(cf)),
                param: cf,
                'model': 'Copernicus'
            })
            fig_cop = px.line(ts_cop, x="lead_time", y=param,
                        title=f"MAE Time‑series for {row['names']} ({copernicus_to_ecmwf[param]['unit']})",
                        line_dash = 'model',
                        color='model')
            
            figures = [fig_score, fig_ecmwf, fig_cop]
            graphs = [dcc.Graph(figure=fig, style={"height": "400px"}) for fig in figures]
            
            return graphs

if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)



