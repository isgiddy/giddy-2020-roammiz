def read_ship_CTD(ctd_file):

    """Read the S.A. Agulhas II full CTD casts.

    Args:
      ctd_file (str): CTD data file

    Return:
        df: dataframe of the data
        date: date of the cast
        lon: longitude of the cast (deg.decimal)
        lat: latitude of the cast (deg.decimal)

    Dependencies:
        pandas
        numpy

    """
    import pandas as pd
    import numpy as np

    filePath = ctd_file
    df = pd.read_csv(filePath, delim_whitespace=True, skiprows=370, header=None)

    names = ['Temperature', 'Salinity', 'Conductivity [S/m]', 'Pressure[db]',
    't190C: Temperature, 2 [ITS-90, deg C]', 'c1S/m: Conductivity, 2 [S/m]',
    'Potential Temperature [ITS-90, deg C]',  'Potential Temperature, 2 [ITS-90, deg C]',
    'Beam Attenuation, WET Labs C-Star [1/m]', 'Beam Attenuation, WET Labs C-Star, 2 [1/m]',
    'sbeox0ML/L: Oxygen, SBE 43 [ml/l]', 'sbeox0PS: Oxygen, SBE 43 [% saturation]', 'sbeox0Mm/Kg: Oxygen, SBE 43 [umol/kg]',
    'Turbidity, WET Labs ECO BB [m^-1/sr]', 'Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]', 'Turbidity, WET Labs ECO [NTU]',
    'Altimeter [m]', 'PAR/Irradiance, Biospherical/Licor', 'SPAR/Surface Irradiance', 'Salinity, Practical [PSU]', 'Sound Velocity [Chen-Millero, m/s]',
    'Density [density, kg/m^3]', 'flag']

    df.columns = names

    x = pd.read_csv(filePath, sep='/t', skiprows=9, header=None)
    
    # read the latitude
    deg=str(x.iloc[0])[23:25]
    minsec=str(x.iloc[0])[26:31]
    
    deg=np.array(deg).astype(float)
    minsec=np.array(minsec).astype(float)
    
    lat = deg+(minsec/60)
    lat = round(lat, 2)
    
    # if cast is in S Hemisphere, make negative
    if str(x.iloc[0])[32]=='S':
        lat = -lat
        
    # read the longitude
    deg=str(x.iloc[1])[24:27]
    minsec=str(x.iloc[1])[28:33]

    deg=np.array(deg).astype(float)
    minsec=np.array(minsec).astype(float)

    lon = deg+(minsec/60)
    lon = round(lon, 2)
    # if cast is in W Hemisphere, make negative
    if str(x.iloc[0])[32]=='W':
        lon = -lon

    # read the date
    date=str(x.iloc[2])[25:45]
    date=pd.to_datetime(date)

    return df, date, lon, lat