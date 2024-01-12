# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 15:18:10 2022

This code imports the pixel data from HLS imagery for v2.0

@author: EmilyMyers
"""

# Imports
import argparse
import os
import numpy as np
import rasterio
from rasterio.windows import Window
import pandas as pd

# Read in command line arguments. Expecting a year value.
parser = argparse.ArgumentParser()
parser.add_argument("phen_name",type=str)
args = parser.parse_args()

# Functions
def return_year_doy(image_name):
    """
    Parameters
    ----------
    image_name : string containing the name of an HLS image

    Returns
    -------
    year : int of the year the image was collected
    doy : int of the day of year the image was collected

    """
    year = int(image_name[15:19])
    doy = int(image_name[19:22])
    return year,doy

def create_image_path(path1,path2,S30_or_L30,year):
    """
    Parameters
    ----------
    path1 : string containing the start of the image path
    path2 : string containing the third part of the image path
    S30_or_L30 : string specifying whether we are looking at L30 or S30 imagery
    year : string or int specifying the year of the imagery we are looking at

    Returns
    -------
    image_path : string containing the path to all L30 or S30 images in a particular year

    """
    assert(S30_or_L30=='S30' or S30_or_L30=='L30'), "Input to S30_or_L30 must be in the form of 'S30' or 'L30'"
    year = str(year)
    image_path = path1+S30_or_L30+'/'+year+path2
    return image_path

def return_phenocam_pixels(row,col,image_path,flatten=True):
    """
    Parameters
    ----------
    row : row of target phenocam
    col : col of target phenocam
    image_path : string containing the path to the image
    flatten : whether or not the returned pixels will be flattened to a single dimension

    Returns
    -------
    w : the twelve-pixel window around the phenocam of interest
    
    """
    with rasterio.open(image_path,driver='GTiff') as src:
        w = src.read(window = Window(col-1,row-2,3,4))
        w = np.squeeze(w)
        w_flat = np.ndarray.flatten(w)
        if flatten == True:
            w = w_flat
        if np.mean(w_flat) == -9999:
            w = None
    return w

def return_mean_std(w):
    """
    Parameters
    ----------
    w : 12-pixel window around a phenocam of interest, numpy array

    Returns
    -------
    mean_north : the mean value of the pixels around and north of the phenocam
    mean_center : the mean value of the pixels around the phenocam
    std_north : the standard deviation of the pixels around and north of the phenocam
    std_center : the standard deviation of the pixels around the phenocam

    """
    if w is None:
        mean_north = None; mean_center = None; std_north = None; std_center = None;
    else:
        w_flat = np.ndarray.flatten(w)
        if np.size(w_flat) != 12:
            raise Exception("Input array is not the expected size")
        mean_north = np.mean(w_flat[0:9])
        mean_center = np.mean(w_flat[3:12])
        std_north = np.std(w_flat[0:9])
        std_center = np.std(w_flat[3:12])
    return mean_north,mean_center,std_north,std_center

def return_phenocam_row_col(phenocam_name):
    """
    Parameters
    ----------
    phenocam_name : string containing phenocam name

    Returns
    -------
    lat : latitude of input phenocam
    lon : longitude of input phenocam

    """
    if phenocam_name == 'jershrubland':
        row,col = 2809, 922
    elif phenocam_name == 'jershrubland2':
        row,col = 2817, 943
    elif phenocam_name == 'jernovel':
        row,col = 2892, 917
    elif phenocam_name == 'jernovel2':
        row,col = 2881, 935
    elif phenocam_name == 'jergrassland':
        row,col = 3112, 930
    elif phenocam_name == 'jergrassland2':
        row,col = 3109, 953
    elif phenocam_name == 'jerbajada':
        row,col = 3137, 1554
    elif phenocam_name == 'jernort':
        row,col = 2987, 1074
    elif phenocam_name == 'ibp':
        row,col = 3088, 894
    elif phenocam_name == 'jernwern':
        row,col = 2957, 1229
    elif phenocam_name == 'NEON.D14.JORN.DP1.00033':
        row,col = 3086, 902
    elif phenocam_name == 'jersand':
        row,col = 3370, 1063
    else:
        raise Exception("Phenocam name not recognized")
    return row,col

"""
Paths and Variables
"""
image_path_1 = 'path_to_hls_imagery' # File location of L30 and S30 folders created by HLS bulk download
image_path_2 = '/13/S/C/S'
phenocam = args.phen_name
row,col = return_phenocam_row_col(phenocam)
years_L30 = [2014,2015,2016,2017,2018,2019,2020,2021,2022]
years_S30 = [2016,2017,2018,2019,2020,2021,2022]
output_dir = 'data/outputs_hls/'

"""
Initialize data frames
"""
empty_data_frame_window = pd.DataFrame(columns=('Year',
                                                'DOY',
                                                'Satellite',
                                                'Phenocam',
                                                'CoastalAerosol',
                                                'Blue',
                                                'Green',
                                                'Red',
                                                'RedEdge1',
                                                'RedEdge2',
                                                'RedEdge3',
                                                'NIRBroad',
                                                'NIRNarrow',
                                                'SWIR1',
                                                'SWIR2',
                                                'WaterVapor',
                                                'Cirrus',
                                                'TIR1',
                                                'TIR2',
                                                'Quality'))

empty_data_frame_meanstd = pd.DataFrame(columns=('Year',
                                                 'DOY',
                                                 'Satellite',
                                                 'Phenocam',
                                                 'CenterOrNorth',
                                                 'CoastalAerosol_mean',
                                                 'CoastalAerosol_std',
                                                 'Blue_mean',
                                                 'Blue_std',
                                                 'Green_mean',
                                                 'Green_std',
                                                 'Red_mean',
                                                 'Red_std',
                                                 'RedEdge1_mean',
                                                 'RedEdge1_std',
                                                 'RedEdge2_mean',
                                                 'RedEdge2_std',
                                                 'RedEdge3_mean',
                                                 'RedEdge3_std',
                                                 'NIRBroad_mean',
                                                 'NIRBroad_std',
                                                 'NIRNarrow_mean',
                                                 'NIRNarrow_std',
                                                 'SWIR1_mean',
                                                 'SWIR1_std',
                                                 'SWIR2_mean',
                                                 'SWIR2_std',
                                                 'WaterVapor_mean',
                                                 'WaterVapor_std',
                                                 'Cirrus_mean',
                                                 'Cirrus_std',
                                                 'Quality'))

"""
Extract band info and save it in pandas dataframe
"""
df_window = empty_data_frame_window
df_center = empty_data_frame_meanstd
df_north = empty_data_frame_meanstd

satellite = 'L30'
for y in years_L30:
    year = y
    image_path = create_image_path(image_path_1, image_path_2, satellite, year)
    image_list = os.listdir(image_path)
    
    #for i in range(3):
    for i in image_list:
        im_name = i
        year,doy = return_year_doy(im_name)
        temp_path = image_path+'/'+im_name+'/'+im_name
        
        coastalaerosol = return_phenocam_pixels(row, col, temp_path+'.B01.tif')
        blue = return_phenocam_pixels(row, col, temp_path+'.B02.tif')
        green = return_phenocam_pixels(row, col, temp_path+'.B03.tif')
        red = return_phenocam_pixels(row, col, temp_path+'.B04.tif')
        rededge1 = None
        rededge2 = None
        rededge3 = None
        nirbroad = None
        nirnarrow = return_phenocam_pixels(row, col, temp_path+'.B05.tif')
        swir1 = return_phenocam_pixels(row, col, temp_path+'.B06.tif')
        swir2 = return_phenocam_pixels(row, col, temp_path+'.B07.tif')
        watervapor = None
        cirrus = return_phenocam_pixels(row, col, temp_path+'.B09.tif')
        tir1 = return_phenocam_pixels(row, col, temp_path+'.B10.tif')
        tir2 = return_phenocam_pixels(row, col, temp_path+'.B11.tif')
        quality = return_phenocam_pixels(row, col, temp_path+'.Fmask.tif')
        temp_data = {'Year':year,
                     'DOY':doy,
                     'Satellite':satellite,
                     'Phenocam':phenocam,
                     'CoastalAerosol':coastalaerosol,
                     'Blue':blue,
                     'Green':green,
                     'Red':red,
                     'RedEdge1':rededge1,
                     'RedEdge2':rededge2,
                     'RedEdge3':rededge3,
                     'NIRBroad':nirbroad,
                     'NIRNarrow':nirnarrow,
                     'SWIR1':swir1,
                     'SWIR2':swir2,
                     'WaterVapor':watervapor,
                     'Cirrus':cirrus,
                     'TIR1':tir1,
                     'TIR2':tir2,
                     'Quality':quality}
        
        df_window = pd.concat([df_window,pd.DataFrame.from_dict([temp_data])],ignore_index=True)
        
        ca_mean_north,ca_mean_center,ca_std_north,ca_std_center = return_mean_std(coastalaerosol)
        b_mean_north,b_mean_center,b_std_north,b_std_center = return_mean_std(blue)
        g_mean_north,g_mean_center,g_std_north,g_std_center = return_mean_std(green)
        r_mean_north,r_mean_center,r_std_north,r_std_center = return_mean_std(red)
        re1_mean_north,re1_mean_center,re1_std_north,re1_std_center = return_mean_std(rededge1)
        re2_mean_north,re2_mean_center,re2_std_north,re2_std_center = return_mean_std(rededge2)
        re3_mean_north,re3_mean_center,re3_std_north,re3_std_center = return_mean_std(rededge3)
        n_b_mean_north,n_b_mean_center,n_b_std_north,n_b_std_center = return_mean_std(nirbroad)
        n_n_mean_north,n_n_mean_center,n_n_std_north,n_n_std_center = return_mean_std(nirnarrow)
        sw1_mean_north,sw1_mean_center,sw1_std_north,sw1_std_center = return_mean_std(swir1)
        sw2_mean_north,sw2_mean_center,sw2_std_north,sw2_std_center = return_mean_std(swir2)
        wv_mean_north,wv_mean_center,wv_std_north,wv_std_center = return_mean_std(watervapor)
        cir_mean_north,cir_mean_center,cir_std_north,cir_std_center = return_mean_std(cirrus)
        qa_center = np.ndarray.flatten(quality); qa_center = qa_center[7]
        
        temp_data = {'Year':year,
                     'DOY':doy,
                     'Satellite':satellite,
                     'Phenocam':phenocam,
                     'CenterOrNorth':'north',
                     'CoastalAerosol_mean':ca_mean_north,
                     'CoastalAerosol_std':ca_std_north,
                     'Blue_mean':b_mean_north,
                     'Blue_std':b_std_north,
                     'Green_mean':g_mean_north,
                     'Green_std':g_std_north,
                     'Red_mean':r_mean_north,
                     'Red_std':r_std_north,
                     'RedEdge1_mean':re1_mean_north,
                     'RedEdge1_std':re1_std_north,
                     'RedEdge2_mean':re2_mean_north,
                     'RedEdge2_std':re2_std_north,
                     'RedEdge3_mean':re3_mean_north,
                     'RedEdge3_std':re3_std_north,
                     'NIRBroad_mean':n_b_mean_north,
                     'NIRBroad_std':n_b_std_north,
                     'NIRNarrow_mean':n_n_mean_north,
                     'NIRNarrow_std':n_n_std_north,
                     'SWIR1_mean':sw1_mean_north,
                     'SWIR1_std':sw1_std_north,
                     'SWIR2_mean':sw2_mean_north,
                     'SWIR2_std':sw2_std_north,
                     'WaterVapor_mean':wv_mean_north,
                     'WaterVapor_std':wv_std_north,
                     'Quality':qa_center}
        
        df_north = pd.concat([df_north,pd.DataFrame.from_dict([temp_data])],ignore_index=True)
        
        temp_data = {'Year':year,
                     'DOY':doy,
                     'Satellite':satellite,
                     'Phenocam':phenocam,
                     'CenterOrNorth':'center',
                     'CoastalAerosol_mean':ca_mean_center,
                     'CoastalAerosol_std':ca_std_center,
                     'Blue_mean':b_mean_center,
                     'Blue_std':b_std_center,
                     'Green_mean':g_mean_center,
                     'Green_std':g_std_center,
                     'Red_mean':r_mean_center,
                     'Red_std':r_std_center,
                     'RedEdge1_mean':re1_mean_center,
                     'RedEdge1_std':re1_std_center,
                     'RedEdge2_mean':re2_mean_center,
                     'RedEdge2_std':re2_std_center,
                     'RedEdge3_mean':re3_mean_center,
                     'RedEdge3_std':re3_std_center,
                     'NIRBroad_mean':n_b_mean_center,
                     'NIRBroad_std':n_b_std_center,
                     'NIRNarrow_mean':n_n_mean_center,
                     'NIRNarrow_std':n_n_std_center,
                     'SWIR1_mean':sw1_mean_center,
                     'SWIR1_std':sw1_std_center,
                     'SWIR2_mean':sw2_mean_center,
                     'SWIR2_std':sw2_std_center,
                     'WaterVapor_mean':wv_mean_center,
                     'WaterVapor_std':wv_std_center,
                     'Quality':qa_center}
        
        df_center = pd.concat([df_center,pd.DataFrame.from_dict([temp_data])],ignore_index=True)

satellite = 'S30'
for y in years_S30:
    year = y
    image_path = create_image_path(image_path_1, image_path_2, satellite, year)
    image_list = os.listdir(image_path)
    
    #for i in range(3):
    for i in image_list:
        im_name = i
        year,doy = return_year_doy(im_name)
        temp_path = image_path+'/'+im_name+'/'+im_name        
        coastalaerosol = return_phenocam_pixels(row, col, temp_path+'.B01.tif')
        blue = return_phenocam_pixels(row, col, temp_path+'.B02.tif')
        green = return_phenocam_pixels(row, col, temp_path+'.B03.tif')
        red = return_phenocam_pixels(row, col, temp_path+'.B04.tif')
        rededge1 = return_phenocam_pixels(row, col, temp_path+'.B05.tif')
        rededge2 = return_phenocam_pixels(row, col, temp_path+'.B06.tif')
        rededge3 = return_phenocam_pixels(row, col, temp_path+'.B07.tif')
        nirbroad = return_phenocam_pixels(row, col, temp_path+'.B08.tif')
        nirnarrow = return_phenocam_pixels(row, col, temp_path+'.B8A.tif')
        swir1 = return_phenocam_pixels(row, col, temp_path+'.B11.tif')
        swir2 = return_phenocam_pixels(row, col, temp_path+'.B12.tif')
        watervapor = return_phenocam_pixels(row, col, temp_path+'.B09.tif')
        cirrus = return_phenocam_pixels(row, col, temp_path+'.B10.tif')
        tir1 = None
        tir2 = None
        quality = return_phenocam_pixels(row, col, temp_path+'.Fmask.tif')
        temp_data = {'Year':year,
                     'DOY':doy,
                     'Satellite':satellite,
                     'Phenocam':phenocam,
                     'CoastalAerosol':coastalaerosol,
                     'Blue':blue,
                     'Green':green,
                     'Red':red,
                     'RedEdge1':rededge1,
                     'RedEdge2':rededge2,
                     'RedEdge3':rededge3,
                     'NIRBroad':nirbroad,
                     'NIRNarrow':nirnarrow,
                     'SWIR1':swir1,
                     'SWIR2':swir2,
                     'WaterVapor':watervapor,
                     'Cirrus':cirrus,
                     'TIR1':tir1,
                     'TIR2':tir2,
                     'Quality':quality}
        
        df_window = pd.concat([df_window,pd.DataFrame.from_dict([temp_data])],ignore_index=True)
        
        ca_mean_north,ca_mean_center,ca_std_north,ca_std_center = return_mean_std(coastalaerosol)
        b_mean_north,b_mean_center,b_std_north,b_std_center = return_mean_std(blue)
        g_mean_north,g_mean_center,g_std_north,g_std_center = return_mean_std(green)
        r_mean_north,r_mean_center,r_std_north,r_std_center = return_mean_std(red)
        re1_mean_north,re1_mean_center,re1_std_north,re1_std_center = return_mean_std(rededge1)
        re2_mean_north,re2_mean_center,re2_std_north,re2_std_center = return_mean_std(rededge2)
        re3_mean_north,re3_mean_center,re3_std_north,re3_std_center = return_mean_std(rededge3)
        n_b_mean_north,n_b_mean_center,n_b_std_north,n_b_std_center = return_mean_std(nirbroad)
        n_n_mean_north,n_n_mean_center,n_n_std_north,n_n_std_center = return_mean_std(nirnarrow)
        sw1_mean_north,sw1_mean_center,sw1_std_north,sw1_std_center = return_mean_std(swir1)
        sw2_mean_north,sw2_mean_center,sw2_std_north,sw2_std_center = return_mean_std(swir2)
        wv_mean_north,wv_mean_center,wv_std_north,wv_std_center = return_mean_std(watervapor)
        cir_mean_north,cir_mean_center,cir_std_north,cir_std_center = return_mean_std(cirrus)
        qa_center = np.ndarray.flatten(quality); qa_center = qa_center[7]
        
        temp_data = {'Year':year,
                     'DOY':doy,
                     'Satellite':satellite,
                     'Phenocam':phenocam,
                     'CenterOrNorth':'north',
                     'CoastalAerosol_mean':ca_mean_north,
                     'CoastalAerosol_std':ca_std_north,
                     'Blue_mean':b_mean_north,
                     'Blue_std':b_std_north,
                     'Green_mean':g_mean_north,
                     'Green_std':g_std_north,
                     'Red_mean':r_mean_north,
                     'Red_std':r_std_north,
                     'RedEdge1_mean':re1_mean_north,
                     'RedEdge1_std':re1_std_north,
                     'RedEdge2_mean':re2_mean_north,
                     'RedEdge2_std':re2_std_north,
                     'RedEdge3_mean':re3_mean_north,
                     'RedEdge3_std':re3_std_north,
                     'NIRBroad_mean':n_b_mean_north,
                     'NIRBroad_std':n_b_std_north,
                     'NIRNarrow_mean':n_n_mean_north,
                     'NIRNarrow_std':n_n_std_north,
                     'SWIR1_mean':sw1_mean_north,
                     'SWIR1_std':sw1_std_north,
                     'SWIR2_mean':sw2_mean_north,
                     'SWIR2_std':sw2_std_north,
                     'WaterVapor_mean':wv_mean_north,
                     'WaterVapor_std':wv_std_north,
                     'Quality':qa_center}
        
        df_north = pd.concat([df_north,pd.DataFrame.from_dict([temp_data])],ignore_index=True)
        
        temp_data = {'Year':year,
                     'DOY':doy,
                     'Satellite':satellite,
                     'Phenocam':phenocam,
                     'CenterOrNorth':'center',
                     'CoastalAerosol_mean':ca_mean_center,
                     'CoastalAerosol_std':ca_std_center,
                     'Blue_mean':b_mean_center,
                     'Blue_std':b_std_center,
                     'Green_mean':g_mean_center,
                     'Green_std':g_std_center,
                     'Red_mean':r_mean_center,
                     'Red_std':r_std_center,
                     'RedEdge1_mean':re1_mean_center,
                     'RedEdge1_std':re1_std_center,
                     'RedEdge2_mean':re2_mean_center,
                     'RedEdge2_std':re2_std_center,
                     'RedEdge3_mean':re3_mean_center,
                     'RedEdge3_std':re3_std_center,
                     'NIRBroad_mean':n_b_mean_center,
                     'NIRBroad_std':n_b_std_center,
                     'NIRNarrow_mean':n_n_mean_center,
                     'NIRNarrow_std':n_n_std_center,
                     'SWIR1_mean':sw1_mean_center,
                     'SWIR1_std':sw1_std_center,
                     'SWIR2_mean':sw2_mean_center,
                     'SWIR2_std':sw2_std_center,
                     'WaterVapor_mean':wv_mean_center,
                     'WaterVapor_std':wv_std_center,
                     'Quality':qa_center}
        
        df_center = pd.concat([df_center,pd.DataFrame.from_dict([temp_data])],ignore_index=True)

df_north = df_north.astype({"Quality":np.uint8})
df_center = df_center.astype({"Quality":np.uint8})

"""
Save results
"""
window_filename = output_dir+phenocam+'_window.pkl'
center_filename = output_dir+phenocam+'_center.csv'
north_filename = output_dir+phenocam+'_north.csv'

df_window.to_pickle(window_filename)
df_center.to_csv(center_filename)
df_north.to_csv(north_filename)

print(phenocam+" complete!")