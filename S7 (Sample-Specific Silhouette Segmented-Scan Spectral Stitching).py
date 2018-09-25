#!/bin/env python3.6
# encoding: utf-8

"""
Created on September - 17 - 2017
S7 is a tool to group signals across scans, reduce noise, and normalize intensities across m/z segments.
@author: Travis Thompson

"""


import numpy
import pandas
import datetime
from MSFileReader.MSFileReader import ThermoRawfile
import multiprocessing
import Multiprocessing_Functions as MP
import wx
import wx.lib.agw.genericmessagedialog as GMD
import wx.lib.dialogs
import re
import os
import sys
import threading
import errno
from collections import OrderedDict
import argparse

__version__ = "1.0.0"

numpy.warnings.filterwarnings('ignore')



##TODO Possibly add option not to normalize by injection time.

EVT_THREAD_ABORTED = wx.NewId()
EVT_THREAD_ERROR = wx.NewId()
EVT_THREAD_COMPLETED = wx.NewId()

digital_precision_64_128 = 2**-17
digital_precision_128_256 = 2**-16
digital_precision_256_512 = 2**-15
digital_precision_512_1024 = 2**-14
digital_precision_1024_2048 = 2**-13




def overlap_alignment(index, mass):
    return index.get_loc(mass, method = "nearest")


class ThreadErrorEvent(wx.PyEvent):
    def __init__(self, message):
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_THREAD_ERROR)
        self.message = message
        
        

class ThreadAbortedEvent(wx.PyEvent):
    def __init__(self, message):
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_THREAD_ABORTED)
        self.message = message
        
        
        
class ThreadCompletedEvent(wx.PyEvent):
    def __init__(self):
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_THREAD_COMPLETED)
        


## Class taken from https://wiki.wxpython.org/LongRunningTasks and modified.
# Thread class that executes processing
class S8Thread(threading.Thread):
    """S8 Thread Class."""
    def __init__(self, RAW_and_User_info):
        """Init S8 Thread Class."""
        threading.Thread.__init__(self)
        self.RAW_and_User_info = RAW_and_User_info
        self.want_abort = False
        self.daemon = True
        
        ## Silence the setting with copy warning because it is spurious and makes it harder
        ## to see actual warnings.
        pandas.set_option('chained_assignment',None)
        numpy.warnings.filterwarnings('ignore')
        
        
        # This starts the thread running on creation, but you could
        # also make the GUI thread responsible for calling this
        self.start() 
    
    
    def run(self):
        """Run S8 Thread."""
        
        try:
        
            self.CheckUserAbort()
            
            self.num_of_cores = multiprocessing.cpu_count()
                    
            self.Build_Writers()
            
            self.CheckUserAbort()
            
            self.ReadInPeakData()
            
            self.CheckUserAbort()
            
            self.Calculate_IC_Correction()
            
            self.CheckUserAbort()
            
            front_segments_to_skip = self.RAW_and_User_info.front_segments_to_skip
            rear_segments_to_skip = self.RAW_and_User_info.rear_segments_to_skip
            scan_data_by_mass = self.scan_data_by_mass
            mass_ranges = self.RAW_and_User_info.mass_ranges
            self.previous_dataframe = 0
            self.previous_num_of_scans = 0
            self.noise_row_count = 0
            self.full_data_row_count = 0
            self.full_segments = []
            self.overlap_segments = []
            self.overlap_ratios = []
            
            for i in range(0+front_segments_to_skip, len(scan_data_by_mass) - rear_segments_to_skip):
                ## If there is no data in the segment then we can skip doing anything,
                ## but the overlap region needs to be managed.
                if len(scan_data_by_mass[i]) == 0:
                    message = "No data for segment" + str(mass_ranges[i])
                    wx.PostEvent(self.RAW_and_User_info, ThreadErrorEvent(message))
                    self.previous_dataframe = pandas.DataFrame()
                    self.previous_num_of_scans = 0
                    continue
                
                
                intensity_dataframe = pandas.DataFrame(scan_data_by_mass[i])
                self.current_mass_range_index = i
                self.intensity_dataframe = intensity_dataframe
                self.current_mass_range = mass_ranges[i]
                self.precision_multiple = int(self.RAW_and_User_info.segment_settings["Bin Widths"][self.current_mass_range])
                
                self.CheckUserAbort()
                
                
                self.GroupPeaksAcrossScans()
                
                self.CheckUserAbort()
                
                number_of_scans = self.number_of_scans
                ## Normalize intensities by injection time.
                ## Find max injection time. Injection time should be in the second row.
                max_injection_time = intensity_dataframe.max(axis=1)[1]
                ## Loop over all columns that are scans.
                for k in range(0,number_of_scans):
                    intensity_dataframe.iloc[5:,k] = max_injection_time/intensity_dataframe.iloc[1,k]*intensity_dataframe.iloc[5:,k]
                    
                
                
                self.CalculateBucketStats()
                
                self.CheckUserAbort()
                
                self.CorrectBasedOnInternalCalibrant()
                
                self.CheckUserAbort()
                
                self.RemoveNoise()
                
                self.CheckUserAbort()
                
                self.CalculateOverlapRegion()
                
                self.CheckUserAbort()
                
                self.AddFinalColumnsAndRearrange()
                
                self.CheckUserAbort()
                
                self.ManageDataframePrinting()
                
                self.CheckUserAbort()
                
                
                
            self.WriteToFile()
            wx.PostEvent(self.RAW_and_User_info, ThreadCompletedEvent())
            
        except Exception as e:
            wx.PostEvent(self.RAW_and_User_info, ThreadAbortedEvent(repr(e)))
            return
            
            

        
    
    
    def abort(self):
        """Abort S8 Thread."""
        # Method for use by main thread to signal an abort.
        self.want_abort = True

        
     
    def CheckUserAbort(self):
        if self.want_abort:
            wx.PostEvent(self.RAW_and_User_info, ThreadAbortedEvent("User Aborted"))
            sys.exit()
            
            
            
    def overlap_alignment(index, mass):
        return index.get_loc(mass, method = "nearest")
    
    
    
    
    def Build_Writers(self):
        output_directory = self.RAW_and_User_info.output_directory
        raw_file_name = self.RAW_and_User_info.raw_file_name
        
        raw_file_name = re.split(r"\.raw", raw_file_name)[0]
        
        
        if self.RAW_and_User_info.flags["Full Data"] or self.RAW_and_User_info.flags["Noise Data"]:
            raw_file_directory = output_directory + "\\" + raw_file_name
            if not os.path.exists(raw_file_directory):
                try:
                    os.makedirs(raw_file_directory)
                    
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
        
            
        if self.RAW_and_User_info.flags["Peak List"]:            
            peak_list_directory = output_directory + "\\Peak_Lists"
            if not os.path.exists(peak_list_directory):
                try:
                    os.makedirs(peak_list_directory)
                    
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
        
        
        if self.RAW_and_User_info.flags["Noise Data"]:
            noise_writers = []
            noise_writer_index = 0
            noise_writer_template = raw_file_directory + "\\Noise_Reduction_"
            noise_writers.append(pandas.ExcelWriter(noise_writer_template + raw_file_name + "_1.xlsx", engine="xlsxwriter"))
        
        if self.RAW_and_User_info.flags["Full Data"]:
            full_data_writers = []
            full_data_writer_index = 0
            full_data_writer_template =  raw_file_directory + "\\Full_Data_"
            full_data_writers.append(pandas.ExcelWriter(full_data_writer_template + raw_file_name + "_1.xlsx", engine="xlsxwriter"))
        
        if self.RAW_and_User_info.flags["Peak List"]:
            peak_list_writer = pandas.ExcelWriter(peak_list_directory + "\\Peak_List_" + raw_file_name + ".xlsx", engine="xlsxwriter")



        if self.RAW_and_User_info.flags["Noise Data"]:
            self.noise_writer_index = noise_writer_index
            self.noise_writers = noise_writers
            self.noise_writer_template = noise_writer_template
        
        if self.RAW_and_User_info.flags["Full Data"]:
            self.full_data_writer_index = full_data_writer_index
            self.full_data_writer_template = full_data_writer_template
            self.full_data_writers = full_data_writers
        
        if self.RAW_and_User_info.flags["Peak List"]:
            self.peak_list_writer = peak_list_writer
        
        








    def ReadInPeakData(self):
        
        num_of_mass_ranges = self.RAW_and_User_info.num_of_mass_ranges
        num_of_scans = self.RAW_and_User_info.num_of_scans
        rawfile = self.RAW_and_User_info.rawfile
        mass_ranges = self.RAW_and_User_info.mass_ranges
        internal_calibration_masses = self.RAW_and_User_info.segment_settings["Internal Calibration Masses"]
        internal_calibration_tolerances = self.RAW_and_User_info.segment_settings["Internal Calibration Tolerances"]
        max_inj_times = self.RAW_and_User_info.segment_settings["Max Injection Time"]
        
        scan_data_by_mass = [{} for _ in range(num_of_mass_ranges)]
        internal_calibrant_data_by_mass_range = [{} for _ in range(num_of_mass_ranges)]
        scans_with_no_internal_calibrant = []
        removed_scans = []
        
        for i in range(1,num_of_scans+1):
            
            trailer_extra = rawfile.GetTrailerExtraForScanNum(i)
            
            current_mass_range = rawfile.GetMassRangeFromScanNum(i, 0)
            
            mass_range_index = mass_ranges.index(current_mass_range)
                
            label_data, label_flag = rawfile.GetLabelData(i)
            
            masses = pandas.Series(label_data.mass)
            
            ## Drop scans with no internal calibrant.
            if internal_calibration_masses[current_mass_range] != 0.0:
                closest_index = pandas.Index(masses).get_loc(internal_calibration_masses[current_mass_range], method="nearest")
                closest_mass = masses.iloc[closest_index]
                correction_factor = internal_calibration_masses[current_mass_range] - closest_mass
                tolerance = internal_calibration_tolerances[current_mass_range]
        ## Commenting this out for now 5-16-2018. Had this in for testing, but really needs to be a strict cut off.
    #            while abs(correction_factor) > tolerance and tolerance < .01:
    #                tolerance += .0001
                if abs(correction_factor) >= tolerance:
                    scans_with_no_internal_calibrant.append([current_mass_range, i, abs(correction_factor)])
                    removed_scans.append(i)
                    continue
            
                
                
                mass_accuracy = correction_factor/internal_calibration_masses[current_mass_range] * 1000000
                
            else:
                mass_accuracy = 1000000
                correction_factor = 0.0
                tolerance = internal_calibration_tolerances[current_mass_range]
        
            std = numpy.std(numpy.array(label_data.intensity), ddof = 1)
                    
            
            ## Drop scans with max injection time.
            if trailer_extra["Ion Injection Time (ms)"] >= max_inj_times[current_mass_range]:
                removed_scans.append(i)
                continue
            
            scan_data_by_mass[mass_range_index][i] = pandas.Series(((std, trailer_extra["Ion Injection Time (ms)"], 
                                                                     correction_factor, mass_accuracy, tolerance) + label_data.intensity),
                                                                   index = ([0, 1, 2, 3, 4] + list(masses)))
        
            if internal_calibration_masses[current_mass_range] != 0.0:
                masses.index = masses
                internal_calibrant_data_by_mass_range[mass_range_index][i] = masses[(masses > internal_calibration_masses[current_mass_range] - .01) & (masses < internal_calibration_masses[current_mass_range] + .01)]



        self.scan_data_by_mass = scan_data_by_mass
        self.internal_calibrant_data_by_mass_range = internal_calibrant_data_by_mass_range








    def Calculate_IC_Correction(self):
        
        front_segments_to_skip = self.RAW_and_User_info.front_segments_to_skip
        rear_segments_to_skip = self.RAW_and_User_info.rear_segments_to_skip
        internal_calibrant_data_by_mass_range = self.internal_calibrant_data_by_mass_range
        precision_multiples = self.RAW_and_User_info.segment_settings["Bin Widths"]
        internal_calibration_masses = self.RAW_and_User_info.segment_settings["Internal Calibration Masses"]
        internal_calibration_tolerances = self.RAW_and_User_info.segment_settings["Internal Calibration Tolerances"]
        mass_ranges = self.RAW_and_User_info.mass_ranges
        scan_data_by_mass = self.scan_data_by_mass
        num_of_cores = self.num_of_cores
        
        
        needs_quit = False
        correction_factors = []
        for i in range(0+front_segments_to_skip, len(scan_data_by_mass) - rear_segments_to_skip):
            
            current_mass_range = mass_ranges[i]
            
            if internal_calibration_masses[current_mass_range] == 0.0:
                correction_factors.append(0.0)
                continue
            
            ## Only take the masses that are around the internal calibrant.
            intensity_dataframe = pandas.DataFrame(internal_calibrant_data_by_mass_range[i])
            number_of_scans = len(intensity_dataframe.columns)
            temp = pandas.DataFrame(intensity_dataframe.index)
            temp.columns = ["m/z"]
            temp["shifted"] = temp.iloc[:,0].shift(-1)
            temp["difference"] = temp.iloc[:,1] - temp.iloc[:,0]
            
            ## Append all bucket beginnings toegether. New buckets start when nearest neighbor difference is 11 times the
            ## digital precision or greater.
            bucket_beginnings = temp[temp["m/z"] < 128].index[temp["difference"][temp["m/z"] < 128] >= digital_precision_64_128*precision_multiples[current_mass_range]]
            bucket_beginnings = bucket_beginnings.append(temp[(temp["m/z"] >= 128) & (temp["m/z"] < 256)].index[temp["difference"][(temp["m/z"] >= 128) & (temp["m/z"] < 256)] >= digital_precision_128_256*precision_multiples[current_mass_range]])
            bucket_beginnings = bucket_beginnings.append(temp[(temp["m/z"] >= 256) & (temp["m/z"] < 512)].index[temp["difference"][(temp["m/z"] >= 256) & (temp["m/z"] < 512)] >= digital_precision_256_512*precision_multiples[current_mass_range]])
            bucket_beginnings = bucket_beginnings.append(temp[(temp["m/z"] >= 512) & (temp["m/z"] < 1024)].index[temp["difference"][(temp["m/z"] >= 512) & (temp["m/z"] < 1024)] >= digital_precision_512_1024*precision_multiples[current_mass_range]])
            bucket_beginnings = bucket_beginnings.append(temp[temp["m/z"] > 1024].index[temp["difference"][temp["m/z"] > 1024] >= digital_precision_1024_2048*precision_multiples[current_mass_range]])
    
            low_indexes = bucket_beginnings + 1
            high_indexes = low_indexes
            low_indexes = numpy.insert(low_indexes, 0, 0)     
            high_indexes = numpy.append(high_indexes, len(intensity_dataframe))
            
            
    
            bucket_high_ends = pandas.DataFrame(intensity_dataframe.index[high_indexes-1])
            bucket_high_ends.columns = ["high_range"]
            bucket_high_ends["high_index"] = high_indexes
            bucket_high_ends["low_index"] = low_indexes
            bucket_high_ends["low_range"] = intensity_dataframe.index[low_indexes]
            bucket_high_ends["m/z Ranges"] = bucket_high_ends["low_range"].map(repr) + " - " + bucket_high_ends["high_range"].map(repr)
            bucket_high_ends.index = bucket_high_ends["low_index"]

            temp_frame = MP.parallelize_dataframe(bucket_high_ends, MP.parallel_range_calcs, num_of_cores, intensity_dataframe)
        
        
            temp_frame.index = bucket_high_ends.index                            
            bucket_high_ends["Median m/z in m/z Range"] = temp_frame["Median m/z in m/z Range"]
            bucket_high_ends["Number of Signals in m/z Range"] = temp_frame["Number of Signals in m/z Range"]         
                   
            ## Put temporary dataframe into the main dataframe.
            ## This aligns the m/z ranges with the starting m/z.
            bucket_high_ends.index = bucket_high_ends["low_range"]
            intensity_dataframe["Median m/z in m/z Range"] = bucket_high_ends["Median m/z in m/z Range"]
            intensity_dataframe["Number of Signals in m/z Range"] = bucket_high_ends["Number of Signals in m/z Range"]
            
            ## Filter down to only the buckets with at least 40% repetition rate.
            intensity_dataframe = intensity_dataframe[intensity_dataframe["Number of Signals in m/z Range"] > number_of_scans*.4]
            if len(intensity_dataframe) == 0:
                message = "No signals in internal calibrant range have a repetition rate greater than 40% for segment", i, "try using a different internal calibrant m/z."
                wx.PostEvent(self.RAW_and_User_info, ThreadErrorEvent(message))
                needs_quit = True
            else:
              
                ## Determine the correction factor for the current segment and add it to the list.
                closest_index = pandas.Index(intensity_dataframe.loc[:, "Median m/z in m/z Range"].dropna(how="all")).get_loc(internal_calibration_masses[current_mass_range], method="nearest")
                closest_mass = intensity_dataframe.loc[:, "Median m/z in m/z Range"].dropna(how="all").iloc[closest_index]
                correction_factor = internal_calibration_masses[current_mass_range] - closest_mass
                tolerance = internal_calibration_tolerances[current_mass_range]
        ## Commenting this out for now 5-16-2018. Had this in for testing, but really needs to be a strict cut off.
#                while abs(correction_factor) > tolerance and tolerance < .01:
#                    tolerance += .0001
                    
                if abs(correction_factor) < tolerance:
                    correction_factors.append(correction_factor)
                else:
                    message = "Couldn't find internal calibrant for segment.", i, "Try increasing tolerance limit.", tolerance, correction_factor
                    wx.PostEvent(self.RAW_and_User_info, ThreadErrorEvent(message))
                    needs_quit = True


        if needs_quit:
            wx.PostEvent(self.RAW_and_User_info, ThreadAbortedEvent(None))
            sys.exit()
        
        self.correction_factors = correction_factors
        return True







    def GroupPeaksAcrossScans(self):
        
        intensity_dataframe = self.intensity_dataframe
        precision_multiple = self.precision_multiple
        
        
        number_of_scans = len(intensity_dataframe.columns)
        intensity_dataframe["m/z Ranges"] = numpy.nan
        intensity_dataframe["Number of Signals in m/z Range"] = numpy.nan
        
        ## Find clusters by looking for local maxima.
        ## This is an experimental algorithm for better bucketing, but it has a few problems 
        ## that need to be worked out. Leaving it for now.
        temp = pandas.DataFrame(intensity_dataframe.index)
        temp.columns = ["m/z"]
        temp["shifted"] = temp.iloc[:,0].shift(-1)
        temp["difference"] = temp.iloc[:,1] - temp.iloc[:,0]
        temp.iloc[-1,2] = 0.0
        temp.iloc[0:5, 2] = 0
        
        ## Append all bucket beginnings toegether. New buckets start when nearest neighbor difference is precision_multiple times the
        ## digital precision or greater.
        bucket_beginnings = temp[temp["m/z"] < 128].index[temp["difference"][temp["m/z"] < 128] >= digital_precision_64_128*precision_multiple]
        bucket_beginnings = bucket_beginnings.append(temp[(temp["m/z"] >= 128) & (temp["m/z"] < 256)].index[temp["difference"][(temp["m/z"] >= 128) & (temp["m/z"] < 256)] >= digital_precision_128_256*precision_multiple])
        bucket_beginnings = bucket_beginnings.append(temp[(temp["m/z"] >= 256) & (temp["m/z"] < 512)].index[temp["difference"][(temp["m/z"] >= 256) & (temp["m/z"] < 512)] >= digital_precision_256_512*precision_multiple])
        bucket_beginnings = bucket_beginnings.append(temp[(temp["m/z"] >= 512) & (temp["m/z"] < 1024)].index[temp["difference"][(temp["m/z"] >= 512) & (temp["m/z"] < 1024)] >= digital_precision_512_1024*precision_multiple])
        bucket_beginnings = bucket_beginnings.append(temp[temp["m/z"] > 1024].index[temp["difference"][temp["m/z"] > 1024] >= digital_precision_1024_2048*precision_multiple])

        low_indexes = bucket_beginnings + 1
        high_indexes = low_indexes
        low_indexes = numpy.insert(low_indexes, 0, 5)     
        high_indexes = numpy.append(high_indexes, len(intensity_dataframe))
        
        

        bucket_high_ends = pandas.DataFrame(intensity_dataframe.index[high_indexes-1])
        bucket_high_ends.columns = ["high_range"]
        bucket_high_ends["high_index"] = high_indexes
        bucket_high_ends["low_index"] = low_indexes
        bucket_high_ends["low_range"] = intensity_dataframe.index[low_indexes]
        bucket_high_ends["m/z Ranges"] = bucket_high_ends["low_range"].map(repr) + " - " + bucket_high_ends["high_range"].map(repr)
        bucket_high_ends.index = bucket_high_ends["low_index"]


        self.number_of_scans = number_of_scans
        self.bucket_high_ends = bucket_high_ends








    def CalculateBucketStats(self):
        
        bucket_high_ends = self.bucket_high_ends
        intensity_dataframe = self.intensity_dataframe
        num_of_cores = self.num_of_cores
        number_of_scans = self.number_of_scans
        
        ## Calculate stats for each group of peaks across scans.
        temp_frame = MP.parallelize_dataframe(bucket_high_ends, MP.parallel_range_calcs, num_of_cores, intensity_dataframe)
        temp_frame.index = bucket_high_ends.index                            
        bucket_high_ends["Number of Signals in m/z Range"] = temp_frame["Number of Signals in m/z Range"]
        bucket_high_ends["(Sum of Intensities)/(Number of Scans)"] = temp_frame["(Sum of Intensities)/(Number of Scans)"]
        bucket_high_ends["RSTD of Intensities/(Number of Scans)"] = temp_frame["RSTD of Intensities/(Number of Scans)"]
        bucket_high_ends["Median m/z in m/z Range"] = temp_frame["Median m/z in m/z Range"]
        bucket_high_ends["m/z RSTD in m/z Range"] = temp_frame["m/z RSTD in m/z Range"]

        ## Put temporary dataframe into the main dataframe for printing to excel.
        ## This aligns the m/z ranges with the starting m/z.
        bucket_high_ends.index = bucket_high_ends["low_range"]
        intensity_dataframe["m/z Ranges"] = bucket_high_ends["m/z Ranges"]
        intensity_dataframe["Number of Signals in m/z Range"] = bucket_high_ends["Number of Signals in m/z Range"]
        intensity_dataframe["Repetition Rate (%)"] = intensity_dataframe["Number of Signals in m/z Range"] / number_of_scans * 100
        intensity_dataframe["(Sum of Intensities)/(Number of Scans)"] = bucket_high_ends["(Sum of Intensities)/(Number of Scans)"]
        intensity_dataframe["RSTD of Intensities/(Number of Scans)"] = bucket_high_ends["RSTD of Intensities/(Number of Scans)"]
        intensity_dataframe["Median m/z in m/z Range"] = bucket_high_ends["Median m/z in m/z Range"]
        intensity_dataframe["m/z RSTD in m/z Range"] = bucket_high_ends["m/z RSTD in m/z Range"]







    def CorrectBasedOnInternalCalibrant(self):
        
        intensity_dataframe = self.intensity_dataframe
        i = self.current_mass_range_index
        correction_factors = self.correction_factors
        mass_ranges = self.RAW_and_User_info.mass_ranges
        rear_segments_to_skip = self.RAW_and_User_info.rear_segments_to_skip
        scan_data_by_mass = self.scan_data_by_mass
        internal_calibration_masses = self.RAW_and_User_info.segment_settings["Internal Calibration Masses"]
        
        ## Add new column for original median m/z.
        intensity_dataframe["Original Median m/z"] = intensity_dataframe.loc[:, "Median m/z in m/z Range"]
        
        
        ## Add correction factors to the average m/z's.
        intensity_dataframe.loc[:, "Median m/z in m/z Range"] = intensity_dataframe.loc[:, "Median m/z in m/z Range"] + correction_factors[i]
        
        current_mass_range = mass_ranges[i]
        
        ## Adjust the right side overlap correction factor.
        if i < len(scan_data_by_mass)-1-rear_segments_to_skip and \
        internal_calibration_masses[current_mass_range] != 0.0 and \
        internal_calibration_masses[mass_ranges[i+1]] != 0.0:
            right_correction_factor = (correction_factors[i] + correction_factors[i+1])/2
            true_correction_factor = right_correction_factor - correction_factors[i]
            overlap_max_mass = mass_ranges[i][1]
            overlap_min_mass = mass_ranges[i+1][0]
            intensity_dataframe.loc[
                    (intensity_dataframe.index >= overlap_min_mass - .1) 
                    & (intensity_dataframe.index <= overlap_max_mass + .1), "Median m/z in m/z Range"] = \
            intensity_dataframe.loc[
                    (intensity_dataframe.index >= overlap_min_mass - .1) 
                    & (intensity_dataframe.index <= overlap_max_mass + .1), "Median m/z in m/z Range"] + true_correction_factor
            
            
        ## Adjust the left side overlap correction factor.
        if i > 0 and \
        internal_calibration_masses[current_mass_range] != 0.0 and \
        internal_calibration_masses[mass_ranges[i-1]] != 0.0:
            left_correction_factor = (correction_factors[i-1] + correction_factors[i])/2
            true_correction_factor = left_correction_factor - correction_factors[i]
            overlap_max_mass = mass_ranges[i-1][1]
            overlap_min_mass = mass_ranges[i][0]
            intensity_dataframe.loc[
                    (intensity_dataframe.index >= overlap_min_mass - .1) 
                    & (intensity_dataframe.index <= overlap_max_mass + .1), "Median m/z in m/z Range"] = \
            intensity_dataframe.loc[
                    (intensity_dataframe.index >= overlap_min_mass - .1) 
                    & (intensity_dataframe.index <= overlap_max_mass + .1), "Median m/z in m/z Range"] + true_correction_factor

        

        ## Add new column for median correction factor.
        intensity_dataframe["Correction Factor"] = intensity_dataframe.loc[:, "Median m/z in m/z Range"] - intensity_dataframe.loc[:, "Original Median m/z"]










    def RemoveNoise(self):
        
        intensity_dataframe = self.intensity_dataframe
        number_of_scans = self.number_of_scans
        repetition_rate = self.RAW_and_User_info.repetition_rate
        
        
        ## Copy values down so we can easily eliminate the noise peaks.
        temp = intensity_dataframe.loc[:,"Repetition Rate (%)"].copy()
        intensity_dataframe["Repetition Rate (%)"].fillna(method="ffill", inplace=True)
        
        
        
        ## Set the ion injection time and standard deviation rows to the noise level
        ## so they won't get removed, and remove noise.
        intensity_dataframe.loc[[0,1,2,3,4],"Repetition Rate (%)"] = repetition_rate
        noise_dataframe = intensity_dataframe[intensity_dataframe["Repetition Rate (%)"] >= repetition_rate]
        ## Set the columns equal to temp so that there are blanks instead of copied down values.
        noise_dataframe.loc[:,"Repetition Rate (%)"] = temp
        intensity_dataframe.loc[:,"Repetition Rate (%)"] = temp
        
        
        ## Clear the first 3 rows so the entries aren't mixed in with the sorted values.
        intensity_dataframe.loc[[0,1,2,3,4],"m/z Ranges"] = numpy.nan
        intensity_dataframe.loc[[0,1,2,3,4],"Number of Signals in m/z Range"] = numpy.nan
        intensity_dataframe.loc[[0,1,2,3,4],"Repetition Rate (%)"] = numpy.nan
        intensity_dataframe.loc[[0,1,2,3,4],"(Sum of Intensities)/(Number of Scans)"] = numpy.nan
        intensity_dataframe.loc[[0,1,2,3,4],"RSTD of Intensities/(Number of Scans)"] = numpy.nan
        intensity_dataframe.loc[[0,1,2,3,4],"Median m/z in m/z Range"] = numpy.nan
        intensity_dataframe.loc[[0,1,2,3,4],"m/z RSTD in m/z Range"] = numpy.nan
        intensity_dataframe.loc[[0,1,2,3,4],"Original Median m/z"] = numpy.nan
        intensity_dataframe.loc[[0,1,2,3,4],"Correction Factor"] = numpy.nan
    
        
        
        ## Create new columns of sorted values so graph looks better.
        intensity_dataframe[["Sorted Number of Signals in m/z Range",
                             "Sorted Repetition Rate (%)",
                             "Sorted (Sum of Intensities)/(Number of Scans)",
                             "Sorted Median m/z"]] = \
                             intensity_dataframe.sort_values(
                                     ["Number of Signals in m/z Range",
                                      "(Sum of Intensities)/(Number of Scans)"]
                                     ).loc[:,
                             ["Number of Signals in m/z Range",
                              "Repetition Rate",
                              "(Sum of Intensities)/(Number of Scans)",
                              "Median m/z in m/z Range"]
                             ].set_index(intensity_dataframe.index)
       
                          
                             
        ## Clear the first 2 rows so the entries aren't mixed in with the sorted values.
        noise_dataframe.loc[[0,1,2,3,4],"m/z Ranges"] = numpy.nan
        noise_dataframe.loc[[0,1,2,3,4],"Number of Signals in m/z Range"] = numpy.nan
        noise_dataframe.loc[[0,1,2,3,4],"Repetition Rate (%)"] = numpy.nan
        noise_dataframe.loc[[0,1,2,3,4],"(Sum of Intensities)/(Number of Scans)"] = numpy.nan
        noise_dataframe.loc[[0,1,2,3,4],"RSTD of Intensities/(Number of Scans)"] = numpy.nan
        noise_dataframe.loc[[0,1,2,3,4],"Median m/z in m/z Range"] = numpy.nan
        noise_dataframe.loc[[0,1,2,3,4],"m/z RSTD in m/z Range"] = numpy.nan
        noise_dataframe.loc[[0,1,2,3,4],"Original Median m/z"] = numpy.nan
        noise_dataframe.loc[[0,1,2,3,4],"Correction Factor"] = numpy.nan
    
    
                         
        ## Create new columns of sorted values so graph looks better.
        noise_dataframe[["Sorted Number of Signals in m/z Range",
                         "Sorted Repetition Rate (%)",
                             "Sorted (Sum of Intensities)/(Number of Scans)",
                             "Sorted Median m/z"]] = \
                             noise_dataframe.sort_values(
                                     ["Number of Signals in m/z Range",
                                      "(Sum of Intensities)/(Number of Scans)"]
                                     ).loc[:,
                             ["Number of Signals in m/z Range",
                              "Repetition Rate (%)",
                              "(Sum of Intensities)/(Number of Scans)",
                              "Median m/z in m/z Range"]
                             ].set_index(noise_dataframe.index)
                                 
                             
        ## Count the number of intensity readings for each m/z.
        intensity_dataframe["Number of Scans m/z Appeared"] = intensity_dataframe.iloc[:,0:number_of_scans].count(axis=1)
        noise_dataframe["Number of Scans m/z Appeared"] = noise_dataframe.iloc[:,0:number_of_scans].count(axis=1)
        
        self.noise_dataframe = noise_dataframe
        
        
        
        
        
        
        
        
        
    def CalculateOverlapRegion(self):
        
        i = self.current_mass_range_index
        mass_ranges = self.RAW_and_User_info.mass_ranges
        previous_dataframe = self.previous_dataframe
        number_of_scans = self.number_of_scans
        noise_dataframe = self.noise_dataframe
        overlap_tolerance = self.RAW_and_User_info.segment_settings["Overlap Tolerances"][self.current_mass_range]
        previous_num_of_scans = self.previous_num_of_scans
        if self.RAW_and_User_info.flags["Noise Data"]:
            noise_writers = self.noise_writers
            noise_writer_index = self.noise_writer_index
        full_segments = self.full_segments
        overlap_ratios = self.overlap_ratios
        overlap_segments = self.overlap_segments
        
        
        
        ## Calculate overlapping region stuff.
        if i != 0 and len(previous_dataframe) > 0:
            ## Determine min and max mass for the overlap region.
            overlap_max_mass = mass_ranges[i-1][1]
            overlap_min_mass = mass_ranges[i][0]
            
            current_num_of_scans = number_of_scans
            
            
            ## Pull out only the rows in the 2 mass ranges that are in the overlap region.
            first = previous_dataframe[
                    (previous_dataframe["Median m/z in m/z Range"] >= overlap_min_mass - .1) 
                    & (previous_dataframe["Median m/z in m/z Range"] <= overlap_max_mass + .1)]
            
            first.dropna(inplace=True)
            
            #first = first[first["Number of Signals in m/z Range"] >= previous_num_of_scans]
            
            
            second = noise_dataframe.loc[:, 
                              ["Number of Signals in m/z Range", 
                               "(Sum of Intensities)/(Number of Scans)",
                               "Median m/z in m/z Range"]]
                    
            second = second[
                     (second["Median m/z in m/z Range"] >= overlap_min_mass - .1) 
                     & (second["Median m/z in m/z Range"] <= overlap_max_mass + .1)]
            
            second.dropna(inplace=True)
            
            #second = second[second["Number of Signals in m/z Range"] >= current_num_of_scans]
            
           
                    
            first.columns = ["First Range Signals", "First Range Intensities", "First Range Median m/z"]
            first = first.sort_values("First Range Median m/z")
            first.index = list(range(len(first)))
            
            
            second.columns = ["Second Range Signals", "Second Range Intensities", "Second Range Median m/z"]
            second = second.sort_values("Second Range Median m/z")
            second.index = list(range(len(second)))
            
            ## Find the indexes of the closest matching masses in first to the masses in second.
            alignment_indexes = second["Second Range Median m/z"].apply(lambda x: overlap_alignment(pandas.Index(first["First Range Median m/z"]), x))
            
            ## Put the matching masses in line for easy subtraction.
            matching_masses = first.loc[list(alignment_indexes)]["First Range Median m/z"]
            matching_masses.index = list(range(len(matching_masses)))
            
            ## Subtract the closest match found in the first from the mass in the second.
            mass_differences = second["Second Range Median m/z"] - matching_masses
            mass_differences = mass_differences.abs()
                    
            
            ## Drop the masses whose closest match has a difference greater than the input tolerance.
            second.drop(alignment_indexes[mass_differences > overlap_tolerance].index, inplace=True)
            alignment_indexes.drop(alignment_indexes[mass_differences > overlap_tolerance].index, inplace=True)
            
            ## Multiple m/z's can have a nearest match with a mass differences of less than .004.
            ## We have to have only one m/z match another, so we find the duplicates and keep the one with the lowest mass difference.
            duplicates = alignment_indexes[alignment_indexes[mass_differences <= overlap_tolerance].duplicated()]
            unique_duplicates = duplicates.iloc[:].unique()
            for g in range(len(unique_duplicates)):
                duplicates_indexes = list(duplicates[duplicates == unique_duplicates[g]].index)
                duplicates_indexes.append(duplicates_indexes[0]-1)
                smallest_diff = mass_differences.iloc[duplicates_indexes].idxmin()
                duplicates_indexes.remove(smallest_diff)
                second.drop(duplicates_indexes, inplace=True)
                alignment_indexes.drop(duplicates_indexes, inplace=True)
    
    
            ## Reindex the second dataframe before merging so the masses line up
            second.index = alignment_indexes
            
            ## Add columns for percent of max signals.
            first["First Percent of Max Signals"] = first["First Range Signals"]/previous_num_of_scans*100
            second["Second Percent of Max Signals"] = second["Second Range Signals"]/current_num_of_scans*100
            
            ## Merge the dataframes and drop rows that don't have a match.
            merged = pandas.concat([first, second], axis =1)
            merged.dropna(inplace = True)
            
            ## If merged has no length because none of the matches are strong enough then it will error.
            ## This check prevents that.
            if len(merged != 0):
                merged.index = list(range(len(merged)))
                ## Calculate ion ratios.
                merged.insert(merged.columns.get_loc("Second Range Signals"), "First Ion Ratio", merged["First Range Intensities"] / merged["First Range Intensities"].max())
                second_matching_index = merged["First Range Intensities"][merged["First Range Intensities"] == merged["First Range Intensities"].max()].index[0]
                merged["Second Ion Ratio"] = merged["Second Range Intensities"] / merged["Second Range Intensities"].iloc[second_matching_index]
                
                ## Commented this out to calculate the ion ratio in a different way, but I think we will go back to this.
                #merged["Second Ion Ratio"] = merged["Second Range Intensities"] / merged["Second Range Intensities"].max()
                merged["Median m/z"] = (merged["First Range Median m/z"] + merged["Second Range Median m/z"]) /2
                
                
                if self.RAW_and_User_info.flags["Noise Data"]:
                    ## Write the merged dataframe out to the Excel file for overlap region alignment.
                    sheet_name = "Overlap " + str(overlap_min_mass) + "-" + str(overlap_max_mass)
                    merged.to_excel(noise_writers[noise_writer_index], sheet_name = sheet_name)
                
                ## Make a copy of the full overlap region for the final peak list.
                ## Had to make this copy to make the 100% overlap matching work, but it is no longer necessary.
                ## Leaving it in for posterity right now. Will be removed in final version more than likely.
                original_merged = merged.copy()
                
    
                
            else:
                message = "No matching masses in overlap region " + str(overlap_min_mass) + "-" + str(overlap_max_mass)
                wx.PostEvent(self.RAW_and_User_info, ThreadErrorEvent(message))
                message = "Some m/z's will not have an intensity in the final peak list because the segment they are in could not be normalized."
                wx.PostEvent(self.RAW_and_User_info, ThreadErrorEvent(message))
    
            
    
            
        ## Save dataframe to use for overlapping regions.
        previous_dataframe = noise_dataframe.loc[:, 
                              ["Number of Signals in m/z Range", 
                               "(Sum of Intensities)/(Number of Scans)",
                               "Median m/z in m/z Range"]]
        
        previous_num_of_scans = number_of_scans
        
    
        if self.RAW_and_User_info.flags["Peak List"]:
            ## Add intensities and m/z's to the final peak list, dropping the overlap region 
            ## peaks to be handled later.
            full_segments.append(previous_dataframe.iloc[:,[1,2]].dropna())
            if i != 0 and len(original_merged) != 0:
                ## Drop the rows in the previous 2 segments that are in the overlap region.
                first_range_index = full_segments[i-1][full_segments[i-1].loc[:,"Median m/z in m/z Range"].isin(original_merged.loc[:,"First Range Median m/z"])].index
                full_segments[i-1].drop(first_range_index, inplace=True)
                second_range_index = full_segments[i][full_segments[i].loc[:,"Median m/z in m/z Range"].isin(original_merged.loc[:,"Second Range Median m/z"])].index
                full_segments[i].drop(second_range_index, inplace=True)
                ## Add the peaks in the overlap region to the overlap list.
                temp = original_merged.iloc[:, [1,6,10]]
                temp.index = original_merged.iloc[:, 10]
                temp.columns = ["First Range Intensities", "Second Range Intensities", "Median m/z in m/z Range"]
                overlap_segments.append(temp)
                overlap_ratios.append((merged["First Range Intensities"]/merged["Second Range Intensities"]).mean())
        
        
        self.previous_dataframe = previous_dataframe
        self.full_segments = full_segments
        self.overlap_segments = overlap_segments
        self.overlap_ratios = overlap_ratios
        self.previous_num_of_scans = previous_num_of_scans







    def AddFinalColumnsAndRearrange(self):
        
        intensity_dataframe = self.intensity_dataframe
        noise_dataframe = self.noise_dataframe
        
        ## Add the total number of buckets column.
        intensity_dataframe["Total Number of m/z Buckets"] = numpy.nan
        intensity_dataframe.loc[0,"Total Number of m/z Buckets"] = intensity_dataframe["Number of Signals in m/z Range"].count()
        noise_dataframe["Total Number of m/z Buckets"] = numpy.nan
        noise_dataframe.loc[0,"Total Number of m/z Buckets"] = noise_dataframe["Number of Signals in m/z Range"].count()
        
        ## Add a neighbor distance column to the intensity_dataframe.
        intensity_dataframe["Neighbor Distance"] = intensity_dataframe.index
        intensity_dataframe.loc[:,"Neighbor Distance"] = intensity_dataframe.loc[:,"Neighbor Distance"].shift(-1)
        intensity_dataframe.loc[:,"Neighbor Distance"] = intensity_dataframe.loc[:,"Neighbor Distance"] - intensity_dataframe.index
        ## Move the numbers down one position.
        intensity_dataframe.loc[:,"Neighbor Distance"] = intensity_dataframe.loc[:,"Neighbor Distance"].shift(1)
        ## Move the Neighbor Distance column to the front.
        row_count = intensity_dataframe["Neighbor Distance"]
        intensity_dataframe.drop("Neighbor Distance", axis = 1, inplace = True)
        intensity_dataframe.insert(0, "Neighbor Distance", row_count)
        
        ## Add a column to show integer value for neighbor distance.
        intensity_dataframe["Multiple of Precision"] = intensity_dataframe["Neighbor Distance"]
        intensity_dataframe.loc[intensity_dataframe.index < 128, "Multiple of Precision"] = intensity_dataframe.loc[intensity_dataframe.index < 128, "Multiple of Precision"] / digital_precision_64_128
        intensity_dataframe.loc[(intensity_dataframe.index >= 128) & (intensity_dataframe.index < 256), "Multiple of Precision"] = intensity_dataframe.loc[(intensity_dataframe.index >= 128) & (intensity_dataframe.index < 256), "Multiple of Precision"] / digital_precision_128_256
        intensity_dataframe.loc[(intensity_dataframe.index >= 256) & (intensity_dataframe.index < 512), "Multiple of Precision"] = intensity_dataframe.loc[(intensity_dataframe.index >= 256) & (intensity_dataframe.index < 512), "Multiple of Precision"] / digital_precision_256_512
        intensity_dataframe.loc[(intensity_dataframe.index >= 512) & (intensity_dataframe.index < 1024), "Multiple of Precision"] = intensity_dataframe.loc[(intensity_dataframe.index >= 512) & (intensity_dataframe.index < 1024), "Multiple of Precision"] / digital_precision_512_1024       
        intensity_dataframe.loc[intensity_dataframe.index > 1024, "Multiple of Precision"] = intensity_dataframe.loc[intensity_dataframe.index > 1024, "Multiple of Precision"] / digital_precision_1024_2048
        ## Move the Multiple of Precision column to the front.
        row_count = intensity_dataframe["Multiple of Precision"]
        intensity_dataframe.drop("Multiple of Precision", axis = 1, inplace = True)
        intensity_dataframe.insert(1, "Multiple of Precision", row_count)

        
#        ## Add the Half Width column to the intensity_dataframe.
#        intensity_dataframe["Half Width"] = half_widths.loc[:, "Half Width"]
#        ## Move the Half Width column to the front.
#        row_count = intensity_dataframe["Half Width"]
#        intensity_dataframe.drop("Half Width", axis = 1, inplace = True)
#        intensity_dataframe.insert(1, "Half Width", row_count)

        
        
        ## Rename the starting indexes so they have the right names.
        intensity_dataframe.rename(index = {0: "Standard Deviation", 
                                            1: "Ion Injection Time (ms)", 
                                            2: "Mass Correction Difference",
                                            3: "Mass Accuracy",
                                            4: "Internal Calibration Tolerance"}, inplace = True)
        ## Move Number of Scans m/z Appeared column to the left.
        row_count = intensity_dataframe["Number of Scans m/z Appeared"]
        intensity_dataframe.drop("Number of Scans m/z Appeared", axis = 1, inplace = True)
        intensity_dataframe.insert(intensity_dataframe.columns.get_loc("m/z Ranges"), "Number of Scans m/z Appeared", row_count)
        ## Clear the Standard Deviation and Ion Injection Time rows for the columns that don't make sense.
        intensity_dataframe.loc[["Standard Deviation", "Ion Injection Time (ms)", "Mass Correction Difference", "Mass Accuracy", "Internal Calibration Tolerance"],["Number of Scans m/z Appeared", "Neighbor Distance", "Multiple of Precision"]] = numpy.nan
        intensity_dataframe.iloc[5,0:2] = numpy.nan





        ## Add a neighbor distance column to the noise_dataframe.
        noise_dataframe["Neighbor Distance"] = noise_dataframe.index
        noise_dataframe.loc[:,"Neighbor Distance"] = noise_dataframe.loc[:,"Neighbor Distance"].shift(-1)
        noise_dataframe.loc[:,"Neighbor Distance"] = noise_dataframe.loc[:,"Neighbor Distance"] - noise_dataframe.index
        noise_dataframe.loc[:,"Neighbor Distance"] = noise_dataframe.loc[:,"Neighbor Distance"].shift(1)
        ## Move the Neighbor Distance column to the front.
        row_count = noise_dataframe["Neighbor Distance"]
        noise_dataframe.drop("Neighbor Distance", axis = 1, inplace = True)
        noise_dataframe.insert(0, "Neighbor Distance", row_count)
        
        ## Add a column to show integer value for neighbor distance.
        noise_dataframe["Multiple of Precision"] = noise_dataframe["Neighbor Distance"]
        noise_dataframe.loc[noise_dataframe.index < 256, "Multiple of Precision"] = noise_dataframe.loc[noise_dataframe.index < 256, "Multiple of Precision"] / 0.0000152587890625
        noise_dataframe.loc[(noise_dataframe.index >= 256) & (noise_dataframe.index < 512), "Multiple of Precision"] = noise_dataframe.loc[(noise_dataframe.index >= 256) & (noise_dataframe.index < 512), "Multiple of Precision"] / 0.000030517578125
        noise_dataframe.loc[(noise_dataframe.index >= 512) & (noise_dataframe.index < 1024), "Multiple of Precision"] = noise_dataframe.loc[(noise_dataframe.index >= 512) & (noise_dataframe.index < 1024), "Multiple of Precision"] / 0.00006103515625       
        noise_dataframe.loc[noise_dataframe.index > 1024, "Multiple of Precision"] = noise_dataframe.loc[noise_dataframe.index > 1024, "Multiple of Precision"] / 0.0001220703125
        ## Move the Multiple of Precision column to the front.
        row_count = noise_dataframe["Multiple of Precision"]
        noise_dataframe.drop("Multiple of Precision", axis = 1, inplace = True)
        noise_dataframe.insert(1, "Multiple of Precision", row_count)


#        ## Add the Half Width column to the noise_dataframe.
#        noise_dataframe["Half Width"] = half_widths.loc[:, "Half Width"]
#        ## Move the Half Width column to the front.
#        row_count = noise_dataframe["Half Width"]
#        noise_dataframe.drop("Half Width", axis = 1, inplace = True)
#        noise_dataframe.insert(1, "Half Width", row_count)



        noise_dataframe.rename(index = {0: "Standard Deviation", 
                                        1: "Ion Injection Time (ms)", 
                                        2:"Mass Correction Difference",
                                        3: "Mass Accuracy",
                                        4: "Internal Calibration Tolerance"}, inplace = True)
        ## Move Number of Scans m/z Appeared column to the left.
        row_count = noise_dataframe["Number of Scans m/z Appeared"]
        noise_dataframe.drop("Number of Scans m/z Appeared", axis = 1, inplace = True)
        noise_dataframe.insert(noise_dataframe.columns.get_loc("m/z Ranges"), "Number of Scans m/z Appeared", row_count)
        ## Clear the Number of Scans for rows that dont make sense.
        noise_dataframe.loc[["Standard Deviation", "Ion Injection Time (ms)", "Mass Correction Difference", "Mass Accuracy", "Internal Calibration Tolerance"],["Number of Scans m/z Appeared", "Neighbor Distance", "Multiple of Precision"]] = numpy.nan
        noise_dataframe.iloc[5,0:2] = numpy.nan
        
        
        
        
        
        
        
        
    def ManageDataframePrinting(self):
        
        intensity_dataframe = self.intensity_dataframe
        
        if self.RAW_and_User_info.flags["Full Data"]:
            full_data_row_count = self.full_data_row_count
            full_data_writers = self.full_data_writers
            full_data_writer_index = self.full_data_writer_index
            full_data_writer_template = self.full_data_writer_template
        
        if self.RAW_and_User_info.flags["Noise Data"]:
            noise_dataframe = self.noise_dataframe
            noise_row_count = self.noise_row_count
            noise_writers = self.noise_writers
            noise_writer_index = self.noise_writer_index
            noise_writer_template = self.noise_writer_template
        
        
        mass_ranges = self.RAW_and_User_info.mass_ranges
        i = self.current_mass_range_index
        raw_file_name = self.RAW_and_User_info.raw_file_name
        
        
        if self.RAW_and_User_info.flags["Noise Data"]:
            ## If there are greater than around one million rows in the data frame
            ## it has to be split up to write to an Excel file.
            if len(noise_dataframe) > 1000000:
                j = 2
                temp = []
                while len(noise_dataframe)/j > 1000000:
                    j += 1
                
                for g, frame in noise_dataframe.groupby(numpy.arange(len(noise_dataframe)) / (len(noise_dataframe) / j + 1)):
                    temp.append(frame)
                
                if noise_row_count + len(temp[0]) > 1000000:
                    for j in range(len(temp)-1):
                        noise_writer_index += 1
                        noise_writers.append(pandas.ExcelWriter(noise_writer_template + raw_file_name + "_" + str(noise_writer_index+1) + ".xlsx", engine="xlsxwriter"))
                        sheet_name = str(mass_ranges[i])
                        temp[j].to_excel(noise_writers[noise_writer_index], sheet_name = sheet_name)
                else:
                    sheet_name = str(mass_ranges[i])
                    temp[0].to_excel(noise_writers[noise_writer_index], sheet_name = sheet_name)
                    for j in range(1,len(temp)-1):
                        noise_writer_index += 1
                        noise_writers.append(pandas.ExcelWriter(noise_writer_template + raw_file_name + "_" + str(noise_writer_index+1) + ".xlsx", engine="xlsxwriter"))
                        sheet_name = str(mass_ranges[i])
                        temp[j].to_excel(noise_writers[noise_writer_index], sheet_name = sheet_name)
                        
                noise_row_count = 0
                noise_writer_index += 1
                noise_writers.append(pandas.ExcelWriter(noise_writer_template + raw_file_name + "_" + str(noise_writer_index+1) + ".xlsx", engine="xlsxwriter"))
                noise_dataframe = temp[-1]
            
        
        if self.RAW_and_User_info.flags["Full Data"]:
        
            if len(intensity_dataframe) > 1000000:
                j = 2
                temp = []
                while len(intensity_dataframe)/j > 1000000:
                    j += 1
                
                for g, frame in intensity_dataframe.groupby(numpy.arange(len(intensity_dataframe)) / (len(intensity_dataframe) / j + 1)):
                    temp.append(frame)
    
                if full_data_row_count + len(temp[0]) > 1000000:
                    for j in range(len(temp)-1):
                        full_data_writer_index += 1
                        full_data_writers.append(pandas.ExcelWriter(full_data_writer_template + raw_file_name + "_" + str(full_data_writer_index+1) + ".xlsx", engine="xlsxwriter"))
                        sheet_name = str(mass_ranges[i])
                        temp[j].to_excel(full_data_writers[full_data_writer_index], sheet_name = sheet_name)
                else:
                    sheet_name = str(mass_ranges[i])
                    temp[0].to_excel(full_data_writers[full_data_writer_index], sheet_name = sheet_name)
                    for j in range(1,len(temp)-1):
                        full_data_writer_index += 1
                        full_data_writers.append(pandas.ExcelWriter(full_data_writer_template + raw_file_name + "_" + str(full_data_writer_index+1) + ".xlsx", engine="xlsxwriter"))
                        sheet_name = str(mass_ranges[i])
                        temp[j].to_excel(full_data_writers[full_data_writer_index], sheet_name = sheet_name)
                        
                full_data_row_count = 0
                full_data_writer_index += 1
                full_data_writers.append(pandas.ExcelWriter(full_data_writer_template + raw_file_name + "_" + str(full_data_writer_index+1) + ".xlsx", engine="xlsxwriter"))
                intensity_dataframe = temp[-1]


        
        
        if self.RAW_and_User_info.flags["Noise Data"]:
            noise_row_count += len(noise_dataframe)
            if noise_row_count > 1000000:
                noise_row_count = 0
                noise_writer_index += 1
                noise_writers.append(pandas.ExcelWriter(noise_writer_template + raw_file_name + "_" + str(noise_writer_index+1) + ".xlsx", engine="xlsxwriter"))
        
        
        if self.RAW_and_User_info.flags["Full Data"]:
            full_data_row_count += len(intensity_dataframe)
            if full_data_row_count > 1000000:
                full_data_row_count = 0
                full_data_writer_index += 1
                full_data_writers.append(pandas.ExcelWriter(full_data_writer_template + raw_file_name + "_" + str(full_data_writer_index+1) + ".xlsx", engine="xlsxwriter"))
        
        
        
        ## Write the dataframe to the excel file.
        sheet_name = str(mass_ranges[i])
        
        if self.RAW_and_User_info.flags["Full Data"]:
            intensity_dataframe.to_excel(full_data_writers[full_data_writer_index], sheet_name = sheet_name)
            
        if self.RAW_and_User_info.flags["Noise Data"]:    
            noise_dataframe.to_excel(noise_writers[noise_writer_index], sheet_name = sheet_name)
        
        
        
        
        
        
        
        
        
    def WriteToFile(self):
        
        num_of_mass_ranges = self.RAW_and_User_info.num_of_mass_ranges
        
        
        if self.RAW_and_User_info.flags["Peak List"]:
            overlap_ratios = self.overlap_ratios
            full_segments = self.full_segments
            overlap_segments = self.overlap_segments
            peak_list_writer = self.peak_list_writer
        
        if self.RAW_and_User_info.flags["Noise Data"]:
            noise_writers = self.noise_writers
        
        if self.RAW_and_User_info.flags["Full Data"]:
            full_data_writers = self.full_data_writers
        
        
        
        if self.RAW_and_User_info.flags["Peak List"]:
            ## Only do this part if there are multiple segments.
            if num_of_mass_ranges > 1:    
                max_overlap_ratio = max(overlap_ratios)
                max_OR_index = overlap_ratios.index(max_overlap_ratio)
                multiplier = 1
                for i in range(max_OR_index, len(overlap_ratios)):
                    multiplier = multiplier*overlap_ratios[i]
                    full_segments[i+1].iloc[:,0] = full_segments[i+1].iloc[:,0]*multiplier
                    overlap_segments[i].iloc[:,1] = overlap_segments[i].iloc[:,1]*multiplier
                  
                divider = overlap_ratios[max_OR_index-1]**2
                for i in range(max_OR_index, 0, -1):
                    divider = divider/overlap_ratios[i-1]
                    full_segments[i-1].iloc[:,0] = full_segments[i-1].iloc[:,0]/divider
                    overlap_segments[i].iloc[:,0] = overlap_segments[i].iloc[:,0]/divider
            
                
            
            ## Combine the segments together
            stitched_df = pandas.concat(full_segments)
            
            ## Only do this part if there are multiple segments
            if num_of_mass_ranges > 1:
                ## Combine the overlap regions together.
                temp = pandas.concat(overlap_segments)
                ## Average the intensities in the overlap region.
                temp["(Sum of Intensities)/(Number of Scans)"] = (temp.iloc[:, 0] + temp.iloc[:, 1] )/ 2
                ## drop the unaveraged columns.
                temp.drop(["First Range Intensities", "Second Range Intensities"], axis=1, inplace = True)
                ## Combine it all together.
                stitched_df = pandas.concat([stitched_df, temp])
            
            stitched_df = stitched_df.sort_values(["Median m/z in m/z Range"])
            stitched_df.columns = ["Intensity", "m/z"]
            stitched_df.index = range(0, len(stitched_df))
            stitched_df.loc[:, "Relative Intensity"] = stitched_df.loc[:, "Intensity"]/stitched_df.loc[:, "Intensity"].max()*100
            row_count = stitched_df["Relative Intensity"]
            stitched_df.drop("Relative Intensity", axis = 1, inplace = True)
            stitched_df.insert(1, "Relative Intensity", row_count)
        
            stitched_df.to_excel(peak_list_writer, sheet_name = "Peak List", index = False)
        
        
        
        
        if self.RAW_and_User_info.flags["Noise Data"]:
            for i in range(0, len(noise_writers)):
                self.CheckUserAbort()
                noise_writers[i].save()
        
        if self.RAW_and_User_info.flags["Full Data"]:
            for i in range(0, len(full_data_writers)):
                self.CheckUserAbort()
                full_data_writers[i].save()
        
        if self.RAW_and_User_info.flags["Peak List"]:
            peak_list_writer.save()
        
        self.RAW_and_User_info.rawfile.Close()
        self.RAW_and_User_info.rawfile = None
        







## Frame to pass errors from the thread.
class Error_Frame(wx.Frame):
    def __init__(self, parent, *args, **kwargs):
        super(Error_Frame, self).__init__(parent, title="Algorithm Errors", size=(1100, 500), *args, **kwargs)
        
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        self.error_txt_ctrl = wx.TextCtrl(panel, style=wx.TE_READONLY|wx.TE_MULTILINE)
        vbox.Add(self.error_txt_ctrl, proportion=1, flag=wx.EXPAND)
        
        panel.SetSizer(vbox)
        panel.Fit()
        

        



class S8_GUI(wx.Frame):
    def __init__(self, parent, title, GUI_show = True, parser_args = None, *args, **kwargs):
        min_size=(1100,550)
        super(S8_GUI, self).__init__(parent, title=title, size=min_size, *args, **kwargs)
        self.SetMinSize(min_size)
        
        self.parser_args = parser_args
        
        if GUI_show:
            self.GUI_mode = True
        else:
            self.GUI_mode = False
            
        self.Connect(-1, -1, EVT_THREAD_ABORTED, self.thread_aborted)
        self.Connect(-1, -1, EVT_THREAD_ERROR, self.thread_error)
        self.Connect(-1, -1, EVT_THREAD_COMPLETED, self.thread_completed)
        
        self.S8_thread = None
        self.output_directory = None
        self.rawfile = None
        self.paths = None
        
        flags = ["Full Data", "Noise Data", "Peak List"]
        self.flags = OrderedDict()
        for flag in flags:
            self.flags[flag] = True
        
        ## Required columns for the input file if used.
        self.required_columns = ["Internal Calibration Masses", 
                                 "Internal Calibration Tolerances", 
                                 "Overlap Tolerances", 
                                 "Bin Widths",
                                 "Max Injection Time"]
        
        self.segment_settings = {}
        for setting in self.required_columns:
            self.segment_settings[setting] = {}
            
        
        self.repetition_rate = 10
        self.rear_segments_to_skip = 0
        self.front_segments_to_skip = 0
            
        if self.GUI_mode:
            self.InitUI()
            self.Centre()
            self.Show()
        else:
            
            ## rep_rate, and skips aren't actually set here so this is superfluous
            self.repetition_rate = self.parser_args.repetition_rate
            self.rear_segments_to_skip = self.parser_args.rear_skip
            self.front_segments_to_skip = self.parser_args.front_skip
            self.flags["Full Data"] = self.parser_args.full_data
            self.flags["Noise Data"] = self.parser_args.noise_data
            self.flags["Peak List"] = self.parser_args.peak_list
            self.output_directory = self.parser_args.output_directory
            self.paths = [self.parser_args.raw_file]
            self.current_raw_file_index = 0
            
            
            if not self.RAW_Open_and_Count(self.parser_args.raw_file):
                exit(1)
            if not self.Determine_Max_Injection_Time():
                exit(1)
            if not self.Load_From_File(None):
                exit(1)
            if not self.CheckOutputDirectory():
                exit(1)
            if not self.ValidateCurrentRawFile():
                exit(1)
            if not self.CheckReportOutputs():
                exit(1)
                
            self.S8_thread = S8Thread(self)
            self.Algorithm_error_printer("Starting run on file " + self.raw_file_name + ":")
            self.Algorithm_error_printer('Timestamp: {:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now()))


    def InitUI(self):
        
        ##############
        ## Error Frame
        ##############
        
        self.error_frame = Error_Frame(self)
        self.error_frame.Centre()
        self.error_frame.Show()
        
        
                
        ##############
        ## Menu Bar
        ##############
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        
        ## Add items to file menu.
        open_item = fileMenu.Append(wx.ID_OPEN, "&Open", "Open RAW File(s)")
        self.Bind(wx.EVT_MENU, self.OnOpen, open_item)
        
        fileMenu.AppendSeparator()
        
        quit_item = fileMenu.Append(wx.ID_EXIT, '&Quit', 'Quit Application')
        self.Bind(wx.EVT_MENU, self.OnQuit, quit_item)
        
        
        ## Add file menu to the menu bar.
        menubar.Append(fileMenu, "&File")
        ## Put menu bar in frame.
        self.SetMenuBar(menubar)
        
        ## Add Help menu.
        helpMenu = wx.Menu()
        
        about_item = helpMenu.Append(100, "&About", "Information about the program")
        self.Bind(wx.EVT_MENU, self.OnAbout, about_item)
        
        settings_item = helpMenu.Append(101, "Settings Description", "Description of program settings")
        self.Bind(wx.EVT_MENU, self.OnSettingsDescription, settings_item)
        
        menubar.Append(helpMenu, "&Help")



        panel = wx.Panel(self)
        self.statusbar = self.CreateStatusBar()
        
        vbox = wx.BoxSizer(wx.VERTICAL)

        
        ##############
        ## Current Raw File Display
        ##############
        current_raw_file_header = wx.StaticText(panel, label = "Current RAW File(s):")
        self.current_raw_file_txt_ctrl = wx.TextCtrl(panel, style=wx.TE_READONLY|wx.TE_MULTILINE)
        
        vbox.Add(current_raw_file_header, flag=wx.ALIGN_LEFT | wx.LEFT | wx.TOP, border=10)
        vbox.Add(self.current_raw_file_txt_ctrl, flag=wx.EXPAND | wx.ALIGN_LEFT | wx.LEFT, border=10)
        
        
        
        
        ##############
        ## Algorithm Settings Panel
        ##############
        
        self.segment_settings_panel = wx.Panel(panel)
        vbox.Add(self.segment_settings_panel, proportion=1, flag=wx.EXPAND)
        
                
        
        
        
        ###############
        ## Report Flags
        ###############
        
        check_boxes_label = wx.StaticText(panel, label="Report Outputs")
        vbox.Add(check_boxes_label, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.TOP, border=10)
        
        hbox_check_boxes = wx.BoxSizer(wx.HORIZONTAL)
        
        self.flag_checkboxes = {}
        for flag_name, flag_value in self.flags.items():
            self.flag_checkboxes[flag_name] = wx.CheckBox(panel, label = flag_name)
            self.flag_checkboxes[flag_name].SetValue(flag_value)
            self.flag_checkboxes[flag_name].Bind(wx.EVT_CHECKBOX, self.Set_Output_Flags)
            hbox_check_boxes.Add(self.flag_checkboxes[flag_name])
            
        
        vbox.Add(hbox_check_boxes, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, border=10)
        
        
        
        ###############
        ## Repetition Rate
        ###############
        
        repetition_rate_label = wx.StaticText(panel, label="Repetition Rate Filter %")
        self.repetition_rate_txt_control = wx.TextCtrl(panel)
        self.repetition_rate_txt_control.SetValue("10")
        self.repetition_rate_txt_control.Bind(wx.EVT_TEXT, self.OnRepetitionRate)
        vbox.Add(repetition_rate_label, flag = wx.TOP|wx.LEFT|wx.ALIGN_LEFT, border=10)
        vbox.Add(self.repetition_rate_txt_control, flag=wx.ALL|wx.ALIGN_LEFT, border=10)
        
        
        
        
        
        ###############
        ## Front and Rear Segments to Skip
        ###############
        
        front_skip_label = wx.StaticText(panel, label="Number of Front Segments to Skip")
        self.front_skip_txt_control = wx.TextCtrl(panel)
        self.front_skip_txt_control.SetValue("0")
        self.front_skip_txt_control.Bind(wx.EVT_TEXT, self.OnFrontSkip)
        vbox.Add(front_skip_label, flag = wx.TOP|wx.LEFT|wx.ALIGN_LEFT, border=10)
        vbox.Add(self.front_skip_txt_control, flag=wx.ALL|wx.ALIGN_LEFT, border=10)
        
        
        rear_skip_label = wx.StaticText(panel, label="Number of Rear Segments to Skip")
        self.rear_skip_txt_control = wx.TextCtrl(panel)
        self.rear_skip_txt_control.SetValue("0")
        self.rear_skip_txt_control.Bind(wx.EVT_TEXT, self.OnRearSkip)
        vbox.Add(rear_skip_label, flag = wx.TOP|wx.LEFT|wx.ALIGN_LEFT, border=10)
        vbox.Add(self.rear_skip_txt_control, flag=wx.ALL|wx.ALIGN_LEFT, border=10)
        
        
        
        
        
        ###############
        ## Choose Output Directory
        ###############
        
        ## Create widgets to select the directory to save final output in.
        dir_label = wx.StaticText(panel, label = "Save Output In:")
        self.dir_text = wx.TextCtrl(panel, style=wx.TE_READONLY)
        self.dir_button = wx.Button(panel, label="Select Folder")
        self.dir_button.Bind(wx.EVT_BUTTON, self.Dir_Select)
        
        ## Add widgets to panel.
        dir_hbox = wx.BoxSizer(wx.HORIZONTAL)
        dir_hbox.Add(dir_label)
        dir_hbox.Add(self.dir_text, proportion=1)
        dir_hbox.Add(self.dir_button)
        vbox.Add(dir_hbox, flag=wx.ALL|wx.EXPAND, border=10)
        
        
        
        
        ## Add space in panel.
        vbox.Add(wx.StaticText(panel, label=""))
        
        ###############
        ## Run Button
        ###############
        
        ## Add run button.
        self.run_button = wx.ToggleButton(panel, label = "Run")
        self.run_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnToggle)
        vbox.Add(self.run_button, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, border=10)
        
          
        panel.SetSizer(vbox)        
        panel.Fit()
        
        self.panel = panel
        
        
        #self.Fit()
        self.Centre()
        self.Show()
        
        
        
        

        
        
        
        
        
    def OnOpen(self, event):
        
        if self.S8_thread != None:
            message = "Cannot open a new file while one is being processed."
            msg_dlg = wx.MessageDialog(None, message, "Warning", wx.OK | wx.ICON_EXCLAMATION)
            msg_dlg.ShowModal()
            msg_dlg.Destroy()
            return
        
        
        dlg = wx.FileDialog(None,
                            message = "Choose one or more RAW file(s)",
                            wildcard = "RAW files (*.raw)|*.raw",
                            style = wx.FD_OPEN | wx.FD_CHANGE_DIR | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE)
        
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            dlg.Destroy()
            if not all([re.match(r".*.raw", path) for path in paths]):
                self.GUI_error_printer("At least one selected file is not a .raw file.")
                return
            else:
                if not self.RAW_Open_and_Count(paths[0]):
                    return
                if not self.Determine_Max_Injection_Time():
                    return
        else:
            dlg.Destroy()
            return
        
        
        self.current_raw_file_index = 0
        
        num_of_mass_ranges = self.num_of_mass_ranges
        mass_ranges = self.mass_ranges
        
        
        for child in self.segment_settings_panel.GetChildren():
            child.Destroy()
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        ##############
        ## Algorithm Settings For Each Segment
        ##############
        sizer = wx.FlexGridSizer(len(self.required_columns)+1, num_of_mass_ranges+1, 10, 10)
        
        st1 = wx.StaticText(self.segment_settings_panel, label="")
        sizer.Add(st1)
        for mass_range in mass_ranges:
            label = "Segment " + str(mass_range)
            st = wx.StaticText(self.segment_settings_panel, label=label)
            sizer.Add(st)
        
        ## Reset the dictionary that holds the values for the segment settings.
        self.segment_settings = {}
        for setting in self.required_columns:
            self.segment_settings[setting] = {}
        
        ## Create a text control for each segment setting.
        self.segment_text_ctrls = {}
        for setting in self.required_columns:
            sizer.Add(wx.StaticText(self.segment_settings_panel, label = setting + ":"))
            self.segment_text_ctrls[setting] = {}
            for mass_range in mass_ranges:
                self.segment_text_ctrls[setting][mass_range] = wx.TextCtrl(self.segment_settings_panel)
                self.segment_text_ctrls[setting][mass_range].setting = setting
                self.segment_text_ctrls[setting][mass_range].mass_range = mass_range
                self.segment_text_ctrls[setting][mass_range].Bind(wx.EVT_TEXT, self.On_Setting_Change)
                if setting == "Max Injection Time":
                    self.segment_text_ctrls[setting][mass_range].SetValue(str(self.max_inj_times[mass_range]))
                else:
                    self.segment_text_ctrls[setting][mass_range].SetValue("")
                sizer.Add(self.segment_text_ctrls[setting][mass_range])

            
        
        ## Make the grid resizeable.
        for i in range(0, num_of_mass_ranges):
            sizer.AddGrowableCol(i, 1)
        
        
        
        vbox.Add(sizer, proportion=1, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, border=10 )
        
        
        ###############
        ## Load and Save Config File
        ###############
        
        ## Create the buttons to load and save segment info.    
        self.load_button = wx.Button(self.segment_settings_panel, label="Load From File")
        self.load_button.Bind(wx.EVT_BUTTON, self.Load_From_File)
        self.save_button = wx.Button(self.segment_settings_panel, label = "Save To File")
        self.save_button.Bind(wx.EVT_BUTTON, self.Save_To_File)
        
        ## Add buttons to the panel.
        hbox_tolerance_buttons = wx.BoxSizer(wx.HORIZONTAL)
        hbox_tolerance_buttons.Add(self.load_button)
        hbox_tolerance_buttons.Add(self.save_button)
        vbox.Add(hbox_tolerance_buttons, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, border=10)
        
        
        
        
        
        self.segment_settings_panel.SetSizer(vbox)
        self.segment_settings_panel.Fit()
        
        self.current_raw_file_txt_ctrl.Clear()
        label = " \n".join(paths)
        self.current_raw_file_txt_ctrl.write(label)
        
        self.paths = paths
        
        self.panel.Fit()
        self.Fit()
        self.Centre()
        
        




    def OnAbout(self, event):
        message = "Author: Travis Thompson \n\n" + \
        "Creation Date: September 2017 \n" + \
        "Version Number: " + __version__
       
        dlg = GMD.GenericMessageDialog(None, message, "About", wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()





    def OnSettingsDescription(self, event):
        message = "Internal Calibration Masses - The m/z of the compound used to correct the m/z's of the segment. Scans without this mass detected in them are also dropped from the analysis. This is typically an internal standard or reliable contaminate with a known m/z. To forgo this correction and the dropping of scans set the value to 0.0 for every segment that it is not desired for.\n\n" + \
            "Internal Calibration Tolerances - The allowable difference between the Internal Calibration Mass and the closest m/z to it. The closest m/z to the Internal Calibration Mass will be found, and if the m/z difference is too large then the analysis will stop because the Internal Calibration Mass could not be found. This should be as small as possible to consistently find the Internal Calibration Mass without being so large that an incorrect m/z is identified as the Internal Calibrant.\n\n" + \
            "Overlap Tolerances - Exactly like the Internal Calibration Tolerances except the m/z's that are being matched are the m/z's in the overlapping regions between segments. All of the peaks in the overlapping region of one segment are matched to the closest peaks in the overlapping region of the other segment. If the difference between the m/z's is greater than the Overlap Tolerance then the peaks aren't considered to be the same compound and the match is ignored.\n\n" + \
            "Bin Widths - The multiple of 32-bit binary floating point precision to group m/z's with. The m/z's across all scans within a segment are arranged from least to greatest and the difference between the m/z and the next closest m/z is computed. For example, 150.0001, 150.0002, 150.0004 would have nearest neighbor differences of .0001 and .0002. Due to how computers store decimal numbers the nearest neighbor differences will be a multiple of the smallest unit of binary precision. A new grouping will be made whenever the nearest neighbor difference exceeds the Bin Width. See the paper for more details, but the recommended setting is 5.\n\n" + \
            "Max Injection Time - The maximum injection time a scan can have before it is dropped from the analysis. If a scan has an injection time greater than or equal to this value then it is dropped from analysis. By default this is set to the maximum injection time of the scan settings in the .raw file. Typically if a scan has an injection time equal to the maximum injection time for the scan then it was a bad scan, so it is dropped from analysis. To not drop any scans simply set this to a value higher than the maximum injection time for the segment.\n\n" + \
            "Repetition Rate Filter - The minimum percentage of scans a peak must appear in to be considered a real peak. If a peak has m/z's less than this percentage then it considered noise and dropped. For example, if a grouping of m/z's has 5 m/z's but there were 50 total scans then the grouping has a repetition rate of 10%, so if the Repetition Rate Filter is set to greater than 10 the grouping will be dropped.\n\n" + \
            "Front Segments To Skip - The number of segments to skip over in analysis, starting from the front. For example, if the m/z segments are (150 - 300), (270 - 600), and (570 - 900) the first one or two segments could be ignored and the S7 analysis can be done ignoring the (150 - 300) segment or the (150 - 300) and (150 - 300) and (270 - 600) segments. Segments in the middle cannot be ignored since there will be no overlapping region between the (150 - 300) and (570 - 900) segments.\n\n" + \
            "Rear Segments To Skip - The same as Front Segments To Skip, but starting from the rear."


        dlg = wx.lib.dialogs.ScrolledMessageDialog(None, message, "Settings Description")
        dlg.ShowModal()
        dlg.Destroy()






    def On_Setting_Change(self, event):
        control = event.GetEventObject()
        
        control_value = control.GetValue()
        setting_name = control.setting
        mass_range = control.mass_range
        
        self.segment_settings[control.setting][control.mass_range] = control_value
        
        if control_value != "":
        
            try:
                float_value = float(control_value)
            except ValueError:
                
                control.SetBackgroundColour("pink")
                control.Refresh()
                            
                self.GUI_error_printer("Error: Invalid entry for " + setting_name + " in segment " + str(mass_range) + ".")
            else:
                if numpy.isnan(float_value) or numpy.isinf(float_value):
                    
                    control.SetBackgroundColour("pink")
                    control.Refresh()
                    
                    self.GUI_error_printer("Error: Invalid number for " + setting_name + " in segment " + str(mass_range) + ".")
                
                elif setting_name == "Internal Calibration Masses" and \
                (float_value > mass_range[1] or float_value < mass_range[0]) and \
                float_value != 0.0:
                    
                    control.SetBackgroundColour("pink")
                    control.Refresh()
                    
                    self.GUI_error_printer("Error: Internal Calibration Mass outside of mass range for segment " + str(mass_range) + ".")
                    
                elif setting_name == "Max Injection Time" and (float_value <= 0):
                    
                    control.SetBackgroundColour("pink")
                    control.Refresh()
                    
                    self.GUI_error_printer("Error: Max Injection Time outside of range for segment " + str(mass_range) + ".")
                
                else:
                    control.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
                    control.Refresh()
                



    def OnRepetitionRate(self, event):
        
        repetition_rate_okay = False
        
        is_event = event != None
        
        
        if self.GUI_mode:
            if is_event:
                text_ctrl = event.GetEventObject()
            else:
                text_ctrl = self.repetition_rate_txt_control
                
            text = text_ctrl.GetValue()
        else:
            text = self.parser_args.repetition_rate
        
        
        try:
            float(text)
        except ValueError:
            if self.GUI_mode:
                text_ctrl.SetBackgroundColour("pink")
                text_ctrl.Refresh()
            
            if is_event:
                self.GUI_error_printer("Error: Please enter a valid number for Repetition Rate Filter.")
            else:
                self.Algorithm_error_printer("Error: Invalid entry entered for Repetition Rate Filter.")
        else:
            if numpy.isnan(float(text)) or numpy.isinf(float(text)):
                
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour("pink")
                    text_ctrl.Refresh()
                
                if is_event:
                    self.GUI_error_printer("Error: Please enter a valid number for Repetition Rate Filter.")
                else:
                    self.Algorithm_error_printer("Error: Invalid number entered for Repetition Rate Filter.")
            
            elif (float(text) >= 100 or float(text) < 0):
                
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour("pink")
                    text_ctrl.Refresh()
                
                if is_event:
                    self.GUI_error_printer("Error: Repetition Rate Filter outside of range.")
                else:
                    self.Algorithm_error_printer("Error: Repetition Rate Filter outside of range.")
            
            else:
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
                    text_ctrl.Refresh()
                
                self.repetition_rate = float(text)
                repetition_rate_okay = True
                    
        return repetition_rate_okay




    def OnFrontSkip(self, event):
        
        front_skip_okay = False
        
        is_event = event != None
        
                
        if self.GUI_mode:
            if is_event:
                text_ctrl = event.GetEventObject()
            else:
                text_ctrl = self.front_skip_txt_control
                
            text = text_ctrl.GetValue()
        else:
            text = self.parser_args.front_skip
        
        
        try:
            int(text)
        except ValueError:
            if self.GUI_mode:
                text_ctrl.SetBackgroundColour("pink")
                text_ctrl.Refresh()
            
            if is_event:
                self.GUI_error_printer("Error: Please enter a valid number for front segments to skip.")
            else:
                self.Algorithm_error_printer("Error: Ivalid entry for front segments to skip.")
        else:
            if numpy.isnan(int(text)) or numpy.isinf(int(text)):
                
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour("pink")
                    text_ctrl.Refresh()
                
                if is_event:
                    self.GUI_error_printer("Error: Please enter a valid number for front segments to skip.")
                else:
                    self.Algorithm_error_printer("Error: Ivalid number for front segments to skip.")
            
            elif (int(text) > self.num_of_mass_ranges or int(text) < 0):
                
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour("pink")
                    text_ctrl.Refresh()
                
                if is_event:
                    self.GUI_error_printer("Error: Front skip segments outside of segment range.")
                else:
                    self.Algorithm_error_printer("Error: Front skip segments outside of segment range.")
            
            else:
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
                    text_ctrl.Refresh()
                
                self.front_segments_to_skip = int(text)
                front_skip_okay = True
                    
        return front_skip_okay
                    
                    
                    
                    
                    
                    
    def OnRearSkip(self, event):
        
        rear_skip_okay = False
        
        is_event = event != None
        
        
        if self.GUI_mode:
            if is_event:
                text_ctrl = event.GetEventObject()
            else:
                text_ctrl = self.rear_skip_txt_control
                
            text = text_ctrl.GetValue()
        else:
            text = self.parser_args.rear_skip
        
        
        try:
            int(text)
        except ValueError:
            if self.GUI_mode:
                text_ctrl.SetBackgroundColour("pink")
                text_ctrl.Refresh()
            
            if is_event:
                self.GUI_error_printer("Error: Please enter a valid number for rear segments to skip.")
            else:
                self.Algorithm_error_printer("Error: Invalid entry for rear segments to skip.")
        else:
            if numpy.isnan(int(text)) or numpy.isinf(int(text)):
                
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour("pink")
                    text_ctrl.Refresh()
                
                if is_event:
                    self.GUI_error_printer("Error: Please enter a valid number for rear segments to skip.")
                else:
                    self.Algorithm_error_printer("Error: Invalid number for rear segments to skip.")
            
            elif (int(text) > self.num_of_mass_ranges or int(text) < 0):
                
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour("pink")
                    text_ctrl.Refresh()
                
                if is_event:
                    self.GUI_error_printer("Error: Rear skip segments outside of segment range.")
                else:
                    self.Algorithm_error_printer("Error: Rear skip segments outside of segment range.")
            
            else:
                if self.GUI_mode:
                    text_ctrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
                    text_ctrl.Refresh()
                    
                self.rear_segments_to_skip = int(text)
                rear_skip_okay = True
                    
        return rear_skip_okay




    def OnQuit(self, e):
        self.Close()





    def Set_Output_Flags(self, event):
        box = event.GetEventObject()
        self.flags[box.GetLabel()] = box.GetValue()




    def RAW_Open_and_Count(self, path):
        if self.rawfile != None:
            self.rawfile.Close()
            self.rawfile = None
            
        try:
            rawfile = ThermoRawfile(path)
        except IOError:
            self.Algorithm_error_printer("Error: Unable to open RAWfile, check the path and try again.")
            return False
        
        mass_ranges = []
        num_of_scans = rawfile.GetLastSpectrumNumber()
        mass_range = rawfile.GetMassRangeFromScanNum(1, 0)
        mass_ranges.append(mass_range)
        for i in range(1,num_of_scans+1):
                
            mass_range = rawfile.GetMassRangeFromScanNum(i, 0)
            
            if not(mass_range in mass_ranges):
                mass_ranges.append(mass_range)
        
        ## Order mass_ranges from least to greatest.
        mass_ranges.sort()
         
        num_of_mass_ranges = len(mass_ranges)
        
        self.rawfile = rawfile
        self.raw_file_name = os.path.split(path)[1]
        self.num_of_mass_ranges = num_of_mass_ranges
        self.mass_ranges = mass_ranges
        self.num_of_scans = num_of_scans
        return True
        
        
        
        
    def Determine_Max_Injection_Time(self):
        
        rawfile = self.rawfile
        mass_ranges = self.mass_ranges
        
        ## Find maximum ion injection time.
        temp = rawfile.GetInstMethod()
        max_inj_times = {}
        for i in range(len(mass_ranges)):
            ## Assuming the instrument method information is always stored in the same way this pattern should find the max
            ## injection time for each segment.
            try:
                max_inj_time = re.search(r"Scan Range \(m/z\) = " + str(int(mass_ranges[i][0])) + "-" + str(int(mass_ranges[i][1])) + r"\n\t\t\tMaximum Injection Time \(ms\) = \d*", temp).group()
                ## If the pattern was found then the injection time should be the 3rd number.
                max_inj_time = int(re.findall(r"\d+", max_inj_time)[2])
                max_inj_times[mass_ranges[i]] = max_inj_time
            except:
                self.Algorithm_error_printer("Error: Could not find Maximum Injection Time for the segment " + str(mass_ranges[i]))
                return False
            
            
        self.max_inj_times = max_inj_times
        return True
    



    def Load_From_File(self, event):
        if self.GUI_mode:
            ## Open a window to select the csv file.
            path = ""
            dlg = wx.FileDialog(None,
                                message = "Choose a CSV File",
                                wildcard = "CSV files (*.csv)|*.csv",
                                style = wx.FD_OPEN | wx.FD_CHANGE_DIR | wx.FD_FILE_MUST_EXIST)
        
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return False
        else:
            path = self.parser_args.settings_file
    
        ## Read in the csv and put the values in the textctrls.
        try:
            df = pandas.read_csv(path)
        except IOError:
            self.GUI_error_printer("Error: Settings file does not exist.")
            return False
        else:
            if set(self.required_columns).issubset(df.columns):
                                
                if self.GUI_mode:
                    for setting in self.required_columns:
                        for i in range(0, len(df)):
                            if i < len(self.mass_ranges):
                                self.segment_text_ctrls[setting][self.mass_ranges[i]].SetValue(df.loc[i, setting].round(12).astype(str))
                                
                    self.ValidateSegmentSettings()
                    return True
                else:
                    if len(df) != len(self.mass_ranges):
                        self.GUI_error_printer("Error: The input csv file does not have the same number of rows as there are mass ranges in the RAW file.")
                        return False
                    else:
                        for setting in self.required_columns:
                            for i in range(0, len((df))):
                                self.segment_settings[setting][self.mass_ranges[i]] = df.loc[i, setting]
                        
                        return True
                        
            else:
                self.GUI_error_printer("Error: Selected csv file does not have the correct columns.")
                return False
 
    
    
    
    def Save_To_File(self, event):
        if self.ValidateSegmentSettings():    
            ## Open a window for the user to select where to save.
            dlg = wx.FileDialog(None,
                        message = "Save file as ...",
                        wildcard = "CSV files (*.csv)|*.csv",
                        style = wx.FD_SAVE | wx.FD_CHANGE_DIR)
    
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
    
            ## Write the txtctrl values out to a csv file.
            df = pandas.DataFrame(columns = self.required_columns, index=self.mass_ranges)
            for setting_name, setting_values in self.segment_settings.items():
                for mass_range, value in setting_values.items():
                    df.loc[mass_range, setting_name] = value
                
            df.to_csv(path, index=False)
    
    
    
    
    def ValidateSegmentSettings(self):
        all_clear = True
        
        ## Only change the txt control color if it is the first raw file in the batch.
        change_txt_cntrls = self.current_raw_file_index == 0
        
        for setting_name, setting_values in self.segment_settings.items():
            for mass_range, value in setting_values.items():
                text = value
                if self.GUI_mode:
                    text_ctrl = self.segment_text_ctrls[setting_name][mass_range]
                try:
                    float(text)
                except ValueError:
                    
                    if self.GUI_mode and change_txt_cntrls:
                        text_ctrl.SetBackgroundColour("pink")
                        text_ctrl.Refresh()
                    
                    self.Algorithm_error_printer("Error: Invalid entry for " + setting_name + " in segment " + str(mass_range) + ".")
                    all_clear = False
                else:
                    if numpy.isnan(float(text)) or numpy.isinf(float(text)):
                        
                        if self.GUI_mode and change_txt_cntrls:
                            text_ctrl.SetBackgroundColour("pink")
                            text_ctrl.Refresh()
                        
                        self.Algorithm_error_printer("Error: Invalid number for " + setting_name + " in segment " + str(mass_range) + ".")
                        all_clear = False
                    
                    elif setting_name == "Internal Calibration Masses" and \
                    (float(text) > mass_range[1] or float(text) < mass_range[0]) and\
                    float(text) != 0.0:
                        
                        if self.GUI_mode and change_txt_cntrls:
                            text_ctrl.SetBackgroundColour("pink")
                            text_ctrl.Refresh()
                        
                        self.Algorithm_error_printer("Error: Internal Calibration Mass outside of mass range for segment " + str(mass_range) + ".")
                        all_clear = False
                        
                    elif setting_name == "Max Injection Time" and (float(text) <= 0):
                        
                        if self.GUI_mode and change_txt_cntrls:
                            text_ctrl.SetBackgroundColour("pink")
                            text_ctrl.Refresh()
                        
                        self.Algorithm_error_printer("Error: Max Injection Time outside of range for segment " + str(mass_range) + ".")
                        all_clear = False
                    
                    else:
                        if self.GUI_mode and change_txt_cntrls:
                            text_ctrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
                            text_ctrl.Refresh()
                        
                        self.segment_settings[setting_name][mass_range] = float(text)
                        
        return all_clear
            
            


    def Dir_Select(self, event):
        dlg = wx.DirDialog(None, "Choose a folder:", style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST | wx.DD_CHANGE_DIR)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.output_directory = dlg.GetPath()
            self.dir_text.SetValue(self.output_directory)
            dlg.Destroy()
        else:
            dlg.Destroy()
            return
        
        
        
        
    def CheckOutputDirectory(self):
        
        if self.output_directory == None:
            if self.GUI_mode:
                self.Algorithm_error_printer("Error: No output directory indicated \n")
                
                message = "Please select an output directory."
                msg_dlg = wx.MessageDialog(None, message, "Warning", wx.OK | wx.ICON_EXCLAMATION)
                msg_dlg.ShowModal()
                msg_dlg.Destroy()
            else:
                self.GUI_error_printer("Error: No output directory indicated.")
            return False
        
        if not os.path.exists(self.output_directory):
            try:
                os.makedirs(self.output_directory)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    self.GUI_error_printer("Error while creating output directory.")
                    return False
            else:
                return True
        else:
            return True

   
    
    
    
    
    def DisableImportantControls(self):
        
        ## Disable setting text controls such as Overlap Tolerances.
        for setting_name, mass_range_dict in self.segment_text_ctrls.items():
            for mass_range_name, text_ctrl in mass_range_dict.items():
                text_ctrl.Disable()
                
                
        ## Disable flag checkboxes.
        for flag_name, flag_checkbox in self.flag_checkboxes.items():
            flag_checkbox.Disable()
            
        ## Disable repetition rate text control.
        self.repetition_rate_txt_control.Disable()
        
        
        ## Disable front and rear skip text controls.
        self.front_skip_txt_control.Disable()
        self.rear_skip_txt_control.Disable()
        
        
        ## Disable Load and Save buttons.
        self.load_button.Disable()
        self.save_button.Disable()
        
        ## Disable selecting output directory button.
        self.dir_button.Disable()
        
        
        
        
        
        
    def EnableImportantControls(self):
        
        ## Enable setting text controls such as Overlap Tolerances.
        for setting_name, mass_range_dict in self.segment_text_ctrls.items():
            for mass_range_name, text_ctrl in mass_range_dict.items():
                text_ctrl.Enable()
                
                
        ## Enable flag checkboxes.
        for flag_name, flag_checkbox in self.flag_checkboxes.items():
            flag_checkbox.Enable()
            
        
        ## Enable repetition rate text control.
        self.repetition_rate_txt_control.Enable()
        
        
        ## Enable front and rear skip text controls.
        self.front_skip_txt_control.Enable()
        self.rear_skip_txt_control.Enable()
        
        
        ## Enable Load and Save buttons.
        self.load_button.Enable()
        self.save_button.Enable()
        
        ## Enable selecting output directory button.
        self.dir_button.Enable()
    
    
    
    
    
    def CheckReportOutputs(self):
        
        if not any([flag_value for flag_name, flag_value in self.flags.items()]):
            if self.GUI_mode:
                self.GUI_error_printer("At least one report output must be selected.")
            
            self.Algorithm_error_printer("Error: At least one report output must be selected \n")
            return False
        
        else:
            return True
            
            
    
    
    
    
    
    def OnToggle(self, event):
        button = event.GetEventObject()
        
        ## Clicking run will always start the batch over.
        self.current_raw_file_index = 0
        
        ## Clear any previous error messages.
        self.statusbar.SetStatusText("")
    
    
        if button.GetValue() == True:
            
            if self.paths == None:
                self.GUI_error_printer("No RAW file selected.")
                button.SetValue(False)
                return
            
            
            ## Make sure the settings are okay, an output directory is chosen, and at least one report is selected before proceeding.
            if self.ValidateCurrentRawFile() and self.CheckOutputDirectory() and self.CheckReportOutputs():
            
                if not self.S8_thread:
                    self.S8_thread = S8Thread(self)
                    self.Algorithm_error_printer("Starting run on file " + self.raw_file_name + ":")
                    self.Algorithm_error_printer('Timestamp: {:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now()))
                    button.SetLabel("Cancel")
                    
                    self.DisableImportantControls()
                    
                else:
                    message = "Thread already started, something seriously wrong happened."
                    msg_dlg = wx.MessageDialog(None, message, "Warning", wx.OK | wx.ICON_EXCLAMATION)
                    msg_dlg.ShowModal()
                    msg_dlg.Destroy()
                    
            else:
                button.SetValue(False)
                return
                    
        
        else:
            
            if self.S8_thread:
                self.S8_thread.abort()
            else:
                message = "Thread already aborted, something seriously wrong happened."
                msg_dlg = wx.MessageDialog(None, message, "Warning", wx.OK | wx.ICON_EXCLAMATION)
                msg_dlg.ShowModal()
                msg_dlg.Destroy()


    
    
    
    def thread_aborted(self, event):
        
        if event.message != None:
            self.Algorithm_error_printer(event.message)
            dlg_message = "Run Aborted! \n\n" + event.message
        else:
            dlg_message = "Run Aborted!" 
        
        self.Algorithm_error_printer("Run Aborted")
        self.Algorithm_error_printer('Timestamp: {:%Y-%b-%d %H:%M:%S} \n'.format(datetime.datetime.now()))
        
        self.S8_thread = None
        
        if self.rawfile != None:
            self.rawfile.Close()
            self.rawfile = None
        
        if self.GUI_mode and (self.current_raw_file_index == len(self.paths) - 1 or event.message == "User Aborted"):
            msg_dlg = wx.MessageDialog(None, dlg_message, "Run Aborted", wx.OK | wx.ICON_EXCLAMATION)
            msg_dlg.ShowModal()
            msg_dlg.Destroy()
        
            self.run_button.SetValue(False)
            self.run_button.SetLabel("Run")
            
            self.EnableImportantControls()
        
            
        elif self.FindNextValidRawFile():    
            
            self.S8_thread = S8Thread(self)
            self.Algorithm_error_printer("Starting run on file " + self.raw_file_name + ":")
            self.Algorithm_error_printer('Timestamp: {:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now()))
            
        else:
            exit(1)
        
        
    
    
    
    def thread_error(self, event):
        
        if event.message != None:
            self.Algorithm_error_printer(event.message)
        



    def thread_completed(self, event):
                
        self.Algorithm_error_printer("Run Completed")
        self.Algorithm_error_printer('Timestamp: {:%Y-%b-%d %H:%M:%S} \n'.format(datetime.datetime.now()))
        
        self.S8_thread = None 
        
        if self.rawfile != None:
            self.rawfile.Close()
            self.rawfile = None
        
        if self.GUI_mode and self.current_raw_file_index == len(self.paths) - 1:
            message = "Run Completed!"
            msg_dlg = wx.MessageDialog(None, message, "Run Complete", wx.OK | wx.ICON_NONE)
            msg_dlg.ShowModal()
            msg_dlg.Destroy()
        
            self.run_button.SetValue(False)
            self.run_button.SetLabel("Run")
            
            self.EnableImportantControls()
            
        
        elif self.FindNextValidRawFile():    
            
            self.S8_thread = S8Thread(self)
            self.Algorithm_error_printer("Starting run on file " + self.raw_file_name + ":")
            self.Algorithm_error_printer('Timestamp: {:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now()))
            
        else:
            exit(0)
               
        
        

    

    def Algorithm_error_printer(self, message):
#        try:
#            ## write will append text, Clear() will clear. I think.
#            self.error_frame.error_txt_ctrl.write(message + "\n")
#        except NameError:
#            print(message)
            
        if self.GUI_mode:
            ## write will append text, flush will clear. I think.
            self.error_frame.error_txt_ctrl.write(message + "\n")
        else:
            print(message)
            if re.match(r"Error.*", message):
                print("Aborting run.")
                exit(1)
    
        
        
    def GUI_error_printer(self, message):
#        try:
#            self.statusbar.SetStatusText(message)
#        except NameError:
#            print(message)
        
        if self.GUI_mode:
            self.statusbar.SetStatusText(message)
        else:
            print(message)
            if re.match(r"Error.*", message):
                print("Aborting run.")
            exit(1)
                




    def ValidateCurrentRawFile(self):
        all_clear = True
        
        raw_file_name = os.path.split(self.paths[self.current_raw_file_index])[1]
        self.Algorithm_error_printer("Validating rawfile " + raw_file_name + ":")
        
        if not self.RAW_Open_and_Count(self.paths[self.current_raw_file_index]):
            all_clear = False
        if not self.Determine_Max_Injection_Time():
            all_clear = False
        
        
        front_skip_okay = self.OnFrontSkip(None)
        rear_skip_okay = self.OnRearSkip(None)
        
        ## Make sure the front and rear number of segments to skip have valid values and that they don't skip too many segments.
        if rear_skip_okay and front_skip_okay and (self.rear_segments_to_skip + self.front_segments_to_skip) > self.num_of_mass_ranges:
            message = "There are more segments to skip than there are total segments. /nDecrease the number of segments to skip."
            self.Algorithm_error_printer(message)
            all_clear = False
        
        
        if not rear_skip_okay :
            all_clear = False
        
        if not front_skip_okay :
            all_clear = False
        
        ## Make sure the repetition rate is okay.
        if not self.OnRepetitionRate(None) :
            all_clear = False
        
        if not self.ValidateSegmentSettings():
            all_clear = False
            
        return all_clear
        





    def FindNextValidRawFile(self):
        """Increment current_raw_file_index until a valid raw file is found or we run out of paths to check.
           Return False if run out of paths and True otherwise."""
        while True:
            if self.current_raw_file_index >= len(self.paths) - 1:
                return False
            
            self.current_raw_file_index += 1
            
            if self.ValidateCurrentRawFile():
                return True
        
            
            
        


#############################    ARGPARSER DEFS   #######################################
## Modified from: 
## https://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin
def repetition_rate(x):
    x = float(x)
    if x < 0.0 or x >= 100.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 100.0)"%(x,))
    return x



def raw_file(x):
    x = os.path.abspath(x)
    if not (re.match(r".*\.raw", x) and os.path.isfile(x)):
        raise argparse.ArgumentTypeError("Invalid entry for raw_file.")
    return x


def settings_file(x):
    x = os.path.abspath(x)
    if not (re.match(r".*\.csv", x) and os.path.isfile(x)):
        raise argparse.ArgumentTypeError("Invalid entry for settings_file.")
    return x


def output_directory(x):
    x = os.path.abspath(x)
    return x



    


#############################    MAIN DEF   #######################################

def main():
    
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser()
        parser.add_argument("raw_file", type=raw_file, help="Path to the Thermo RAW file.")
        parser.add_argument("settings_file", type=settings_file, help="Path to the .csv file with the settings for each segment in the RAW file.")
        parser.add_argument("output_directory", type=output_directory, help="Path to the directory where the output files will be saved. \
                            The directory is created if it does not exist.")
        parser.add_argument("-r", "--repetition_rate", type=repetition_rate, default=10.0, help="Value used for the repetition rate filter. \
                            Default = 10.0. Any signals that show up less than this rate will be dropped. \
                            Ex. If a signal is in 1 out of 150 scans its repetition rate is 1/150*100 = 0.6666 \
                            which is less than 10, so it would be dropped.")
        parser.add_argument("-fs", "--front_skip", type=int, default=0, help="The number of segments counting from the lowest mass range to skip over in analysis.\
                            Ex. If there are 3 mass ranges (150, 300), (270,500), (470,800) the first 2 ranges can be skipped by entering values of 1 or 2.")
        parser.add_argument("-rs", "--rear_skip", type=int, default=0, help="The number of segments counting from the highest mass range to skip over in analysis.\
                            Ex. If there are 3 mass ranges (150, 300), (270,500), (470,800) the last 2 ranges can be skipped by entering values of 1 or 2.")
        parser.add_argument("-fd", "--full_data", action="store_true", help="If selected prints out the \"full data\" report.")
        parser.add_argument("-nd", "--noise_data", action="store_true", help="If selected prints out the \"noise data\" report.")
        parser.add_argument("-pl", "--peak_list", action="store_false", help="If selected does not print out the \"peak list\" report.")
        
        args = parser.parse_args()
        
        if not (args.full_data or args.noise_data or args.peak_list):
            parser.error("At least one report type must be selected as an output.")


    
        ex = wx.App(False)
        S8_GUI(None, title = "S8 Analysis", GUI_show = False, parser_args = args)
        ex.MainLoop()
        
    else:
        ex = wx.App(False)
        S8_GUI(None, title = "S8 Analysis")
        ex.MainLoop() 
        
        
    




if __name__ == "__main__":
    main()
              
                
