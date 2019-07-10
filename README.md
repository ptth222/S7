S7 (Sample-Specific Silhouette Segmented-Scan Spectral Stitching)
=============

S7 is an implementation of the paper Improved Segmented-Scan Spectral Stitching for Stable Isotope Resolved Metabolomics (SIRM) by Ultra-High-Resolution Fourier Transform Mass Spectrometry.
It is a program to normalize segment intensities by injection time, group m/z's across scans, correct m/z's based on internal calibrants, reduce the number of noise peaks, and normalize intensities between segments.

Scans without any internal calibrant are dropped. The closest m/z to the Internal Calibration Mass is found in each scan, and if the difference between it and the Internal Calibration Mass is greater than the Internal Calibration Tolerance the scan is dropped from the analysis.

Intensities are normalized by the ratio of the segment's highest injection time with each scan's injection time. For example, if the highest injection time for the scans in a segment is 50 ms and another scan has an injection time of 40 ms then that scan's intensities are multiplied by 50/40 = 1.25. Every scan's intensities are multiplied by the ratio of the highest injection time and that scan's injection time, 50/(scan injection time).

m/z's are grouped across scans using minimum neighbor distances. Whenever a gap greater than the minimum distance is seen between adjacent m/z's a new grouping is created. The minimum neighbor distance is a multiple of the smallest unit of 32-bit digital precision for the m/z, see the paper for more detail. 

Noise peaks are removed by determining a peak's repetition rate and dropping it if it is less than the user defined repetition rate filter. The repetition rate is the percentage of scans a peak appeared in. For example if a grouping of m/z's has signals from 5 scans and the segment had a total of 50 scans then the repetition rate for the peak would be 5/50*100 = 10.

The median m/z value of the m/z grouping is used as the m/z of the peak. This peak m/z is then used to find the closest peak to the indicated internal calibrant for the segment and correct all of the peak m/z's in that segment by the difference between the closest peak and the internal calibrant. If the closest match is not within the indicated internal calibrant tolerance then a message is displayed and the analysis stopped.

Intensities are normalized between segments by comparing the ratios of the peak intensities of peaks found in the overlapping regions of adjacent segments. Peaks are matched in overlapping regions based on the indicated overlap tolerance and the ratio of the intensities for each match are computed. The average of the ratios is then used to normalize the intensities of the higher m/z segment by multiplying them by the ratio.

There are up to 3 different outputs to the program depending on user input. There is a "Full Data" output which is an Excel file with the peak information for every segment. The "Noise Reduction" output is an Excel file with the peak information for every segment but with the noise peaks removed. The file also includes information for the overlapping regions. Both the "Full Data" and "Noise Reduction" files may be split into multiple files depending on the size of the data. The "Peak List" output is an Excel file that simply lists the final corrected m/z's and intensities for each peak.

Reference
-----------
If you use the software for a publication please use DOI [10.5281/zenodo.1434958](https://zenodo.org/badge/latestdoi/150272490) as a reference.
If you would like to reference the paper please use DOI (paper DOI here).

Quick Start
-----------
 - If not already installed, install Thermo's MSFileReader.
 - Download [S7-v\*.\*.*.zip](https://zenodo.org/badge/latestdoi/150272490).
 - Unzip the file.
 - Run the installer to install the program.
 - Double click the shortcut or the exe in the installation directory to start the program in GUI mode.
 - Open one or more *.raw files from File -> Open.
 - Fill in the settings for each segment and click "Save to File" so they can be "Load from File" in the future.
 - Select the outputs you would like to see.
 - Select a directory to save the outputs in.
 - Click "Run" to start the analysis.

 A "Peak Lists" folder will be created in the output directory and all of the peak lists for each .raw file will be in this folder, but Full Data and Noise Reduction reports will be in a folder with the same name as the .raw file for each .raw file.

If you see an error like "The procedure entry point ucrtbase.terminate could not be located in the dynamic link library api-ms-win-crt-runtime-l1-1-0.dll", or "Error loading Python DLL 'C:\Users\...\python36.dll'. LoadLibrary: The specified procedure could not be found." on startup it is because Windows is missing an update. The specific update is the Universal C Runtime and can be obtained through Windows updates or [here](https://support.microsoft.com/en-us/help/2999226/update-for-universal-c-runtime-in-windows).

## **Important Notes**

Currently S7 only supports .raw files with m/z ranges from 64 to 2048, so if you are looking for compounds outside of that range the binning method will not work appropriately. If you desperately need to use this program on m/z's outside of that range contact us and we can see about extending the range of the program.

Depending on the density of peaks in the .raw file the program may take a long time to complete, and the time can be highly variable. For our data, wide scan full spectrum .raw files with only one m/z segment are the fastest and can take as little as 10 minutes to generate all reports, but the .raw files with many (10-12) m/z segments could take close to an hour to generate all reports. A bigger factor than the density of peaks, however, is what reports are generated. The "Full Data" and "Noise Reduction" reports are generally quite large, so the bottle neck in program completion time ends up being the creation of these reports. For example, the previously mentioned .raw files that would take an hour to complete if all reports were generated would only take around 10 minutes to complete if only the "Peak List" report was generated.

The program is quite memory intensive, so it is advised to have as much RAM as possible. The program has been tested to run on machines with as little as 8 GB of RAM successfully. Lower amounts of RAM have not been tested, but at least 8 GB is the recommendation.

Disclaimer
-----------
This software is supplied "as is" and "with all faults". There are no representations or warranties of any kind concerning the safety, suitability, lack of viruses, inaccuracies, typographical errors, or other harmful components of this software. The user is solely responsible for determining whether this software is compatible with your equipment and other software installed on your equipment. The user is also solely responsible for the protection of their equipment and backup of their data. We will not be liable for any damages you may suffer in connection with using, modifying, or distributing this software. We also reserve the right to make changes in the software without notification, and make no representations or warranties that the software will be suitable for use without further testing or modification. Use the software at your own risk.

License
--------
This program is a free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. See the file LICENSE for details.

Contact
--------
We can be contacted over the Issues, or directly at ptth222@uky.edu.


Using S7
=================

You can use the software as an application or run it from the command line. Both have similar functionality, but the application accepts multiple .raw files while the command line only accepts one .raw file at a time. The limitation is easily overcome with a .bat script. If running multiple .raw files from the application it is important to know that the program assumes the settings are the same for each .raw file, so if different .raw files need to use different settings they cannot be ran in a batch together.

It is important to understand the underlying assumption that the software makes on the structure of the .raw file. It is assumed that the .raw file will have different m/z segments that overlap with the next closest segments on either side. For example, (150 - 300), (270 - 600), (570 - 900), etc. If adjacent segments do not overlap then the analysis cannot be done. One exception is if there is only one segment. If there is only one segment then the program will complete as normal, but there will be no intensity normalization from overlap peak matching since there is no overlapping region.

#### **Terminology:**

segment - One m/z range within the .raw file. For example (150 - 300). 

#### **Settings Explained:**

**Internal Calibration Masses** - The m/z of the compound used to correct the m/z's of the segment. Scans without this mass detected in them are also dropped from the analysis. This is typically an internal standard or reliable contaminate with a known m/z. To forgo this correction and the dropping of scans set the value to 0.0 for every segment that it is not desired for.

**Internal Calibration Tolerances** - The allowable difference between the Internal Calibration Mass and the closest m/z to it. The closest m/z to the Internal Calibration Mass will be found, and if the m/z difference is too large then the analysis will stop because the Internal Calibration Mass could not be found. This should be as small as possible to consistently find the Internal Calibration Mass without being so large that an incorrect m/z is identified as the Internal Calibrant.

**Overlap Tolerances** - Exactly like the Internal Calibration Tolerances except the m/z's that are being matched are the m/z's in the overlapping regions between segments. All of the peaks in the overlapping region of one segment are matched to the closest peaks in the overlapping region of the other segment. If the difference between the m/z's is greater than the Overlap Tolerance then the peaks aren't considered to be the same compound and the match is ignored.

**Bin Widths** - The multiple of 32-bit binary floating point precision to group m/z's with. The m/z's across all scans within a segment are arranged from least to greatest and the difference between the m/z and the next closest m/z is computed. For example, 150.0001, 150.0002, 150.0004 would have nearest neighbor differences of .0001 and .0002. Due to how computers store decimal numbers the nearest neighbor differences will be a multiple of the smallest unit of binary precision. A new grouping will be made whenever the nearest neighbor difference exceeds the Bin Width. See the paper for more details, but the recommended setting is 5.

**Max Injection Time** - The maximum injection time a scan can have before it is dropped from the analysis. If a scan has an injection time greater than or equal to this value then it is dropped from analysis. By default this is set to the maximum injection time of the scan settings in the .raw file. Typically if a scan has an injection time equal to the maximum injection time for the scan then it was a bad scan, so it is dropped from analysis. To not drop any scans simply set this to a value higher than the maximum injection time for the segment.

**Repetition Rate Filter** - The minimum percentage of scans a peak must appear in to be considered a real peak. If a peak has m/z's less than this percentage then it considered noise and dropped. For example, if a grouping of m/z's has 5 m/z's but there were 50 total scans then the grouping has a repetition rate of 10%, so if the Repetition Rate Filter is set to greater than 10 the grouping will be dropped.

**Front Segments To Skip** - The number of segments to skip over in analysis, starting from the front. For example, if the m/z segments are (150 - 300), (270 - 600), and (570 - 900) the first one or two segments could be ignored and the S7 analysis can be done ignoring the (150 - 300) segment or the (150 - 300) and (150 - 300) and (270 - 600) segments. Segments in the middle cannot be ignored since there will be no overlapping region between the (150 - 300) and (570 - 900) segments.

**Rear Segments To Skip** - The same as Front Segments To Skip, but starting from the rear.