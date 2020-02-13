Distinguishing between Flaring and Non-Flaring Active Regions- A Machine Learning Perspectiver

In this study, our aim is to distinguish flaring active region from the non-flaring one. Active region which generates a flare during their transit over visible solar disk is called flaring region. Active region which does not generate any flare during their transit is known as non-flaring region.

We have used SDO/HMI active region magnetogram data for this study. One can download the data from JSOC website. One can develop their own program to create the database of magnetic parameters (one can see our python script for generation of magnetic parameters) or, one can download magnetic parameter keywords from JSOC directly. HMI team provides us directly the magnetic parameter keywords which can also be used for the study. One can download SHARPs data keyword using drms client.  


First, we have to generate a database for flaring and non-flaring active region. One can use flare_noflare_database.py for this purpose. Our database also only considers the data which is within +/- 72 degree.

Second, we have to go back 24 hr, and we have to cut a 12 hr window. For this purpose we have used final_positive_data_different_span.m for flaring active region; and final_negative_data_different_span.m for non-flaring active region. We use MATLAB programming language for that purpose. If non-flaring active region exists for 14 days in the solar surface. We use day number 7 as a reference., i.e., middle point.

We use positive_zero_span_data.m for flaring active regions and negative_zero_span_data.m for non-flaring active region when we want to generate the data with 24 hour loopback but zero span.

Third, we saved all .csv generated files in a folder. Then we use some R-program to remove the empty .csv files if any. Finally, we listed all non-empty .csv files into a .txt file.

R program for removing empty .csv files in a directory:

library(R.utils)

lapply(Filter(function(x) countLines(x)<=3, list.files(pattern='.csv')), unlink)

Using this two command line in R, we can remove all empty .csv files from a directory. It is necessary for this project because we have to list all non-empty .csv files of a directory into a .txt file.

Finally, we have to list all non empty .csv files of the directory into a .txt file.

Command: ls | grep '.csv' > sax.txt
