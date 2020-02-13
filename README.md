Distinguishing between Flaring and Non-Flaring Active Regions- A Machine Learning Perspectiver

In this study, our aim is to distinguish flaring active region from the non-flaring one. Active region which generates a flare during their transit over visible solar disk is called flaring region. Active region which does not generate any flare during their transit is known as non-flaring region.

We have used SDO/HMI active region magnetogram data for this study. One can download the data from JSOC website. One can develop their own program to create the database of magnetic parameters (one can see our python script for generation of magnetic parameters) or, one can download magnetic parameter keywords from JSOC directly. HMI team provides us directly the magnetic parameter keywords which can also be used for the study. One can download SHARPs data keyword using drms client.  


First, we have to generate a database for flaring and non-flaring active region. One can use flare_noflare_database.py for this purpose.

Second, we have to go back 24 hr, and we have to cut a 12 hr window. For this purpose we have used final_positive.m for flaring active region; and final_negative.m for non-flaring active region. We use MATLAB programming language for that purpose. If non-flaring active region exists for 14 days in the solar surface. We use day number 7 as a reference., i.e., middle point.

We use positive_zero_span_data.m for flaring active regions and negative_zero_span_data.m for non-flaring active region when we want to generate the data with 24 hour loopback but zero span.
