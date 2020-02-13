Distinguishing between Flaring and Non-Flaring Active Regions- A Machine Learning Perspectiver

In this study, our aim is to distinguish flaring active region from the non-flaring one. Active region which generates a flare during their transit over visible solar disk is called flaring region. Active region which does not generate any flare during their transit is known as non-flaring region.

We have used SDO/HMI active region magnetogram data for this study. One can download the data from JSOC website. One can develop their own program to create the database of magnetic parameters (one can see our python script for generation of magnetic parameters) or, one can download magnetic parameter keywords from JSOC directly. HMI team provides us directly the magnetic parameter keywords which can also be used for the study.  


First, we have to generate a database for flaring and non-flaring active region. One can use flare_noflare_database.py for this purpose.
