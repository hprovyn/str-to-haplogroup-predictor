STR TO HAPLOGROUP PREDICTOR

This project is financially supported by YSEQ and is the work of Hunter Provyn with some ideas contributed by Thomas Krahn.

WHAT IT DOES

Processes sample STR and SNP results along with a SNP phylogenetic tree to run an experiment that predicts SNP branch from a set of STRs absent any SNP info. 

REQUIREMENTS

. Python 3.6 (2 may work)#Python dependencies
  . packages numpy, pandas, time, sys, sklearn, random, json, os, operator, shutil

. CSV file containing STR and SNP results where each line contains a single ID, MARKER, ALLELE
. JSON file containing YFull tree

EXPERIMENT PROCESS

There are two shell files:

A) runExperiment.sh

Does everything to prepare data and run experiment.
Takes no command line arguments.

B) runPredict.sh <sampleId>

This can only be run after runExperiment.sh because it requires the csv output containing STRs and SNPs for each sample as input.

Takes a single command line argument representing a sample id in your data directory.

It will train a model using as many STRs as present in prediction sample and print a result to an output file in the data directory called "prediction". This contains the set of STRs used to train the model, model train time, and the ranked order of predicted classes.

Description of Scripts

Four scripts are ran in order as part of the experiment:

1) parseResults.py
    for each sample id, creates a directory containing "str", "pos", "neg" csv files containing STR allele-value pairs, SNP positives and negatives in a single-line CSV format

2) cladeFinder.py
    runs cladeFinder algorithm on each sample, taking positives and negatives as inputs to determine terminal subclade
    writes out terminal subclade to cladeFinderOutputFile

3) createCSVinputForRF.py
    creates CSV file containing for each sample, the STRs, SNP-calculated terminal subclade, downstream branches not confirmed negative

4) RF_predictor.py
    increases resolution of training data by inferring a downstream subclade for low SNP resolution samples that match another sample of a SNP-confirmed downstream branch within a STR GD cutoff and which the sample is not confirmed negative for
    creates random forest model with a holdout to use for testing
    reports accuracy per class along with ids which were incorrectly predicted over course of 25 runs

    these variables may be overridden by setting them in config.txt:
          
        iterations=25 #number of model runs to execute for calculating average accuracy
        refine_hg_based_on_strs_gd_cutoff=12 #Genetic distance cutoff used to increase SNP resolution of samples with closest matches
        model=ab #any combination of a,b,c,d - each corresponds to a set of STRs

        a = ["DYS391","DYS389I","DYS437","DYS439","DYS389II","DYS438","DYS426","DYS393","YCAII","DYS390","DYS385","Y-GATA-H4","DYS388","DYS447","DYS19","DYS392"]
        b = ["DYS458","DYS455","DYS454","DYS464","DYS448","DYS449","DYS456","DYS576","CDY","DYS460","DYS459","DYS570","DYS607","DYS442"]
        c = ["DYS728","DYS723","DYS711","DYR76","DYR33","DYS727","DYR157","DYS713","DYS531","DYS578","DYF395","DYS590","DYS537","DYS641","DYS472","DYF406S1","DYS511","DYS557","DYS490","DYS446","DYS481","DYS413","DYS534","DYS450","DYS425","DYS594","DYS444","DYS520","DYS436","DYS565","DYS572","DYS617","DYS568","DYS487","DYS640","DYS492"]
        d = ["DYR112","DYS518","DYS614","DYS626","DYS644","DYS684","DYS710","DYS485","DYS632","DYS495","DYS540","DYS714","DYS716","DYS717","DYS505","DYS556","DYS549","DYS589","DYS522","DYS494","DYS533","DYS636","DYS575","DYS638","DYS462","DYS452","DYS445","Y-GATA-A10","DYS463","DYS441","Y-GGAAT-1B07","DYS525","DYS712","DYS593","DYS650","DYS532","DYS715","DYS504","DYS513","DYS561","DYS552","DYS726","DYS635","DYS587","DYS643","DYS497","DYS510","DYS434","DYS461","DYS435"]
