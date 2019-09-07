STR TO HAPLOGROUP PREDICTOR

This project is financially supported by YSEQ and is the work of Hunter Provyn with some ideas contributed by Thomas Krahn.

WHAT IT DOES

Processes sample STR and SNP results along with a SNP phylogenetic tree to run an experiment that predicts SNP branch from a set of STRs absent any SNP info. 

REQUIREMENTS

. Python 3.6 #Python dependencies
  . packages numpy, pandas, time, sys, sklearn, random, json, os, operator, shutil, pickle

. CSV file containing STR and SNP results where each line contains a single ID, MARKER, ALLELE
. JSON file containing YFull tree

INSTALLATION

sudo apt-get install python3

sudo apt-get install python3-pandas
sudo apt-get install python3-sklearn
sudo apt-get install python3-numpy

pytabix is not part of python3 distribution libraries. So you'll need to install it using pip.

pip3 install --user pytabix  (will install only for --user but is ok because web-data won't need this package)

EXPERIMENT vs PREDICT MODES

The runExperiment.sh script runs an experiment that trains an optimized model that can be used to predict SNP Haplogroup from a set of input STRs.

Files are created in two different directories in order to prevent exposure of the sample data to the webserver, which only needs to be able to deserialize a model to run a prediction:

$predictionDir - This contains only those files necessary to run a prediction
$experimentDir - This contains tsv formatted samples, its gz and tabix indexed versions and other experiment meta data

Example Secure Deployment Strategy

Install on two machines:

A) Experiment machine - has private sample data, no web-server
B) Web Server machine

1. Run the experiment on the Experiment machine.
2. scp the $predictionDir onto the Web Server machine.
3. Set the $predictionDir in config.txt to point to this

Now you have a webserver that can predict STR haplogroup that cannot possibly expose sample data, because it doesn't have it.


EXPLANATION OF SHELL FILES:

A) runExperiment.sh

Does everything to prepare data and run experiment.
Takes no command line arguments.
For all STR group combinations, trains and serializes a model and optimizes and persists :a prediction policy that in some cases changes prediction to a less specific haplogroup panel according to utility weights (if any are changed from default of zero in config file) or to maximize percent correct prediction.

The serialized models and other files necessary to make a prediction are stored in the $predictionDir.

B) runPredict.sh <$STR1=$ALLELE1,$STR2=$ALLELE2...>

This can only be run after runExperiment.sh because it requires the models that were trained during this process.

Takes a single command line argument representing the STR allele and value pairs you wish to predict the SNP haplogroup of.

Loads model from pickle file trained to as many STRs as present in the query, applies prediction policy computed in experiment step, and returns the ranked order of predicted classes in an html table.

The predict.php script calls runPredict.sh in order to display the predicted haplogroup in the browser. 

Description of Scripts

Four scripts are ran in order as part of the experiment:

1) parseResults.py
    for each sample id, creates a directory containing "str", "pos", "neg" csv files containing STR allele-value pairs, SNP positives and negatives in a single-line CSV format

2) cladeFinder.py
    runs cladeFinder algorithm on each sample, taking positives and negatives as inputs to determine terminal subclade
    writes out terminal subclade to cladeFinderOutputFile

3) createCSVinputForRF.py
    creates CSV file containing for each sample, the STRs, SNP-calculated terminal subclade, downstream branches not confirmed negative

4) RF_experiment.py
    creates random forest model with a holdout to use for testing
    reports accuracy per class along with ids which were incorrectly predicted over course of 25 runs
    these variables may be overridden by setting them in config.txt:
          
        iterations=25 #number of model runs to execute for calculating average accuracy
        cladeFinderMaxThreads= #number of parallel threads allowed to run clade finder against all samples
	createCSVinputMaxThreads= #number of parallel threads allowed to create CSV input file
	percentMissingSTRThreshold= #percentage (in int form) of STRs allowed to be missing in a sample for it to be considered valid for a model trained on a specific set of [a,b,c,d] STR sets. Missing values are set to zero.

	model=ab #any combination of a,b,c,d - each corresponds to a set of STRs
        a = ["DYS391","DYS389I","DYS437","DYS439","DYS389II","DYS438","DYS426","DYS393","YCAII","DYS390","DYS385","Y-GATA-H4","DYS388","DYS447","DYS19","DYS392"]
        b = ["DYS458","DYS455","DYS454","DYS464","DYS448","DYS449","DYS456","DYS576","CDY","DYS460","DYS459","DYS570","DYS607","DYS442"]
        c = ["DYS728","DYS723","DYS711","DYR76","DYR33","DYS727","DYR157","DYS713","DYS531","DYS578","DYF395","DYS590","DYS537","DYS641","DYS472","DYF406S1","DYS511","DYS557","DYS490","DYS446","DYS481","DYS413","DYS534","DYS450","DYS425","DYS594","DYS444","DYS520","DYS436","DYS565","DYS572","DYS617","DYS568","DYS487","DYS640","DYS492"]
        d = ["DYR112","DYS518","DYS614","DYS626","DYS644","DYS684","DYS710","DYS485","DYS632","DYS495","DYS540","DYS714","DYS716","DYS717","DYS505","DYS556","DYS549","DYS589","DYS522","DYS494","DYS533","DYS636","DYS575","DYS638","DYS462","DYS452","DYS445","Y-GATA-A10","DYS463","DYS441","Y-GGAAT-1B07","DYS525","DYS712","DYS593","DYS650","DYS532","DYS715","DYS504","DYS513","DYS561","DYS552","DYS726","DYS635","DYS587","DYS643","DYS497","DYS510","DYS434","DYS461","DYS435"]
