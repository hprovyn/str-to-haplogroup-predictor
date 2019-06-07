#PARAMETERS
refine_hg_based_on_strs_gd_cutoff=12
model=ab
iterations=25

CFG_FILE=config.txt
CFG_CONTENT=$(cat $CFG_FILE | sed -r '/[^=]+=[^=]+/!d' | sed -r 's/\s+=\s/=/g')
eval "$CFG_CONTENT"

PATH=$PATH:$pythonPath

#DO NOT TOUCH
parseResults_py="$pythonScriptsDir\parseResults.py"
createCSVinputForRF_py="$pythonScriptsDir\createCSVinputForRF.py"
rfPredictor_py="$pythonScriptsDir\RF_predictor.py"
cladeFinder_py="$pythonScriptsDir\cladeFinder.py"

dataDir="$experimentDir\data"
cladeFinderOutputFile="$experimentDir\sampleCladesFound"
csvOutputForRFfile="$experimentDir\csv_for_rf.csv"
modelOutputFile="$experimentDir\modelOutput.csv"

#python "$parseResults_py" "$resultsFile" "$dataDir"
#python "$cladeFinder_py"  "$treeFile" "$dataDir" "$cladeFinderOutputFile"
#python "$createCSVinputForRF_py" "$treeFile" "$dataDir" "$cladeFinderOutputFile" "$csvOutputForRFfile"
python "$rfPredictor_py"  "$csvOutputForRFfile"  "$modelOutputFile" "$model" "$refine_hg_based_on_strs_gd_cutoff" "$iterations"
