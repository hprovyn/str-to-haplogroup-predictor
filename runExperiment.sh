#SET PARAMS FROM CONFIG
CFG_FILE=config.txt
CFG_CONTENT=$(cat $CFG_FILE | sed -r '/[^=]+=[^=]+/!d' | sed -r 's/\s+=\s/=/g')
eval "$CFG_CONTENT"

PATH=$PATH:$pythonPath

parseResults_py="$pythonScriptsDir\parseResults.py"
createCSVinputForRF_py="$pythonScriptsDir\createCSVinputForRF.py"
rfExperiment_py="$pythonScriptsDir\RF_experiment.py"
cladeFinder_py="$pythonScriptsDir\cladeFinder.py"

dataDir="$experimentDir\data"
cladeFinderOutputFile="$experimentDir\sampleCladesFound"
csvOutputForRFfile="$experimentDir\csv_for_rf.csv"
modelOutputFile="$experimentDir\modelOutput.csv"
panelHierarchyFile="$experimentDir\panelHierarchy.csv"
predictionPolicyFileStem="$experimentDir\predictionPolicy"
modelPickleFileStem="$experimentDir\model"
experimentMapFileStem="$experimentDir\experimentMap"

python "$parseResults_py" "$resultsFile" "$dataDir"
python "$cladeFinder_py"  "$treeFile" "$dataDir" "$cladeFinderOutputFile"
python "$createCSVinputForRF_py" "$treeFile" "$dataDir" "$cladeFinderOutputFile" "$csvOutputForRFfile" "$panelHierarchyFile"
python "$rfExperiment_py" "$csvOutputForRFfile" "$modelOutputFile" "$panelHierarchyFile" "$predictionPolicyFileStem" "$modelPickleFileStem" "$utilityUnderSpecificityError" "$utilityOverSpecificityError" "$utilityCompletelyWrongError" "$utilityCorrect" "$experimentMapFileStem"