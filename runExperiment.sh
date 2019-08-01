#SET PARAMS FROM CONFIG
CFG_FILE=config.txt
CFG_CONTENT=$(cat $CFG_FILE | sed -r '/[^=]+=[^=]+/!d' | sed -r 's/\s+=\s/=/g')
eval "$CFG_CONTENT"

PATH=$PATH:$pythonPath

parseResults_py="$pythonScriptsDir${pathSeparator}parseResults.py"
createCSVinputForRF_py="$pythonScriptsDir${pathSeparator}createCSVinputForRF.py"
rfExperiment_py="$pythonScriptsDir${pathSeparator}RF_experiment.py"
cladeFinder_py="$pythonScriptsDir${pathSeparator}cladeFinder.py"

dataDir="$experimentDir${pathSeparator}data"
cladeFinderOutputFile="$experimentDir${pathSeparator}sampleCladesFound"
csvOutputForRFfile="$experimentDir${pathSeparator}csv_for_rf.csv"
modelOutputFile="$experimentDir${pathSeparator}modelOutput.csv"
panelHierarchyFile="$experimentDir${pathSeparator}panelHierarchy.csv"
predictionPolicyFileStem="$experimentDir${pathSeparator}predictionPolicy"
modelPickleFileStem="$experimentDir${pathSeparator}model"
experimentMapFileStem="$experimentDir${pathSeparator}experimentMap"

python "$parseResults_py" "$resultsFile" "$dataDir"
python "$cladeFinder_py"  "$treeFile" "$dataDir" "$cladeFinderOutputFile"
python "$createCSVinputForRF_py" "$treeFile" "$dataDir" "$cladeFinderOutputFile" "$csvOutputForRFfile" "$panelHierarchyFile"
python "$rfExperiment_py" "$csvOutputForRFfile" "$modelOutputFile" "$panelHierarchyFile" "$predictionPolicyFileStem" "$modelPickleFileStem" "$utilityUnderSpecificityError" "$utilityOverSpecificityError" "$utilityCompletelyWrongError" "$utilityCorrect" "$experimentMapFileStem"
