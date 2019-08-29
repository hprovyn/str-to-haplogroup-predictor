#SET PARAMS FROM CONFIG
CFG_FILE=config.txt
CFG_CONTENT=$(cat $CFG_FILE | sed -r '/[^=]+=[^=]+/!d' | sed -r 's/\s+=\s/=/g')
eval "$CFG_CONTENT"

PATH=$PATH:$pythonPath:$htslibPath

parseResults_py="$pythonScriptsDir${pathSeparator}parseResults.py"
createCSVinputForRF_py="$pythonScriptsDir${pathSeparator}createCSVinputForRF.py"
rfExperiment_py="$pythonScriptsDir${pathSeparator}RF_experiment.py"
cladeFinder_py="$pythonScriptsDir${pathSeparator}cladeFinder.py"

samplestsv="$experimentDir${pathSeparator}samplestsv"
cladeFinderOutputFile="$experimentDir${pathSeparator}sampleCladesFound"
csvOutputForRFfile="$experimentDir${pathSeparator}csv_for_rf.csv"
modelOutputFile="$predictionDir${pathSeparator}modelOutput.csv"
panelHierarchyFile="$predictionDir${pathSeparator}panelHierarchy.csv"
predictionPolicyFileStem="$predictionDir${pathSeparator}predictionPolicy"
modelPickleFileStem="$predictionDir${pathSeparator}model"
experimentMapFileStem="$experimentDir${pathSeparator}experimentMap"

rm $experimentDir${pathSeparator}*
rm $predictionDir${pathSeparator}*
python3 "$parseResults_py" "$resultsFile" "$samplestsv"
sort "$samplestsv" "-k1V" "-k2n" "-k3n" | "bgzip" > "$samplestsv.bgz"
tabix "-s" "1" "-b" "2" "-e" "3" "$samplestsv.bgz"
python3 "$cladeFinder_py"  "$treeFile" "$samplestsv.bgz" "$cladeFinderOutputFile" "$cladeFinderMaxThreads"
python3 "$createCSVinputForRF_py" "$treeFile" "$samplestsv.bgz" "$cladeFinderOutputFile" "$csvOutputForRFfile" "$panelHierarchyFile" "$createCSVinputMaxThreads"
python3 "$rfExperiment_py" "$csvOutputForRFfile" "$modelOutputFile" "$panelHierarchyFile" "$predictionPolicyFileStem" "$modelPickleFileStem" "$utilityUnderSpecificityError" "$utilityOverSpecificityError" "$utilityCompletelyWrongError" "$utilityCorrect" "$experimentMapFileStem" "$percentMissingSTRThreshold"
