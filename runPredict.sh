#SET PARAMS FROM CONFIG
CFG_FILE=config.txt
CFG_CONTENT=$(cat $CFG_FILE | sed -r '/[^=]+=[^=]+/!d' | sed -r 's/\s+=\s/=/g')
eval "$CFG_CONTENT"

PATH=$PATH:$pythonPath

rfPredict_py="$pythonScriptsDir\RF_predict.py"

dataDir="$experimentDir\data"
panelHierarchyFile="$experimentDir\panelHierarchy.csv"
predictionPolicyFileStem="$experimentDir\predictionPolicy"
modelPickleFileStem="$experimentDir\model"

python "$rfPredict_py" "$dataDir" "$panelHierarchyFile" "$predictionPolicyFileStem" "$modelPickleFileStem" "$percentMissingSTRThreshold" $1