#SET PARAMS FROM CONFIG
#THIS MUST BE ABSOLUTE PATH TO CONFIG FILE
CFG_FILE=/var/lib/str-to-haplogroup-predictor/config.txt
CFG_CONTENT=$(cat $CFG_FILE | sed -r '/[^=]+=[^=]+/!d' | sed -r 's/\s+=\s/=/g')
eval "$CFG_CONTENT"

PATH=$PATH:$pythonPath

rfPredict_py="$pythonScriptsDir${pathSeparator}RF_predict.py"

dataDir="$experimentDir${pathSeparator}data"
panelHierarchyFile="$experimentDir${pathSeparator}panelHierarchy.csv"
predictionPolicyFileStem="$experimentDir${pathSeparator}predictionPolicy"
modelPickleFileStem="$experimentDir${pathSeparator}model"

python "$rfPredict_py" "$dataDir" "$panelHierarchyFile" "$predictionPolicyFileStem" "$modelPickleFileStem" "$percentMissingSTRThreshold" $1
