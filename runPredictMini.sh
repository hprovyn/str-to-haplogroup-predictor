#SET PARAMS FROM CONFIG
#THIS MUST BE ABSOLUTE PATH TO CONFIG FILE
CFG_FILE=/var/lib/str-to-haplogroup-predictor/config.txt
CFG_CONTENT=$(cat $CFG_FILE | sed -r '/[^=]+=[^=]+/!d' | sed -r 's/\s+=\s/=/g')
eval "$CFG_CONTENT"

PATH=$PATH:$pythonPath

rfPredict_py="$pythonScriptsDir${pathSeparator}RF_predict.py"

panelHierarchyFile="$predictionDir${pathSeparator}panelHierarchy.csv"
predictionPolicyFileStem="$predictionDir${pathSeparator}predictionPolicy"
modelPickleFileStem="$predictionDir${pathSeparator}model"

outputVersion="mini"

python3 "$rfPredict_py" "$panelHierarchyFile" "$predictionPolicyFileStem" "$modelPickleFileStem" "$percentMissingSTRThreshold" "$haplogroupClassConfigPath" "$outputVersion" $1
