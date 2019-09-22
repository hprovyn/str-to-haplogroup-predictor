#SET PARAMS FROM CONFIG
CFG_FILE=/var/lib/str-to-haplogroup-predictor/config.txt
CFG_CONTENT=$(cat $CFG_FILE | sed -r '/[^=]+=[^=]+/!d' | sed -r 's/\s+=\s/=/g')
eval "$CFG_CONTENT"

PATH=$PATH:$pythonPath:$htslibPath

getModelMetadata_py="$pythonScriptsDir${pathSeparator}getModelMetadata.py"
generalModelMetadataFile="$predictionDir${pathSeparator}model_general_metadata"

python3 "$getModelMetadata_py" "$generalModelMetadataFile"