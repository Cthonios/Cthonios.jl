filename="$1"
parameters="$2"
basename="${filename%.*}"

echo "Running python esp with file $filename"
echo $ESP_ROOT

CASREV=$CASREV \
ESP_ARCH=$ESP_ARCH \
ESP_ROOT=$ESP_ROOT \
LD_LIBRARY_PATH=$PYESPDEPLIBS:$LD_LIBRARY_PATH \
PATH=$ESP_ROOT/bin/:$PATH \
PYTHONPATH=$ESP_ROOT/lib/pyESP/:$PYTHONPATH \
UDUNITS2_XML_PATH=$UDUNITS2_XML_PATH \
python3 $ESP_JULIA_ROOT/esp/scripts/my_script.py $filename > esp.out
# python ~/dev/plato/platoengine/share/plato-cli geometry esp \
#     --input $filename \
#     --output-model "${basename}_opt.csm" \
#     --output-mesh "${basename}_temp.exo" \
#     --parameters $parameters \
#     --tesselation "${basename}_temp.eto" \
#     --workflow aflr2
