<!DOCTYPE html>
<html>
<head>
        <link rel="stylesheet" type="text/css" href="https://www.yseq.net/ext/jquery/ui/redmond/jquery-ui-1.8.22.css" />
        <link rel="stylesheet" type="text/css" href="https://www.yseq.net/ext/jquery/fancybox/jquery.fancybox-1.3.4.css" />
        <link rel="stylesheet" type="text/css" href="https://www.yseq.net/ext/960gs/960_24_col.css" />
        <link rel="stylesheet" type="text/css" href="https://www.yseq.net/stylesheet.css" />

</head>
<body>
        <div id="bodyWrapper" class="container_24">


                <div id="header" class="grid_24">
                  <div id="storeLogo"><a href="https://www.yseq.net/index.php"><img src="https://www.yseq.net/images/store_logo.png" alt="YSEQ DNA Shop" title="YSEQ DNA Shop" width="552" height="175" /></a></div>
                

                
                <div class="grid_24 ui-widget infoBoxContainer">
                  <div class="ui-widget-header infoBoxHeading">&nbsp;&nbsp;<a href="http://www.yseq.net" class="headerNavigation">Top</a> &raquo; <a href="https://www.yseq.net/index.php" class="headerNavigation">Catalog</a></div>
                </div>


<div style="display:inline-block;text-align: left;">
        <h1>YSEQ Haplogroup Predictor (beta version)</h1><br>
    <div style="display:inline-block;padding-top: 50px;padding-right: 30px;padding-bottom: 50px;padding-left: 30px;"><img src="RandomForest.jpg" alt="Random Forest" title="Random Forest" width="300" />
</div><div style="display:inline-block; width:50%;padding-top: 50px;padding-right: 30px;padding-bottom: 50px;padding-left: 80px;">
Developed by Hunter Provyn with input and support from Thomas Krahn (2019).<br>
The haplogroup prediction is backed by a Random Forest model implemented in python using <a href="https://scikit-learn.org/stable/">sklearn</a>.<br><br>

This YSEQ haplogroup predictor is open source software and can be cloned from GitHub: <a href="https://github.com/hprovyn/str-to-haplogroup-predictor">https://github.com/hprovyn/str-to-haplogroup-predictor</a><br><br>

Please always give a link to <a href="http://predict.yseq.net">this</a> original website as a reference.<br><br><br>

</div></div>
<br>
<div style="padding-right: 30px;padding-bottom: 50px;padding-left: 30px;text-align: left">
<h1>Enter STR Alleles for Y-Haplogroup Prediction</h1><br>
Enter STRs in one of two formats, then press ENTER:<br><br>
<br>
<form action="<?php echo htmlentities($_SERVER['PHP_SELF']); ?>" method="POST">
<input type="radio" id="standard" name="format" value="standard" checked> Standard Format: $STR1=$ALLELE1,$STR2=$ALLELE2, ... <b>OR</b> $STR1 $ALLELE1 $STR2 $ALLELE2 ... <br>
<input type="radio" id="ftdna" name="format" value="ftdna"> FTDNA Tab Separated Format (copy results table row): $ALLELE1 TAB $ALLELE2 TAB $ALLELE3... in default FTDNA order containing 12, 25, 37, 67 or 111 STR Alleles <br>

<?php if(isset($_POST['input'])) { ?>
        <input name=input type="text" value="<?php echo $_POST['input']; ?>" maxlength="1500" size="135"></input>
        <script>
            var format = "<?php echo $_POST['format'] ?>"
            if (format == "ftdna") {
                document.getElementById("ftdna").checked = true;
            }
        </script>
        <?php
} else { ?>
        <input name=input type="text" maxlength="1500" size="135"></input><?php
} ?>
</form>

<br><br>
<b>Experiment Information</b><br>
<iframe src="modelInfo.php" width="880" height="280"></iframe><br>
<br>
<i>The YSEQ Haplogroup Predictor was influenced by the ideas of <a href="http://www.hprg.com/hapest5/hapest5b/hapest5.htm">Whit Athey's Haplogroup Predictor</a> (2006) and the <a href="http://www.nevgen.org">NevGen Predictor</a> from Milos Cetkovic Gentula and Aco Nevski (2014).<br><br>
While both mentioned haplogroup predictors use a Bayesian-Allele-Frequency approach, this YSEQ predictor uses machine learning with the random forest technique. <br>
Machine learning is in its infancy, so this predictor is unlikely to give you better or more precise results than the Bayesian predictor types, but at least you can consider it as a <b>second opinion</b> 
with an independent method.<br> 
Note that the YSEQ predictor is based on results of YSEQ customers, but it doesn't reveal the underlying STR profiles and original sample donors (they are not even stored on this web server at all). 
A computer based random number generator is the origin for creating random decision trees 
which are then just tested with a real life truth set. The teaching process simply selects the best decision trees and uses them for the prediction process.<br>
Please consider this beta version as an experiment which is largely untested and which needs a lot of improvements for reliable haplogroup predictions. We hope that when the number of samples increases, more and more 
outlier cases will be covered and considered for the prediction. We hope that this tool will become useful for the genetic genealogy community. The YSEQ Haplogroup Predictor comes with <b>no warranty</b>, explicit or implied.</i><br>
<br>


<?php if(isset($_POST['input'])) { 
        $corrected = "";
        ?>
        <?php $input = trim($_POST["input"]);
        $parsed=$input;
        
        if ($_POST["format"] == "standard") {
            $alleles = explode(",", $input);
            if (strpos($input, '=') === false) {
                $pieces = explode(" ", $input);
                $alleles = array();
                for ($x = 0; $x < count($pieces) / 2; $x++) {
                array_push($alleles, $pieces[$x * 2] . "=" . $pieces[$x * 2 + 1]);
            }
            $parsed = join(",", $alleles);
            }
        } else {
            $ftdnaHeaderOrder = array("DYS393","DYS390","DYS19","DYS391","DYS385","DYS426","DYS388","DYS439","DYS389I","DYS392","DYS389II","DYS458","DYS459","DYS455","DYS454","DYS447","DYS437","DYS448","DYS449","DYS464","DYS460","Y-GATA-H4","YCAII","DYS456","DYS607","DYS576","DYS570","CDY","DYS442","DYS438","DYS531","DYS578","DYF395","DYS590","DYS537","DYS641","DYS472","DYF406S1","DYS511","DYS425","DYS413","DYS557","DYS594","DYS436","DYS490","DYS534","DYS450","DYS444","DYS481","DYS520","DYS446","DYS617","DYS568","DYS487","DYS572","DYS640","DYS492","DYS565","DYS710","DYS485","DYS632","DYS495","DYS540","DYS714","DYS716","DYS717","DYS505","DYS556","DYS549","DYS589","DYS522","DYS494","DYS533","DYS636","DYS575","DYS638","DYS462","DYS452","DYS445","Y-GATA-A10","DYS463","DYS441","Y-GGAAT-1B07","DYS525","DYS712","DYS593","DYS650","DYS532","DYS715","DYS504","DYS513","DYS561","DYS552","DYS726","DYS635","DYS587","DYS643","DYS497","DYS510","DYS434","DYS461","DYS435");
        
            $strs = preg_split("/[\t]/", $input);
            $alleles = array();
            for ($x = 0; $x < count($strs); $x++) {
                if ($ftdnaHeaderOrder[$x] == "Y-GATA-H4") {
                    $wrongYGataH4 = intval($strs[$x]);
                    $corrected = $wrongYGataH4 + 1;
                    array_push($alleles, $ftdnaHeaderOrder[$x] . "=" . $corrected);
                    $corrected = "<br>NOTE: corrected FTDNA's Y-GATA-H4 value to proper NIST format, from " . $wrongYGataH4 . " to " . $corrected . "<br>";
                } else {
                    array_push($alleles, $ftdnaHeaderOrder[$x] . "=" . trim($strs[$x]));
                }
            }
            $parsed = join(",", $alleles);
        }
        $parsedTable = "<table border=&quot;1&quot;>";
        $strRow = "<tr><td>STR</td>";
        $alleleRow = "<tr><td>Allele</td>";
        for ($x = 0; $x < count($alleles); $x++) {
        $alleleExploded = explode("=", $alleles[$x]);
        $strRow = $strRow . "<td>" . $alleleExploded[0] . "</td>";
        $alleleRow = $alleleRow . "<td>" . $alleleExploded[1] . "</td>";
        }
        $parsedTable = $parsedTable . $strRow . "</tr>" . $alleleRow . "</tr></table>";
        $strInputHTML = $parsedTable . $corrected;
        echo '<br><b>STR Input</b><div style="outline: 1px solid black" width="880" height="150"><div id="strinput" style="padding:10px;overflow-x:scroll;max-width:860px"></div></div><br><br>';
        
        $message = exec('/var/lib/str-to-haplogroup-predictor/runPredict.sh ' . $parsed);
        $outputsplit = explode("Model Information", $message);
        $predsplit = str_replace("&", "&amp;", $outputsplit[0]);
        $predsplit = str_replace("\"", "&quot;", $predsplit);
        echo '<b>Prediction</b><div style="outline: 1px solid black" width="880" height="250"><div id="pred" style="padding:10px;overflow-y: scroll;max-height:230px" id="pred" width="860" height="230"></div></div><br><br>';
        $modelsplit = str_replace("&", "&amp;", $outputsplit[1]);
        $modelsplit = str_replace("\"", "&quot;", $modelsplit);
        echo '<b>Model Information</b><div style="outline: 1px solid black" width="880" height="400"><div id="model" style="padding:10px;overflow-y: scroll;max-height:380px"></div></div>';
        $updateSrc = true;
} ?>

</div>
<script>
function removeHTMLCodes(html) {
    var a = html.replace(/&amp;#9608;/g, "█")
    a = a.replace(/&quot;/g, "\"")
    a = a.replace(/&amp;nbsp;/g, " ")
    return a
}

var content = " <?php echo $predsplit ?> "
content = removeHTMLCodes(content)
document.getElementById("pred").innerHTML = content;
var content = " <?php echo $strInputHTML ?> "
content = removeHTMLCodes(content)
document.getElementById("strinput").innerHTML = content;
var content = " <?php echo $modelsplit ?> "
content = removeHTMLCodes(content)
document.getElementById("model").innerHTML = content;
</script>

</body>
</html>
