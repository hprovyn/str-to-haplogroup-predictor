<html>
<body>

STRs input: <?php echo $_POST["input"]; ?>
<br><br>
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
            echo "<br>NOTE: corrected FTDNA's Y-GATA-H4 value to proper NIST format, from " . $wrongYGataH4 . " to " . $corrected . "<br>";
        } else {
            array_push($alleles, $ftdnaHeaderOrder[$x] . "=" . trim($strs[$x]));
        }
    }
    $parsed = join(",", $alleles);
}
$parsedTable = "<table border=\"1\">";
$strRow = "<tr><td>STR</td>";
$alleleRow = "<tr><td>Allele</td>";
for ($x = 0; $x < count($alleles); $x++) {
$alleleExploded = explode("=", $alleles[$x]);
$strRow = $strRow . "<td>" . $alleleExploded[0] . "</td>";
$alleleRow = $alleleRow . "<td>" . $alleleExploded[1] . "</td>";
}
$parsedTable = $parsedTable . $strRow . "</tr>" . $alleleRow . "</tr></table>";
echo $parsedTable . "<br><br>";

$message = exec('/var/lib/str-to-haplogroup-predictor/runPredict.sh ' . $parsed);
echo $message; ?>

</body>
</html>
