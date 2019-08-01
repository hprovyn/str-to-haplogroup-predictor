<html>
<body>

STRs input: <?php echo $_POST["input"]; ?>
<br><br>
<?php $input = trim($_POST["input"]);
$parsed=$input;
$alleles = explode(",", $input);
if (strpos($input, '=') === false) {
    
$pieces = explode(" ", $input);
$alleles = array();
for ($x = 0; $x < count($pieces) / 2; $x++) {
array_push($alleles, $pieces[$x * 2] . "=" . $pieces[$x * 2 + 1]);
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
