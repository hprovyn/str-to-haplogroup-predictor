<html>
<body>

STRs input: <?php echo $_POST["input"]; ?>
<br><br>
<?php $input = trim($_POST["input"]);
$parsed=$input;
if (strpos($input, '=') === false) {
    
$pieces = explode(" ", $input);
$alleles = array();
for ($x = 0; $x < count($pieces) / 2; $x++) {
array_push($alleles, $pieces[$x * 2] . "=" . $pieces[$x * 2 + 1]);
}
$parsed = join(",", $alleles);
}
echo "STRs parsed: " . $parsed . "<br><br>";
$message = exec('/var/lib/str-to-haplogroup-predictor/runPredict.sh ' . $parsed);
echo "Your predicted haplogroup is " . $message; ?>

</body>
</html>
