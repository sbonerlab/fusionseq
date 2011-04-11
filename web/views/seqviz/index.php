		<h1>SeqViz: Identification of all reads which align to one region</h1>
		<br /><br />
		<form action="<? echo WEB_BASE_URL; ?>/seqviz/index.php" method="get">
		<label>Data prefix:</label>
  		<input type="text" name="prefix">
  		<br><br><br>
  		<label>Locus:</label>
		<input type="text" name="location">
		Example: chr21:2449394-3499493
 		<br><br><br>
		<label>Limit number of reads:</label>
		<input type="text" name="readlim">
		<br><br><br>
 		<label>Minimum insert size:</label>
		<input type="text" name="minspan">
		<br><br><br>
		<input type="submit" value="Submit">
		<input type="reset" value="Reset">
		</form>
