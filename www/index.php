
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

 <p><strong>The package is useful if</strong></p>
 <ul>
         <li align="justify">
         <p align="justify">one wants to perform model selection or model averaging on multiply imputed data and the analysis model of interest is either the linear model, the logistic model, the Poisson model, or the Cox proportional hazards model, possibly with a random intercept.</p>
         </li>
         <li align="justify">
         <p align="justify">one wants to obtain bootstrap confidence intervals for model selection or model averaging estimators (with or without missing data/imputation) -- to address model selection uncertainty and to discover relationships of a small effect size.</p>
         </li>
         <li align="justify">
         <p align="justify">one wants to compare different model selection and averaging techniques, easily with the same syntax.</p>
         </li>
 </ul>
 
  <p><strong>The package is of limited use under the following circumstances:</strong></p>
 <ul>
         <li align="justify">
         <p align="justify">if one is interested in model selection or averaging for models other than those listed above, for example parametric survival models, additive models, time-series analysis, and many others. </p>
         </li>
         <li align="justify">
         <p align="justify">if one decides for a specific model selection or averaging technique not provided by the package, look at the manual for more details.</p>
         </li>
          <li align="justify">
         <p align="justify">if the model selection/averaging problem is computationally too intensive, see Section 6.1 from the manual for more details.</p>
         </li> </ul>
 
<p><strong>Manual</strong></p>
The package manual can be found  <a href="MAMI_manual.pdf"><strong>here</strong></a>.
 
 <p><strong>References</strong></p>
<p>Schomaker, M., Heumann, C. (2014) Model Selection and Model Averaging after Multiple Imputation,
Computational Statistics & Data Analysis, 71:758-770</p>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

<p> The package can be downloaded <a href="https://r-forge.r-project.org/R/?group_id=2152"><strong>here</strong></a> (do not forget to also install the dependencies).</p> <p>Or simply type <tt> install.packages("MAMI", repos=c("http://R-Forge.R-project.org","http://cran.at.r-project.org"), dependencies=TRUE)</tt> in <i>R</i>. </p>

</body>
</html>
