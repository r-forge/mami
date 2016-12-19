
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

 <p><strong>When is the package useful for me?</strong></p>
 <ul>
         <li align="justify">
         <p align="justify">If you want to perform model selection or model averaging on multiply imputed data and your analysis model
  of interest is either the linear model, the logistic model, the Poisson model, or the Cox proportional hazards
  model, possibly with a random intercept.</p>
         </li>
         <li align="justify">
         <p align="justify">If you want to obtain bootstrap confidence intervals for model selection or model averaging estimators.</p>
         </li>
 </ul>
 
  <p><strong>When is the package not useful for me?</strong></p>
 <ul>
         <li align="justify">
         <p align="justify">If you are interested in model selection or averaging for models other than those listed above
  e.g. parametric survival models, time series analysis, ... </p>
         </li>
         <li align="justify">
         <p align="justify">If you decide for a specific model selection or averaging technique not provided by the package,
  e.g. recent shrinkage methods, boosting, p-value based approaches, in depth Bayesian Model Averaging, etc. You may however still want to consult the reference below.</p>
         </li>
 </ul>
 
<p><strong>How well is the package developed until now?</strong></p>
The package is generally stable and can be used. However, some parts of the code are still not very efficient. Future releases may also 
offer parallelization. The choice of models is also still limited. Those interested in model selection/averaging after multiple imputation 
for settings which are not contained in the package (yet) are referred to the reference below.
 
 <p><strong>References</strong></p>
<p>Schomaker, M., Heumann, C. (2014) Model Selection and Model Averaging after Multiple Imputation,
Computational Statistics & Data Analysis, 71:758-770</p>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

<p> The package can be downloaded <a href="https://r-forge.r-project.org/R/?group_id=2152"><strong>here</strong></a>.</p> <p>Or simply type <tt> install.packages("MAMI", repos=c("http://R-Forge.R-project.org","http://cran.at.r-project.org"), dependencies=TRUE)</tt> in <i>R</i>. </p>

</body>
</html>
