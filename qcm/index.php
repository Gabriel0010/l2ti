<?php
session_start();
require_once './classes/Pages.class.php';
require_once './classes/Page.class.php';
require_once './classes/Question.class.php';
require_once './classes/Qcm.class.php';
require_once './classes/Qcms.class.php';


# Modif FORTIER
#require_once './controllers/Controlpages.php';
require_once './controllers/ControlPages.php';
?>

<!DOCTYPE html>
<html lang="fr">
    <head>
        <meta charset="utf-8">
        <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <meta name="description" content="">
        <meta name="author" content="">

        <!-- <link rel="icon" href="../../favicon.ico"> -->

        <title>QCM Traitement du signal - <?php echo $page->getTitle(); ?></title>
<?php include_once('partials/_stylesheets.php'); ?>
    </head>

    <body>

    <div class="container">

        <?php include_once('partials/_menu.php'); ?>
                    <!-- Example row of columns -->

        <?php include_once($page->getContent()); ?>

                <!-- Site footer -->
        <?php include_once('partials/_footer.php') ?>

    </div> <!-- /container -->


    <!-- Bootstrap core JavaScript
    ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->
<?php include_once('partials/_javascripts.php') ?>

</body>

</html>
