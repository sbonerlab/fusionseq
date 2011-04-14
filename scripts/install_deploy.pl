#!/usr/bin/env perl

use strict;
use warnings;

use Confp;
use Cwd;

my $path = getcwd() . "/../default.fusionseqrc";
my $config = Confp->new($path);

die("Cannot open default.fusionseqrc") if !$config;

my $cgis = qw(geneFusions_cgi showDetails_cgi seqViz_cgi findFusionPartner_cgi);
my $cmd = "mv " . $cgis . " " . $config->get('WEB_DATA_DIR');

system($cmd);
