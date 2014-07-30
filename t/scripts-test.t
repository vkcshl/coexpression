use strict;
use warnings;
use Bio::KBase::workspaceService::Client;
use Bio::KBase::AuthToken;
use Bio::KBase::CoExpression::Client;
use JSON::XS;
use Test::More;
use Data::Dumper;
use File::Temp qw(tempfile);
use Getopt::Long::Descriptive;
use Text::Table;
use Test::Cmd;


#by Fei, July 29, 2014.
#use kbase-login as kbasetest before running this script


my $server="http://140.221.85.182:7063";
my $ws_id='coexpression_randomtest_';
my $series_id='test1';
my $series_fid='test1_filted_by_scripttest';
my $series_netid='test1_net1';
my $ws_url='https://kbase.us/services/ws';

my $tmp;
my $bin="scripts";


my $tes = Test::Cmd->new(prog => "$bin/coex-filter-genes.py", workdir => '', interpreter => '/usr/bin/python');
ok($tes, "creating Test::Cmd object for coex-filter-genes.py");
$tes->run(args => '-h');
ok($? == 0,"Running coex-filter-genes: produing help file");

$tes->run(args => "--ws_url=$ws_url  --ws_id=$ws_id   --in_id=$series_id --out_id=$series_fid --filter_method=anova --p_value=0.01");
ok($? == 0,"Running coex-filter-genes, anova, 0.01 ");
$tmp=$tes->stdout;
ok($tmp =~ /359/, "identified 359 differential genes using anova and p=0.01");

$tes->run(args => "--ws_url=$ws_url  --ws_id=$ws_id   --in_id=$series_id --out_id=$series_fid --filter_method=anova --p_value=0.05");
ok($? == 0,"Running coex-filter-genes, anova, 0.05 ");
 $tmp=$tes->stdout;
ok($tmp =~ /1131/, "identified 1131 differential genes using anova and p=0.05");


$tes = Test::Cmd->new(prog => "$bin/coex-net-clust.py", workdir => '', interpreter => '/usr/bin/python');
ok($tes, "creating Test::Cmd object for coex-filter-genes.py");

$tes->run(args => '-h');
ok($? == 0,"Running coex-filter-genes: produing help file");

$tes->run(args => "--ws_url=$ws_url  --ws_id=$ws_id   --in_id=$series_fid --out_id=$series_netid -c 0.8 -n simple -m hclust  -k 30");
ok($? == 0,"Running coex-filter-genes, 0.8, 30 ");
$tmp=$tes->stdout;
ok($tmp=~/modulenames/, "running coex-net-clust: produing network clusters. hclust");


$tes->run(args => "--ws_url=$ws_url  --ws_id=$ws_id   --in_id=$series_fid --out_id=$series_netid -c 0.8 -n WGCNA  -m WGCNA -k 30");
ok($? == 0,"Running coex-filter-genes, 0.8, 30 ");
$tmp=$tes->stdout;
ok($tmp=~/modulenames/, "running coex-net-clust: produing network clusters. wgcna");








