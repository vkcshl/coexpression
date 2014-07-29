use strict;
use Bio::KBase::CoExpression::Client;
use Test::More tests => 11;
use Data::Dumper;
use Test::Cmd;
use JSON;


my $ws_id="kbasetest:home";
#my $cc = Bio::KBase::CoExpression::Client->new("http://localhost:7063");

print "This test requires the kb|series.436 is stored in kbasetest:home(workspace)\n";
print "It also requires read and write ability in the kbasetest:home(workspace)\n";
print "use Shinjae's server for client testing\n";
my $cc = Bio::KBase::CoExpression::Client->new("http://140.221.85.182:7063");
ok( defined $cc, "Check if the server is working" );

my $job_id = $cc->filter_genes({ws_id => $ws_id, inobj_id => 'kb|series.436', outobj_id => 'kb|series.436.ftd', p_value => "0.05", method => "anova"});
ok(ref($job_id) eq "ARRAY","filter_genes returns an array");
ok(@{$job_id} eq 2, "returns two job ids for filter_genes");


$job_id = $cc->const_coex_net_clust({ws_id => $ws_id, inobj_id => 'kb|series.436.ftd', outobj_id => 'net_by_kb|series.436', cut_off => "0.75", net_method => "simple", clust_method => "hclust", num_modules => "10"});
ok(ref($job_id) eq "ARRAY","const_coex_net_clust returns an array");
ok(@{$job_id} eq 2, "returns two job ids for const_coex_net_clust");

print "download the newly-generated network object from workspace and check the number of clusters in it\n";
`kbws-get 'net_by_kb|series.436' >tmp_nto`;
open II,"tmp_nto" or die "no NTO file";
my @all=<II>;
close II;
my $data=decode_json($all[0]);
print "get the name of cluster in the NTO\n";

my $no_cluster=0;
my $i=0;
for($i=0;$i<@{$data->{'nodes'}};$i++){
	next if $data->{'nodes'}->[$i]->{'entity_id'} !~/cluster/;
	$no_cluster++;
}
print "No. of cluster expected:10\n";
print "No. of cluster generated:$no_cluster\n";
ok($no_cluster==10, "Generating 10 clusters in the networks\n");


`rm tmp_nto`;


print "testing different parameters\n";
$job_id = $cc->filter_genes({ws_id => $ws_id, inobj_id => 'kb|series.436', outobj_id => 'kb|series.436.ftd',num_genes=>39 , method => "lor"});
ok(ref($job_id) eq "ARRAY","filter_genes returns an array");
ok(@{$job_id} eq 2, "returns two job ids for filter_genes");

$job_id = $cc->const_coex_net_clust({ws_id => $ws_id, inobj_id => 'kb|series.436.ftd', outobj_id => 'net_by_kb|series.436', cut_off => "0.75", net_method => "w", clust_method => "w", num_modules => "50"});
ok(ref($job_id) eq "ARRAY","const_coex_net_clust returns an array");
ok(@{$job_id} eq 2, "returns two job ids for const_coex_net_clust");



$i=int rand(10000);
my $newname="kb|series.436.client-testing".$i;
$job_id = $cc->filter_genes({ws_id => $ws_id, inobj_id => 'kb|series.436', outobj_id => $newname, num_genes=>379 , method => "lor"});

print "download the newly-generated network object from workspace and check if it is existed\n";
`kbws-get $newname >tmp_to`;
open II,"tmp_to" or die "no TO file";
@all=<II>;
close II;
ok($all[0],"new object has been generated in workspace");
#print Dumper($data);



















