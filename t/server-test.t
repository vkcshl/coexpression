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
use Bio::KBase::workspace::ScriptHelpers qw(get_ws_client parseWorkspaceInfo);

#this is a client testing script for coexpression service

my $para;
my $info;
my $ws_name;


my $ws_client = get_ws_client();
#print Dumper($ws_client);
#die;

=head
#Creating workspace for tests
my $num=int rand(100000);
#my $ws_name="coexpression_randomtest_".$num;
my $ws_name="coexpression_randomtest_";

 $para={
                workspace => $ws_name,
                globalread => "n"
};
my $info;
ok($info = $ws_client->create_workspace($para),"creating a new workspace for testing");
#print Dumper($info);
ok( $info->[0]=~/\d/ && $info->[0] !~/[a-z]/i, "a new workspace has been created");
print "new workspace id: $info->[0]\n";



#"uploading testing expression series object into $ws_name\n";
open MM,"expression_series_testing_data" or die "can not find testing file";
my @all=<MM>;
chomp $all[0];
close MM;
my $dat=from_json($all[0]);


undef $para;
$para={
	id=>"test1",
	type=>"KBaseExpression.ExpressionSeries-1.0",
	data=>$dat,
	workspace=>$ws_name
};

ok($info=$ws_client->save_object($para), "uploading testing expression series object into $ws_name");
#print Dumper($info);
=cut






#after created a workspace and uploaded expression series data object
#Now, start testing coexpression scripts
#Instantiating client object for coexpression service
#my $obj = Bio::KBase::CoExpression::Client->new("http://140.221.85.182:7063");
my $cc;
ok( $cc = Bio::KBase::CoExpression::Client->new("http://localhost:7063"), "get coexpression client from localhost");
#print Dumper($obj);



my $ws_id='coexpression_randomtest_';
my $job_id = $cc->filter_genes({ws_id => 'coexpression_randomtest_', inobj_id => 'test1', outobj_id => 'test1.filtered', p_value => "0.05", method => "anova"});

ok(ref($job_id) eq "ARRAY","const_coex_net_clust returns an array");
ok(@{$job_id} eq 2, "returns two job ids for const_coex_net_clust");

print "download the newly-generated network object from workspace and check the number of clusters in it\n";


die;

undef $para;
$para->{'id'}="test1.filtered";
$para->{'workspace'}="coexpression_randomtest_";
my $filtered_ob=$ws_client->get_object($para);

print Dumper($filtered_ob);






























die;

#delete this testing workspace
print "deleting this workspace\n";
$info=$ws_client->delete_workspace({
	workspace =>$ws_name
});






