use strict;
use Bio::KBase::CoExpression::Client;


my $ws_id="kbasetest:home";
my $cc = Bio::KBase::CoExpression::Client->new("http://localhost:7063");
my $job_id = $cc->filter_genes({ws_id => $ws_id, inobj_id => 'kb|series.436', outobj_id => 'kb|series.436.ftd', p_value => "0.05", method => "anova", num_genes => "30"});
print "$$job_id[0] $$job_id[1]\n";

$job_id = $cc->const_coex_net_clust({ws_id => $ws_id, inobj_id => 'kb|series.436.ftd', outobj_id => 'net_by_kb|series.436', cut_off => "0.75", net_method => "simple", clust_method => "hclust", num_modules => "3"});
print "$$job_id[0] $$job_id[1]\n";

