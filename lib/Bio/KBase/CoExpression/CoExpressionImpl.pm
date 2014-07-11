package Bio::KBase::CoExpression::CoExpressionImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

CoExpression

=head1 DESCRIPTION

Co-Expression Service APIs 

 This module provides services for plant expression data in support of the coexpression
 network and ontology driven data needs of the plant sciences community. This version of
 the modules supports retrieval of the following information:
 1. Retrieval of GEO sample ID list for given EO (environmental ontology) and/or PO (plant ontology -plant tissues/organs of interest).
 2. Retrieval of the expression values for given GEO sample ID list.  
 3. For given expression values tables, it computes co-expression clusters or network (CLI only).

It will serve queries for tissue and condition specific gene expression co-expression network for biologically interesting genes/samples. Users can search differentially expressed genes in different tissues or in numerous experimental conditions or treatments (e.g various biotic or abiotic stresses). Currently the metadata annotation is provided for a subset of gene expression experiments from the NCBI GEO microarray experiments for Arabidopsis and Poplar. The samples of these experiments are manually annotated using plant ontology (PO) [http://www.plantontology.org/] and environment ontology (EO) [http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/phenotype/environment/environment_ontology.obo]

=cut

#BEGIN_HEADER
use Bio::KBase::Workflow::KBW;
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
    my %params;
    my @list = qw(ujs_url awe_url shock_url ws_url);
    if ((my $e = $ENV{KB_DEPLOYMENT_CONFIG}) && -e $ENV{KB_DEPLOYMENT_CONFIG}) {  
      my $service = $ENV{KB_SERVICE_NAME};
      if (defined($service)) {
        my $c = Config::Simple->new();
        $c->read($e);
        for my $p (@list) {
          my $v = $c->param("$service.$p");
          if ($v) {
            $params{$p} = $v;
          }
        }
      }
    }
 
    # set default values for testing
    $params{'ujs_url'} = 'http://localhost:7083' if ! defined $params{'ujs_url'};
    $params{'awe_url'} = 'http://localhost:7080' if ! defined $params{'awe_url'};
    $params{'shock_url'} = 'http://kbase.us/services/shock-api' if ! defined $params{'shock_url'};
    $params{'ws_url'} = 'https://kbase.us/services/ws/' if ! defined $params{'ws_url'};
    $self->{_config} = \%params;
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 filter_genes

  $job_id = $obj->filter_genes($args)

=over 4

=item Parameter and return types

=begin html

<pre>
$args is a FilterGenesParams
$job_id is a reference to a list where each element is a string
FilterGenesParams is a reference to a hash where the following keys are defined:
	ws_id has a value which is a string
	inobj_id has a value which is a string
	outobj_id has a value which is a string
	p_value has a value which is a string
	method has a value which is a string
	num_genes has a value which is a string

</pre>

=end html

=begin text

$args is a FilterGenesParams
$job_id is a reference to a list where each element is a string
FilterGenesParams is a reference to a hash where the following keys are defined:
	ws_id has a value which is a string
	inobj_id has a value which is a string
	outobj_id has a value which is a string
	p_value has a value which is a string
	method has a value which is a string
	num_genes has a value which is a string


=end text



=item Description



=back

=cut

sub filter_genes
{
    my $self = shift;
    my($args) = @_;

    my @_bad_arguments;
    (ref($args) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"args\" (value was \"$args\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to filter_genes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'filter_genes');
    }

    my $ctx = $Bio::KBase::CoExpression::Service::CallContext;
    my($job_id);
    #BEGIN filter_genes
    $job_id = Bio::KBase::Workflow::KBW::run_async($self, $ctx, $args);
    #END filter_genes
    my @_bad_returns;
    (ref($job_id) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"job_id\" (value was \"$job_id\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to filter_genes:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'filter_genes');
    }
    return($job_id);
}




=head2 const_coex_net_clust

  $job_id = $obj->const_coex_net_clust($args)

=over 4

=item Parameter and return types

=begin html

<pre>
$args is a ConstCoexNetClustParams
$job_id is a reference to a list where each element is a string
ConstCoexNetClustParams is a reference to a hash where the following keys are defined:
	ws_id has a value which is a string
	inobj_id has a value which is a string
	outobj_id has a value which is a string
	cut_off has a value which is a string
	net_method has a value which is a string
	clust_method has a value which is a string
	num_modules has a value which is a string

</pre>

=end html

=begin text

$args is a ConstCoexNetClustParams
$job_id is a reference to a list where each element is a string
ConstCoexNetClustParams is a reference to a hash where the following keys are defined:
	ws_id has a value which is a string
	inobj_id has a value which is a string
	outobj_id has a value which is a string
	cut_off has a value which is a string
	net_method has a value which is a string
	clust_method has a value which is a string
	num_modules has a value which is a string


=end text



=item Description



=back

=cut

sub const_coex_net_clust
{
    my $self = shift;
    my($args) = @_;

    my @_bad_arguments;
    (ref($args) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"args\" (value was \"$args\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to const_coex_net_clust:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'const_coex_net_clust');
    }

    my $ctx = $Bio::KBase::CoExpression::Service::CallContext;
    my($job_id);
    #BEGIN const_coex_net_clust
    $job_id = Bio::KBase::Workflow::KBW::run_async($self, $ctx, $args);
    #END const_coex_net_clust
    my @_bad_returns;
    (ref($job_id) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"job_id\" (value was \"$job_id\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to const_coex_net_clust:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'const_coex_net_clust');
    }
    return($job_id);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 FilterGenesParams

=over 4



=item Description

series object id


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
ws_id has a value which is a string
inobj_id has a value which is a string
outobj_id has a value which is a string
p_value has a value which is a string
method has a value which is a string
num_genes has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
ws_id has a value which is a string
inobj_id has a value which is a string
outobj_id has a value which is a string
p_value has a value which is a string
method has a value which is a string
num_genes has a value which is a string


=end text

=back



=head2 ConstCoexNetClustParams

=over 4



=item Description

series object id


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
ws_id has a value which is a string
inobj_id has a value which is a string
outobj_id has a value which is a string
cut_off has a value which is a string
net_method has a value which is a string
clust_method has a value which is a string
num_modules has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
ws_id has a value which is a string
inobj_id has a value which is a string
outobj_id has a value which is a string
cut_off has a value which is a string
net_method has a value which is a string
clust_method has a value which is a string
num_modules has a value which is a string


=end text

=back



=cut

1;
