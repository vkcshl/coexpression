package Bio::KBase::CoExpression::Client;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

Bio::KBase::CoExpression::Client

=head1 DESCRIPTION


Co-Expression Service APIs 

 This module provides services for plant expression data in support of the coexpression
 network and ontology driven data needs of the plant sciences community. This version of
 the modules supports retrieval of the following information:
 1. Retrieval of GEO sample ID list for given EO (environmental ontology) and/or PO (plant ontology -plant tissues/organs of interest).
 2. Retrieval of the expression values for given GEO sample ID list.  
 3. For given expression values tables, it computes co-expression clusters or network (CLI only).

It will serve queries for tissue or condition specific co-expression network for biologically interesting genes/samples. Users can search differentially expressed genes in different tissues or in numerous experimental conditions or treatments (e.g various biotic or abiotic stresses). Currently the metadata annotation is provided for a subset of gene expression experiments from the NCBI GEO microarray experiments for Arabidopsis and Poplar. The samples of these experiments are manually annotated using plant ontology (PO) [http://www.plantontology.org/] and environment ontology (EO) [http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/phenotype/environment/environment_ontology.obo]


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => Bio::KBase::CoExpression::Client::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my $token = Bio::KBase::AuthToken->new(@args);
	
	if (!$token->error_message)
	{
	    $self->{token} = $token->token;
	    $self->{client}->{token} = $token->token;
	}
        else
        {
	    #
	    # All methods in this module require authentication. In this case, if we
	    # don't have a token, we can't continue.
	    #
	    die "Authentication failed: " . $token->error_message;
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




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

Description of filter_genes: 
filter_genes provides the function to identify differentially expressed genes given an expression series/experiment. An expression series/experiment contains a list of expression samples. A expression sample is the measurement of mRNA abundance in a biological sample. The design of expression profiling usually includes replicates. The replicates allows us to differ the non-relevent expression variation and the relevent expression variation.
The replicate information is manully extracted by KBase developers. Only a part of samples has been assigned to a replicate group. For those samples without an assignment, the variation of its expression abundance is used directly.
filter_genes now has two methods to identify differentially expressed genes: ANOVA and lor(from limma r package). The output of this function is a list of genes

=back

=cut

sub filter_genes
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function filter_genes (received $n, expecting 1)");
    }
    {
	my($args) = @args;

	my @_bad_arguments;
        (ref($args) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"args\" (value was \"$args\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to filter_genes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'filter_genes');
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "CoExpression.filter_genes",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'filter_genes',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method filter_genes",
					    status_line => $self->{client}->status_line,
					    method_name => 'filter_genes',
				       );
    }
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

Description of const_coex_net_clust
const_coex_net_clust provides the function to build coexpression network and identify the functional modules among it.
A functional module is a network cluster with enrichment of certain biological function. const_coex_net_clust first construct coexpression network. Then, it identifys the clusters among the network. Finally, it identifys the GeneOntology enrichment for the genes in each cluster.

=back

=cut

sub const_coex_net_clust
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function const_coex_net_clust (received $n, expecting 1)");
    }
    {
	my($args) = @args;

	my @_bad_arguments;
        (ref($args) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"args\" (value was \"$args\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to const_coex_net_clust:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'const_coex_net_clust');
	}
    }

    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
	method => "CoExpression.const_coex_net_clust",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'const_coex_net_clust',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method const_coex_net_clust",
					    status_line => $self->{client}->status_line,
					    method_name => 'const_coex_net_clust',
				       );
    }
}



sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "CoExpression.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'const_coex_net_clust',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method const_coex_net_clust",
            status_line => $self->{client}->status_line,
            method_name => 'const_coex_net_clust',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for Bio::KBase::CoExpression::Client\n";
    }
    if ($sMajor == 0) {
        warn "Bio::KBase::CoExpression::Client version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 FilterGenesParams

=over 4



=item Description

num_gene is user for specify how many differentially expressed genes are needed


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

num_modules is used to define the number of modules need to be identified from the network


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

package Bio::KBase::CoExpression::Client::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
