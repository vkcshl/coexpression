use Bio::KBase::CoExpression::CoExpressionImpl;

use Bio::KBase::CoExpression::Service;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = Bio::KBase::CoExpression::CoExpressionImpl->new;
    push(@dispatch, 'CoExpression' => $obj);
}


my $server = Bio::KBase::CoExpression::Service->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");
