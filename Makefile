TOP_DIR = ../..
DEPLOY_RUNTIME = /kb/runtime
TARGET = /kb/deployment
SERVICE_SPEC = CoExpression.spec
SERVICE_NAME = CoExpression
SERVICE_PSGI_FILE = $(SERVICE_NAME).psgi
SERVICE_DIR = $(TARGET)/services/$(SERVICE_NAME)
SERVER_MODULE = lib/Bio/KBase/$(SERVICE_NAME)/Service.pm
#SERVICE = CoExpressionService
SERVICE_PORT = 7063

TPAGE = $(DEPLOY_RUNTIME)/bin/tpage
TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(DEPLOY_RUNTIME) --define kb_service_name=$(SERVICE_NAME) \
        --define kb_service_port=$(SERVICE_PORT)

#include $(TOP_DIR)/tools/Makefile.common

# to wrap scripts and deploy them to $(TARGET)/bin using tools in
# the dev_container. right now, these vars are defined in
# Makefile.common, so it's redundant here.
TOOLS_DIR = $(TOP_DIR)/tools
WRAP_PERL_TOOL = wrap_perl
WRAP_PERL_SCRIPT = bash $(TOOLS_DIR)/$(WRAP_PERL_TOOL).sh
SRC_PERL = $(wildcard scripts/*.pl)

WRAP_RSCRIPT_TOOL = wrap_rscript
WRAP_RSCRIPT_SCRIPT = bash $(TOOLS_DIR)/$(WRAP_RSCRIPT_TOOL).sh
SRC_RSCRIPT = $(wildcard scripts/*.r)

WRAP_PYTHON_TOOL = wrap_python
WRAP_PYTHON_SCRIPT = bash $(TOOLS_DIR)/$(WRAP_PYTHON_TOOL).sh
SRC_PYTHON = $(wildcard scripts/*.py)

# You can change these if you are putting your tests somewhere
# else or if you are not using the standard .t suffix
CLIENT_TESTS = $(wildcard t/client-tests/*.t)
SCRIPT_TESTS = $(wildcard t/script-tests/*.t)
SERVER_TESTS = $(wildcard t/server-tests/*.t)

# This is a very client centric view of release engineering.
# We assume our primary product for the community is the client
# libraries and command line interfaces on which specific 
# science applications can be built.
#
# A service is composed of a client and a server, each of which
# should be independently deployable. Clients are composed of
# an application programming interface and a command line
# interface. In our make targets, the deploy-service deploys
# the server, the deploy-client deploys the application
# programming interface libraries, and the deploy-scripts deploys
# the command line interface (usually scripts written in a
# scripting language but java executables also qualify), and the
# deploy target would be equivelant to deploying a service (client
# libs, scripts, and server).
#
# Because the deployment of the server side code depends on the
# specific software module being deployed, the strategy needs
# to be one that leaves this decision to the module developer.
# This is done by having the deploy target depend on the
# deploy-service target. The module developer who chooses for
# good reason not to deploy the server with the client simply
# manages this dependancy accordingly. One option is to have
# a deploy-service target that does nothing, the other is to
# remove the dependancy from the deploy target.
#
# A smiliar naming convention is used for tests. 


default:

# Test Section

test: test-client test-scripts 
	echo "running client and script tests"

# test-all is deprecated. 
# test-all: test-client test-scripts test-service
#
# What does it mean to test a client. This is a test of a client
# library. If it is a client-server module, then it should be
# run against a running server. You can say that this also tests
# the server, and I agree. You can add a test-service dependancy
# to the test-client target if it makes sense to you. This test
# example assumes there is already a tested running server.
test-client:
	# run each test
	for t in $(CLIENT_TESTS) ; do \
		if [ -f $$t ] ; then \
			$(DEPLOY_RUNTIME)/bin/perl $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

# What does it mean to test a script? A script test should test
# the command line scripts. If the script is a client in a client-
# server architecture, then there should be tests against a 
# running server. You can add a test-service dependancy to the
# test-client target. You could also add a deploy-service and
# start-server dependancy to the test-scripts target if it makes
# sense to you. Future versions of the make files for services
# will move in this direction.
test-scripts:
	# run each test
	for t in $(SCRIPT_TESTS) ; do \
		if [ -f $$t ] ; then \
			$(DEPLOY_RUNTIME)/bin/perl $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

# What does it mean to test a server. A server test should not
# rely on the client libraries or scripts in so far as you should
# not have a test-service target that depends on the test-client
# or test-scripts targets. Otherwise, a circular dependency
# graph could result.
test-server:
	# run each test
	for t in $(SERVER_TESTS) ; do \
		if [ -f $$t ] ; then \
			$(DEPLOY_RUNTIME)/bin/perl $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

# Deployment:
# We are assuming our primary product to the community are
# client side application programming interface libraries and
# command line interface (scripts). The deployment of client
# artifacts should not be dependent on deployment of a server,
# although we recommend deploying the server code with the
# client code if it is useful. We will assume it is useful
# in this target, just delete the dependancy on deploy-service
# if you don't want the server code deployed with the client
# code.
#
# When deploying the client side artifacts, deployment of the
# server is optional. For standard kbase services that implement
# the server stubs generated by the type compiler it can make a
# lot of sense to deploy the server files along with the client
# files. Ultimately, this is decision depends on the specifics
# of the server and it's related architecture. For illustrative
# purposes, we include the dependency in the deploy target as we
# prefer this when it is reasonable."
deploy: deploy-client deploy-service
deploy-all: deploy-client deploy-service

#
# Deploy client should deploy the client artifacts, mainly
# the application programming interface libraries, command
# line scripts, and associated reference documentation.
deploy-client: deploy-dir deploy-libs deploy-scripts deploy-docs

# The deploy-libs and deploy-scripts targets are used to recognize
# and delineate the client types, mainly a set of libraries that
# implement an application programming interface and a set of 
# command line scripts that provide command based execution of
# individual api functions and aggregated sets of api functions.
deploy-libs: build-libs

# Deploying scripts needs some special care. They need to run
# in a certain runtime environment. Users should not have
# to modify their user environments to run kbase scripts other
# than just sourcing a single user-env script. The creation
# of this user-env script is the responsibility of the code
# that builds all the kbase modules. In the code below, we
# run a script in the dev_container tools directory that 
# wraps perl scripts. The name of the perl wrapper script is
# kept in the WRAP_PERL_SCRIPT make variable. This script
# requires some information that is passed to it by way
# of exported environment variables in the bash script below.
#
# What does it mean to wrap a perl script? To wrap a perl
# script means that a bash script is created that sets
# all required envirnment variables and then calls the perl
# script using the perl interperter in the kbase runtime.
# For this to work, both the actual script and the newly 
# created shell script have to be deployed. When a perl
# script is wrapped, it is first copied to TARGET/plbin.
# The shell script can now be created because the necessary
# environment variables are known and the location of the
# script is known. 
deploy-scripts:
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PERL_PATH=$(TARGET)/lib bash ; \
	for src in $(SRC_PERL) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .pl`; \
		echo install $$src $$base ; \
		cp $$src $(TARGET)/plbin ; \
		$(WRAP_PERL_SCRIPT) "$(TARGET)/plbin/$$basefile" $(TARGET)/bin/$$base ; \
	done
	for src in $(SRC_PYTHON) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .py`; \
		echo install $$src $$base ; \
		cp $$src $(TARGET)/pybin ; \
		$(WRAP_PYTHON_SCRIPT) "$(TARGET)/pybin/$$basefile" $(TARGET)/bin/$$base ; \
	done
	for src in $(SRC_RSCRIPT) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .r`; \
		echo install $$src $$base ; \
		cp $$src $(TARGET)/rsbin ; \
		$(WRAP_RSCRIPT_SCRIPT) "$(TARGET)/rsbin/$$basefile" $(TARGET)/bin/$$base ; \
	done

# Deploying a server refers to the deployment of ...{TODO}
deploy-service: deploy-dir deploy-libs deploy-scripts deploy-services deploy-monit deploy-docs

# TODO: need to be updated later
deploy-libs:
	rsync -arv lib/. $(TARGET)/lib/.

deploy-dir:
	mkdir -p $(SERVICE_DIR) 
	mkdir -p $(SERVICE_DIR)/webroot

deploy-monit:
	$(TPAGE) $(TPAGE_ARGS) service/process.tt > $(TARGET)/services/$(SERVICE_NAME)/process.pid

deploy-services:
	$(TPAGE) $(TPAGE_ARGS) service/start_service.tt > $(TARGET)/services/$(SERVICE_NAME)/start_service
	chmod +x $(TARGET)/services/$(SERVICE_NAME)/start_service
	$(TPAGE) $(TPAGE_ARGS) service/stop_service.tt > $(TARGET)/services/$(SERVICE_NAME)/stop_service
	chmod +x $(TARGET)/services/$(SERVICE_NAME)/stop_service


# Deploying docs here refers to the deployment of documentation
# of the API. We'll include a description of deploying documentation
# of command line interface scripts when we have a better understanding of
# how to standardize and automate CLI documentation.
deploy-docs: build-docs
	cp docs/*.html $(TARGET)/services/$(SERVICE_NAME)/webroot/.

# The location of the Client.pm file depends on the --client param
# that is provided to the compile_typespec command. The
# compile_typespec command is called in the build-libs target.
build-docs: compile-docs
	mkdir -p docs; pod2html --infile=lib/Bio/KBase/$(SERVICE_NAME)/Client.pm --outfile=docs/$(SERVICE_NAME).html
	for src in $(SRC_PERL) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .pl`; \
		pod2html $$src > docs/$$base.html ; \
	done
	for src in $(SRC_RSCRIPT) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .r`; \
	done
# Not sure what to do with RScript

# Use this if you want to unlink the generation of the docs from
# the generation of the libs. Not recommended, but could be a
# reason for it that I'm not seeing.
# The compile-docs should depend on build-libs so that we are ensured
# of having a set of documentation that is based on the latest
# type spec.
compile-docs: build-libs

# Build libs should be dependent on the type specification and the
# type compiler. Building the libs in this way means that you don't
# need to put automatically generated code in a source code version
# control repository (ie cvs, git). It also ensures that you always
# have the most  up-to-date libs and documentation if your compile
# docs depends on the compiled libs.
build-libs:
	mkdir -p scripts; compile_typespec \
		--psgi $(SERVICE_PSGI_FILE) \
		--impl Bio::KBase::$(SERVICE_NAME)::$(SERVICE_NAME)Impl \
		--service Bio::KBase::$(SERVICE_NAME)::Service \
		--client Bio::KBase::$(SERVICE_NAME)::Client \
		--py biokbase/$(SERVICE_NAME)/Client \
		--js javascript/$(SERVICE_NAME)/Client \
		$(SERVICE_SPEC) lib

#		--scripts scripts \ # script is not working
