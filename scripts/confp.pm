#!/usr/bin/env perl

package Confp;

use strict;
use warnings;

sub new {
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self = {};

	if (@_) {
		if (!$self->read(shift))
	} else {
		$self->{PATH} = "";
	}

	$self->{PATH} = "";
	$self->{CONFIG}   = {};
	bless($self, $class);
	return $self;
}

sub read {
	my $self = shift;
	my $path = shift;

	# Open file
	# If cannot open {
	# 	return 0;
	# }
	
	# For each line {
	# 	$self->parse(line);
	# }
	
	# return 1;
}

sub parse { 
	my $self = shift;
	my $line = shift;

	# Parse line
	my $key;
	my $value;



	$self->config{$key} = $value;
}


sub get { 
	my $self = shift;
	my $key  = shift;

	return $self->config{$key};
}

1;
