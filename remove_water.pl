#!/usr/bin/perl -w
use strict;
use File::Slurp;

# Declaring file names for the out files
my $output = 'solv.gro';
my $output1 = 'solv1.gro';
my $output2 = 'sol_removed.gro';

# Opening file handlers for output files
open my $outfile,'>', $output or die "Cant write to $output: $!"; 
open my $outfile1,'>',$output1 or die "Cant write to $output: $!"; 
open my $outfile2,'>',$output2 or die "Cant write to $output: $!";
# Reading all the lines from the file and storing into an array
my @array = read_file('solvate.gro');

# Avoiding all the lines that contains SOL atoms and printing the rest to an outfile named solv.gro. This stores the number of atoms line, all the lines of lipid atoms, and box di
mension.
for my $line(@array){
	next if($line =~ m/^.*SOL.*$/);
	print $outfile $line; 
}

# All the lines that contain SOL is searched using regular expression to identify the z coordinates that is within the intial water box. If the water atom line has z coordinate less than 15.8 nm or greater than 5.4 nm, it is pushed to a hash. All other water atom lines are printed to another file named sol_removed.gro.

my %count;
for my $line(@array){
	if($line =~ m/^\s*(\d+)SOL\s+[HW|OW].*\s+-?\d+\.\d+\s+-?\d+\.\d+\s+(-?\d+\.\d+)$/){ 
		if(($2 <= 15.8) &&($2 >= 5.4)){
			push(@{$count{$1}}, $line); 
		}
		else{
		print $outfile2 $line;
		}

	 }
}

# From the Hash structure that stored all the water atom line within the set of z planes, water molecules that have three atoms is printed to an outfile named solv1.gro, the rest is printed to an outfile names sol_removed.gro. This make sure that there are no broken water molecules for next steps.
my $count = 0;
foreach my $key(keys %count){ 
	if(scalar(@{$count{$key}})==3){
		print $outfile1 @{$count{$key}};
		$count++; 
	}
	elsif(scalar(@{$count{$key}}) < 3){ 
		print $outfile2 @{$count{$key}};
	} 
}
# Printing the total number of water molecule atoms that is kept. i.e the number of lines in solv1.gro
print $count;