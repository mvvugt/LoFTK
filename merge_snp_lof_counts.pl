#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use feature qw/ switch say /;
use Getopt::Long;
use Carp;
use IO::File;

my $USAGE = << "And be among her cloudy trophies hung.";

merge_snp_lof_counts: a tool to merge SNP LoF counts files from multiple studies for mega-analysis.
Copyright 2016 Brian Sebastian Cole, PhD. (Adapted for SNP-level data)

./merge_snp_lof_counts -i first_snp_lof_counts,second_snp_lof_counts -o merged_snp_lof_counts

Options:
-i  --input_files        Comma-separated list of input file names to merge. Required.
-o  --output_file        Output file to write merged data to. Required.
-u  --union              Compute union instead of intersection. Default: disabled.
-w  --wing_value         Value to supply for SNPs in the union but not the intersection when --union is enabled. Default: 0
-c  --clobber            Flag to overwrite ("clobber") the output file if it exists. Default: disabled.
-h  --help               Print this message and exit.

And be among her cloudy trophies hung.

main(); #Run program.

sub open_input_files {
  #Given a comma-separated string of input files, return opened filehandles for reading.
  #Croak if there are one or fewer file names provided.
  #Additionally validate the existence of each of the files, croak if any are not files.
  my $input_files_string = shift;
  my @input_file_names = split /,/ , $input_files_string;
  scalar @input_file_names < 2 and croak "Need more than one input file to process: received the string $input_files_string which had " , scalar @input_file_names , " files.\n";
  map { croak "Invalid file provided: $_\n" unless -f } @input_file_names;
  return map IO::File->new( $_ ) , @input_file_names;
}

sub array_count {
  #Given a reference to an array and a scalar number to count, count the number of occurrences of that scalar with the array and return the count as a scalar.
  my ( $array_reference , $to_count ) = @_;
  my $count = 0;
  map { $_ == $to_count and $count++ } @$array_reference;
  return $count;
}

sub parse_snp_lofs_line {
  #Given a reference to a string containing a raw data line from a SNP LoF file, return a key-value pair where the key is a tab-joined 
  #SNP_ID, Allele, Consequence, gene_ID, and gene_symbol, and the value is a reference to a hash containing het/hom LoF counts, 
  #total samples, and the actual genotype data as a tab-delimited string.
  my $snp_lofs_line_reference = shift;
  chomp $$snp_lofs_line_reference;
  my @fields = split /\t/ , $$snp_lofs_line_reference;
  carp "Failed to split a line into fields. Line was this: $$snp_lofs_line_reference\n" unless scalar @fields >= 9; #Need at least 9 fields.
  
  my $snp_id       = shift @fields;
  my $allele       = shift @fields;
  my $consequence  = shift @fields;
  my $gene_id      = shift @fields;
  my $gene_symbol  = shift @fields;
  my $het_freq     = shift @fields; # Skip - will be recomputed
  my $hom_freq     = shift @fields; # Skip - will be recomputed
  my $het_carriers = shift @fields; # Skip - will be recomputed
  my $hom_carriers = shift @fields; # Skip - will be recomputed
  
  my $key = join "\t" , ( $snp_id , $allele , $consequence , $gene_id , $gene_symbol );
  my $value = { heterozygous_lof_count => array_count( \@fields , 1 ) ,
                homozygous_lof_count   => array_count( \@fields , 2 ) ,
                total_samples          => scalar @fields              ,
                genotypes              => join "\t" , @fields         ,
              };
  return ( $key , $value );
}

sub get_samples_string_from_header_line {
  #Given a reference to a raw header line, chomp it, extract the samples, and return a tab-delimited string of the samples.
  my $line_ref = shift;
  chomp $$line_ref;
  my @fields = split /\t/ , $$line_ref;
  # Skip first 9 columns: SNP_ID, Allele, Consequence, gene_ID, gene_symbol, het_freq, hom_freq, het_carriers, hom_carriers
  splice @fields, 0, 9;
  return join "\t" , @fields;
}

sub read_snp_lofs_file {
  #Given an opened filehandle to a SNP LoF file, return a string of the samples and a hash keyed by SNP info
  #and valued by a reference to a hash that contains genotypes as well as frequency of het and hom LoF.
  my $filehandle = shift;
  my ( %genotypes , );
  my $header_line = $filehandle->getline();
  my $samples_string = get_samples_string_from_header_line( \$header_line );
  while ( defined( my $raw_line = $filehandle->getline() ) ) {
    my ( $key , $value ) = parse_snp_lofs_line( \$raw_line );
    $genotypes{ $key } = $value;
  }
  return ( \$samples_string , \%genotypes );
}

sub get_intersection_of_hash_keys {
  #Given an arbitrary number of references to hashes, return an array reference containing the intersection of the keys.
  my @hash_references = @_;
  croak "Need more than one hash reference to get the intersection of hash keys from.\n" unless scalar @hash_references > 1;
  my @intersection_of_hash_keys;
  my $first_hash_reference = shift @hash_references;
  KEY: for my $key ( keys %$first_hash_reference ) {
    HASH: for my $hash_reference ( @hash_references ) {
      next KEY unless exists $hash_reference->{ $key };
    }
    push @intersection_of_hash_keys , $key;
  }
  return \@intersection_of_hash_keys;
}

sub get_union_of_hash_keys {
  #Given an arbitrary number of references to hashes, return an array reference containing the union of the keys.
  my @hash_references = @_;
  croak "Need more than one hash reference to get the union of hash keys from.\n" unless scalar @hash_references > 1;
  my %union_of_hash_keys;
  for my $hash_reference ( @hash_references ) {
    for my $key ( keys %$hash_reference ) {
      $union_of_hash_keys{ $key } = 1;
    }
  }
  my @union_of_hash_keys = keys %union_of_hash_keys;
  return \@union_of_hash_keys;
}

sub get_sample_size {
  #Return the count of "\t" plus 1.
  my $sample_string_reference = shift;
  my $tab_hits = () = $$sample_string_reference =~ /\t/g;
  return $tab_hits + 1;
}

sub repeat_string {
  my ( $string , $times , $spacer ) = @_;
  my $repeated_string = $string;
  for ( 2 .. $times ) {
    $repeated_string .= $spacer . $string;
  }
  return $repeated_string;
}

sub divide {
  my ( $numerator , $denominator ) = @_;
  if ( $denominator == 0 ) {
    carp "Illegal division by zero. This shouldn't happen. Returning NA.";
    return 'NA';
  }
  else {
    if ( $numerator == 0 ) {
      return 0;
    }
    else {
      return $numerator / $denominator;
    }
  }
}

sub main {
  #Collect and validate options.
  my ( $input_files , );
  my ( $output_file , $wing_value , ) = ( 'merged_snp_lofs' , 0 , );
  my ( $union , $clobber , $help , );

  GetOptions( "input_files=s" => \$input_files ,
              "output_file=s" => \$output_file ,
              "wing_value=s"  => \$wing_value  ,
              "union"         => \$union       ,
              "clobber"       => \$clobber     ,
              "help"          => \$help        ,
            );

  $help and print $USAGE and exit;
  die "Need input files.\n$USAGE"                 unless $input_files;
  die "Need an output file to write to.\n$USAGE"  unless $output_file;
  die "Need a wing value for union mode.\n$USAGE" if ( $union and not defined $wing_value );
  die "Not clobbering $output_file.\n$USAGE"      if ( -f $output_file and not $clobber );

  #Read in all the input files.
  my @contents;
  map { push @contents , [ read_snp_lofs_file( $_ ) ] } open_input_files( $input_files );

  #Write header for the output.
  my $output = IO::File->new( $output_file , 'w' );
  my $header = join "\t" , ( "SNP_ID" , "Allele" , "Consequence" , "gene_ID" , "gene_symbol" , 
                             "heterozygous_LoF_frequency" , "homozygous_LoF_frequency" , 
                             "heterozygous_LoF_carriers" , "homozygous_LoF_carriers" );
  for ( @contents ) {
    $header .= "\t" . ${ $_->[0] };
  }
  $output->say( $header );

  #Simplify the merge operation by populating arrays of the datasets.
  my ( @list_of_genotypes_hash_refs , @list_of_samples_strings );
  for ( @contents ) {
    push @list_of_genotypes_hash_refs , $_->[1];
    push @list_of_samples_strings     , $_->[0];
  }
  my @study_sizes = map { get_sample_size( $_ ) } @list_of_samples_strings;

  if ( $union ) {
    my $union_of_hash_keys = get_union_of_hash_keys( @list_of_genotypes_hash_refs );
    for my $hash_key ( @$union_of_hash_keys ) {
      my $output_string;
      my ( $total_samples , $total_merged_het_lofs , $total_merged_hom_lofs ) = ( 0 , 0 , 0 );
      for my $study_index ( 0 .. $#contents ) {
        my ( $study , $study_size ) = ( $contents[ $study_index ][1] , $study_sizes[ $study_index ] );
        $total_samples += $study_size;
        if ( exists $study->{ $hash_key } ) {
          if ( $output_string ) {
            $output_string .= "\t" . $study->{ $hash_key }{ genotypes };
          }
          else {
            $output_string  = $study->{ $hash_key }{ genotypes };
          }
          $total_merged_het_lofs += $study->{ $hash_key }{ heterozygous_lof_count };
          $total_merged_hom_lofs += $study->{ $hash_key }{ homozygous_lof_count };
        }
        else {
          my $wing_values_string = repeat_string( $wing_value , $study_size , "\t" );
          if ( $output_string ) {
            $output_string .= "\t" . $wing_values_string;
          }
          else {
            $output_string = $wing_values_string;
          }
        }
      }
      my ( $het_frequency , $hom_frequency ) = ( divide( $total_merged_het_lofs , $total_samples ) , 
                                                 divide ( $total_merged_hom_lofs , $total_samples ) );
      $output->say( join "\t" , ( $hash_key , $het_frequency , $hom_frequency , 
                                  $total_merged_het_lofs , $total_merged_hom_lofs , $output_string ) );
    }
  }
  else {
    my $intersection_of_hash_keys = get_intersection_of_hash_keys( @list_of_genotypes_hash_refs );
    for my $hash_key ( @$intersection_of_hash_keys ) {
      my $output_string;
      my ( $total_samples , $total_merged_het_lofs , $total_merged_hom_lofs ) = ( 0 , 0 , 0 );
      for my $study_index ( 0 .. $#contents ) {
        my ( $study , $study_size ) = ( $contents[ $study_index ][1] , $study_sizes[ $study_index ] );
        $total_samples += $study_size;
        $output_string .= "\t" if $output_string;
        $output_string .= $study->{ $hash_key }{ genotypes };
        $total_merged_het_lofs += $study->{ $hash_key }{ heterozygous_lof_count };
        $total_merged_hom_lofs += $study->{ $hash_key }{ homozygous_lof_count };
      }
      my ( $het_frequency , $hom_frequency ) = ( divide( $total_merged_het_lofs , $total_samples ) , 
                                                 divide ( $total_merged_hom_lofs , $total_samples ) );
      $output->say( join "\t" , ( $hash_key , $het_frequency , $hom_frequency , 
                                  $total_merged_het_lofs , $total_merged_hom_lofs , $output_string ) );
    }
  }
}
